#include <feel/feel.hpp>

int main(int argc, char**argv)
{
    using namespace Feel;
    using Feel::cout;
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ;

    Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="laplacian",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    auto thechecker = checker( _name= "L1/H1 convergence", 
                            _solution_key="toto"
                            );

    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");

    tic();
    auto Vh = Pch<1>( mesh ); 
    auto u = Vh->element("u"); 
    auto mu = expr(soption(_name="functions.mu")); // diffusion term 
    auto f = expr( soption(_name="functions.f"), "f" ); 
    auto r_1 = expr( soption(_name="functions.a"), "a" ); // Robin left hand side expression 
    auto r_2 = expr( soption(_name="functions.b"), "b" ); // Robin right hand side expression 
    auto n = expr( soption(_name="functions.c"), "c" ); // Neumann expression 
    auto solution = expr( thechecker.solution(), "solution" ); 
    auto g = thechecker.check()?solution:expr( soption(_name="functions.g"), "g" ); 
    auto v = Vh->element( g, "g" ); 
    toc("Vh");

    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=f*id(v));
    l+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_2*id(v));
    l+=integrate(_range=markedfaces(mesh,"Neumann"), _expr=n*id(v));
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    tic();
    a = integrate(_range=elements(mesh),
                  _expr=mu*inner(gradt(u),grad(v)) );
    toc("a");
    a+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_1*idt(u)*id(v));
    a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
    //! if no markers Robin Neumann or Dirichlet are present in the mesh then
    //! impose Dirichlet boundary conditions over the entire boundary
    if ( !mesh->hasAnyMarker({"Robin", "Neumann","Dirichlet"}) )
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
    toc("a");

    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=u);
    toc("a.solve");

    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", u );
    if ( thechecker.check() )
    {
        v.on(_range=elements(mesh), _expr=solution );
        e->add( "solution", v );
    }
    e->save();
    toc("Exporter");

    // compute l2 and h1 norm of u-u_h where u=solution
    auto norms = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            double l2 = normL2(_range=elements(mesh), _expr=idv(u)-expr(solution) );
            toc("L2 error norm");
            tic();
            double h1 = normH1(_range=elements(mesh), _expr=idv(u)-expr(solution), _grad_expr=gradv(u)-grad<2>(expr(solution)) );
            toc("H1 error norm");
            return { { "L2", l2 }, {  "H1", h1 } };
        };

    int status = thechecker.runOnce( norms, rate::hp( mesh->hMax(), Vh->fe()->order() ) );

    // exit status = 0 means no error
    return !status;

}