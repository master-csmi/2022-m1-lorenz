#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feel.hpp>

int main(int argc, char **argv){

    using namespace Feel;

    Environment env(_argc = argc, _argv = argv,
                    _desc = feel_options(),
                    _about = about(_name = "laplacian.e",
                                   _author = "Feel++ Consortium",
                                   _email = "feelpp-devel@feelpp.org"));
    
    const unsigned short int FEELPP_DIM = 2; //modif

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);


    auto f = expr( soption(_name="functions.f"), "f" ); 
    auto r = expr( soption(_name="functions.r"), "r" ); // Robin right hand side expression 
    auto h = expr( soption(_name="functions.h"), "h" ); // Neumann expression  
    auto g = expr( soption(_name="functions.g"), "g" ); 
    

    auto Vh = Pch<1>( mesh ); 
    auto u = Vh->element("u"); 
    auto v = Vh->element( g, "g" );


    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh), _expr=f*id(v));
    l += integrate(_range=markedfaces(mesh,"Robin"), _expr=r*id(v));
    l += integrate(_range=markedfaces(mesh,"Neumann"), _expr=h*id(v));


    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh), _expr=inner(gradt(u),grad(v)) );
    a += integrate(_range=markedfaces(mesh,"Robin"), _expr=idt(u)*id(v));
    a += on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );


    //! if no markers Robin Neumann or Dirichlet are present in the mesh then
    //! impose Dirichlet boundary conditions over the entire boundary
    if ( !mesh->hasAnyMarker({"Robin", "Neumann","Dirichlet"}) )
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );

    
    a.solve(_rhs = l, _solution = u);

    auto e = exporter(_mesh = mesh, _name = "laplacian");
    e->add("uh", u);
    e->save();

    return 0;
};