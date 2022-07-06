#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/exporter.hpp>

int main(int argc, char **argv)
{
    using namespace Feel;
    Environment env(_argc = argc, _argv = argv,
                    _desc = feel_options(),
                    _about = about(_name = "laplacian.e",
                                   _author = "Feel++ Consortium",
                                   _email = "feelpp-devel@feelpp.org"));

    auto mesh = unitSquare();
    auto Vh = Pch<1>(mesh);
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1(_test = Vh);
    l = integrate(_range = elements(mesh),
                  _expr = id(v));

    auto a = form2(_trial = Vh, _test = Vh);
    a = integrate(_range = elements(mesh),
                  _expr = gradt(u) * trans(grad(v)));
    a += on(_range = boundaryfaces(mesh), _rhs = l, _element = u,
            _expr = constant(0.));
    a.solve(_rhs = l, _solution = u);

    auto e = exporter(_mesh = mesh, _name = "laplacian");
    e->add("u", u);
    e->save();
    return 0;
}