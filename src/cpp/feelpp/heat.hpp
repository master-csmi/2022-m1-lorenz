/**
 * @file rht.cpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief Radiative heat transfer
 * @version 0.1
 * @date 2022-07-22
 *
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options("rht options");
    options.add_options()
        // mesh parameters
        ("specs", Feel::po::value<std::string>(),
         "json spec file for rht");

     return options.add(Feel::feel_options());
}

template<int Dim, int Order>
class Heat
{
public:

using mesh_t = MeshType;
using mesh_ptr_t = std::shared_ptr<mesh_t>;
using space_ptr_t = Pch_ptr_t<mesh_t,Order>;

Heat(nl::json const& specs)
{

    auto mesh_ = loadMesh(_mesh = new Mesh<Simplex<2>>, _filename = specs["/Meshes/heat/Import/filename"_json_pointer].get<std::string>());
    Xh_ = Pch<2>(mesh_);
    e_ = exporter(_mesh = mesh_);
    auto v_ = Xh_->element();
    auto a = form2(_test = Xh_, _trial = Xh_);
    auto l = form1(_test = Xh_);

    for (auto [key, material] : specs["/Models/heat/materials"_json_pointer].items())
    {
        LOG(INFO) << fmt::format("material {}", material);
        std::string mat = fmt::format("/Materials/{}/k", material.get<std::string>());
        auto k = specs[nl::json::json_pointer(mat)].get<std::string>();
        a += integrate(_range = markedelements(mesh_, material.get<std::string>()), _expr = expr(k) * gradt(v_) * trans(grad(v_)));
    }

    // BC Neumann
    for (auto &[bc, value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items())
    {
        LOG(INFO) << fmt::format("flux {}: {}", bc, value.dump());
        auto flux = value["expr"].get<std::string>();
        l += integrate(_range = markedfaces(mesh_, bc), _expr = expr(flux) * id(v_));
    }

    // BC Robin
    for (auto &[bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items())
    {
        LOG(INFO) << fmt::format("convective_heat_flux {}: {}", bc, value.dump());
        auto h = value["h"].get<std::string>();
        auto Text = value["Text"].get<std::string>();
        a += integrate(_range = markedfaces(mesh_, bc), _expr = expr(h) * id(v) * idt(v_));
        l += integrate(_range = markedfaces(mesh_, bc), _expr = expr(h) * expr(Text) * id(v_));
    }
}
int run()
{
    // Solve
    a.solve(_rhs = l, _solution = v);
}
void postProcess()
{
    // compute outputs
    auto m = mean(_range = elements(mesh_), _expr = idv(v_));
    auto m_root = mean(_range = markedfaces(mesh_, "Gamma_root"), _expr = idv(v_));
    std::cout << fmt::format("- mean value: {}", m) << std::endl;
    std::cout << fmt::format("-  min value: {}", v_.min()) << std::endl;
    std::cout << fmt::format("-  max value: {}", v_.max()) << std::endl;
    std::cout << fmt::format("-  max deviation: {}", v_.max() - v_.min()) << std::endl;
    std::cout << fmt::format("-  mean root: {}", m_root) << std::endl;
    // Export
    auto 
    e->addRegions();
    e->add("T", v_);
    e->save();

}
private:
    mesh_ptr_t mesh_;
    space_ptr_t Xh_;
    element_t v_;
    Exporter<mesh_t> e_;
};

} // namespace Feel