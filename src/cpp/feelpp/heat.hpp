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

    using mesh_t = Mesh<Simplex<Dim>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    using space_t = Pch_type<mesh_t, Order>;
    using space_ptr_t = Pch_ptrtype<mesh_t,Order>;
    using element_t = typename space_t::element_type;
    using a_t = form2_t<space_t,space_t>;
    using l_t = form1_t<space_t>;
    using exporter_ptr_t = std::shared_ptr<Exporter<mesh_t>>;

    Heat(std::string const &prefix, nl::json const &specs, mesh_ptr_t mesh, double dt, double t0, double T)
        : prefix_(prefix),
          specs_(specs),
          mesh_(mesh),
          Xh_(Pch<Order>(mesh_)),
          v_(Xh_->element()),
          c_(Xh_->element()),
          e_(nullptr),
          a_(form2(_test = Xh_, _trial = Xh_)),
          l_(form1(_test = Xh_)),
          lt_(form1(_test = Xh_)),
          dt_(dt),
          t_(t0),
          T_(T)

    {
        // initial condition (TODO use the json specs)
        v_.on( _range = elements(mesh_), _expr = cst(0.) );

        // TODO : need to add source term here using json specs if there is one

        // assemble the matrix
        for (auto [key, material] : specs["/Models/heat/materials"_json_pointer].items())
        {
            LOG(INFO) << fmt::format("material {}", material);
            std::string mat_k = fmt::format("/Materials/{}/k", material.get<std::string>());
            auto k = specs[nl::json::json_pointer(mat_k)].get<std::string>();
            std::string mat_rho = fmt::format("/Materials/{}/rho", material.get<std::string>());
            auto rho = specs[nl::json::json_pointer(mat_rho)].get<std::string>();
            std::string mat_cp = fmt::format("/Materials/{}/Cp", material.get<std::string>());
            auto cp = specs[nl::json::json_pointer(mat_cp)].get<std::string>();
            a_ += integrate(_range = markedelements(mesh_, material.get<std::string>()), _expr = expr(rho) * expr(cp) * idt(v_) * id(v_)/dt_);
            a_ += integrate(_range = markedelements(mesh_, material.get<std::string>()), _expr = expr(k) * gradt(v_) * trans(grad(v_)));

            
        }

        // BC Neumann
        for (auto &[bc, value] : specs["/BoundaryConditions/heat/flux"_json_pointer].items())
        {
            LOG(INFO) << fmt::format("flux {}: {}", bc, value.dump());
            auto flux = value["expr"].get<std::string>();
            l_ += integrate(_range = markedfaces(mesh_, bc), _expr = expr(flux) * id(v_));
        }

        // BC Robin
        for (auto &[bc, value] : specs["/BoundaryConditions/heat/convective_heat_flux"_json_pointer].items())
        {
            LOG(INFO) << fmt::format("convective_heat_flux {}: {}", bc, value.dump());
            auto h = value["h"].get<std::string>();
            auto Text = value["Text"].get<std::string>();
            a_ += integrate(_range = markedfaces(mesh_, bc), _expr = expr(h) * id(v_) * idt(v_));
            l_ += integrate(_range = markedfaces(mesh_, bc), _expr = expr(h) * expr(Text) * id(v_));
        }
    }
    /**
     * @brief get a new element of Xh
     * 
     * @return element_t 
     */
    element_t newElement() const
    {
        return Xh_->element();
    }
    /**
     * @brief get the solution at the current time
     * 
     * @return element_t const& the solution
     */
    element_t const& solution() const
    {
        return v_;
    }

    /**
     * @brief get the solution at the current time (read-write)
     *
     * @return element_t const& the solution
     */
    element_t & solution() 
    {
        return v_;
    }
    /**
     * @brief get the correction at the current time
     *
     * @return element_t const& the correction
     */
    element_t const &correction() const
    {
        return c_;
    }

    /**
     * @brief get the correction at the current time (read-write)
     *
     * @return element_t const& the correction
     */
    element_t &correction()
    {
        return c_;
    }

    /**
     * @brief set the initial condition
     * 
     */
    void initTimeStep()
    {
        // initial condition: TODO change according to json specs
        v_.zero();
        vold_.zero();
    }

    /**
     * @brief run the time step [t0, t0+dt_]
     * 
     * @return int 
     */
    void run( double t0, element_t const& w )
    {
        t_ = t0;
        // zero out right hand side
        lt_.zero();
        lt_ += l_;
        // update right hande side terms from the time derivative
        for (auto [key, material] : specs_["/Models/heat/materials"_json_pointer].items())
        {
            LOG(INFO) << fmt::format("material {}", material);
            std::string mat_k = fmt::format("/Materials/{}/k", material.get<std::string>());
            auto k = specs_[nl::json::json_pointer(mat_k)].get<std::string>();
            std::string mat_rho = fmt::format("/Materials/{}/rho", material.get<std::string>());
            auto rho = specs_[nl::json::json_pointer(mat_rho)].get<std::string>();
            std::string mat_cp = fmt::format("/Materials/{}/Cp", material.get<std::string>());
            auto cp = specs_[nl::json::json_pointer(mat_cp)].get<std::string>();
            lt_ += integrate(_range = markedelements(mesh_, material.get<std::string>()), _expr = expr(rho) * expr(cp) * idv(w) * id(v_)/dt_);
        }
        
        // Solve
        a_.solve(_rhs = lt_, _solution = v_);
    }
    /**
     * @brief save solution at time @p t and @p iteration 
     * 
     * @param e the element to save
     * @param t time index
     * @param iteration iteration number
     * @paaram type type of solution to save (solution, correction, ...)
     */
    void save( element_t const& e, double t, int iteration, std::string type="solution" )
    {
        e.save(_path=fs::current_path().string(),_name=fmt::format("{}-{}-{:.2f}-{}",prefix_,type,t,iteration));
    }
    /**
     * @brief save solution at time @p t and @p iteration
     *
     * @param e the element to load
     * @param t time index
     * @param iteration iteration number
     * @paaram type type of solution to load (solution, correction, ...)
     */
    element_t load( double t, int iteration, std::string type = "solution" )
    {
        element_t v = Xh_->element();
        v.load(_path = fs::current_path().string(), _name = fmt::format("{}-{}-{:.2f}-{}", prefix_, type,t, iteration));
        return v;
    }
    /**
     * @brief get the exporter
     * 
     * @return exporter_ptr_t& 
     */
    void resetExporter( std::string const& prefix ) 
    {
        e_=exporter(_mesh = mesh_, _name = prefix);
    }
    /**
     * @brief postprocess result at current time
     * 
     */
    void postProcess()
    {
#if 0
        // compute outputs
        auto m = mean(_range = elements(mesh_), _expr = idv(v_));
        auto m_root = mean(_range = markedfaces(mesh_, "Gamma_root"), _expr = idv(v_));
        std::cout << fmt::format("- mean value: {}", m) << std::endl;
        std::cout << fmt::format("-  min value: {}", v_.min()) << std::endl;
        std::cout << fmt::format("-  max value: {}", v_.max()) << std::endl;
        std::cout << fmt::format("-  max deviation: {}", v_.max() - v_.min()) << std::endl;
        std::cout << fmt::format("-  mean root: {}", m_root) << std::endl;
#endif        
        // Export
        e_->step(t_)->setMesh(mesh_);
        e_->step(t_)->addRegions();
        e_->step(t_)->add("T", v_);
        e_->step(t_)->add("correction", c_);
        e_->save();

    }
private:
    std::string prefix_;
    nl::json specs_;
    mesh_ptr_t mesh_;
    space_ptr_t Xh_;
    element_t v_,vold_,c_;
    exporter_ptr_t e_;
    a_t a_;
    l_t l_, lt_;
    // time step
    double dt_;
    // current time
    double t_; 
    // final time
    double T_;
};

} // namespace Feel