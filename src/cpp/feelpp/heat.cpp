#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/utility.hpp>
#include <feelpp/heat.hpp>

int main(int argc, char** argv) {
    using namespace Feel;
    try{
        
        Environment env(_argc = argc, _argv = argv,
                        _desc = makeOptions(),
                        _about = about(_name = "heat",
                                    _author = "Feel++ Consortium",
                                    _email = "feelpp@cemosis.fr"));
        auto jsonfile = removeComments(readFromFile(Environment::expand(soption("specs"))));
        std::istringstream istr(jsonfile);
        json specs = json::parse(istr);
        
        // number of time interval for parareal
        int time_partitions = specs["/Time/partitions"_json_pointer].get<int>();

        // split the main worldcomm into P+1 subworldcomms
        auto [color, w, wglob] = Environment::worldCommPtr()->split(time_partitions+1);

        // coarse time step
        ///double dt = specs["/Time/dt"_json_pointer].get<double>();
        //double t0 = specs["/Time/t0"_json_pointer].get<double>();
        double T = specs["/Time/final_time"_json_pointer].get<double>();
        double t0 = 0.;
        // on global rank 0, we have the coarse integrator
        if ( w->globalRank() == 0 )
        {
            double t0_coarse = t0;
            double dt_coarse = T / time_partitions;
            double T_coarse = T;
            // coarse heat mesh: unfortunately we duplicate on the P processor the coarse integrator
            // need to fix that
            Heat<2, 1> heatcoarse( fmt::format("heat-coarse-{}",w->globalRank()), specs, w, dt_coarse, t0_coarse, T_coarse );
            
            bool work = true;
            while ( work )
            {
                for (double t = dt_coarse; t < T_coarse; t += dt_coarse) {
                    std::cout << "====================================" << std::endl;
                    std::cout << fmt::format("t = {}", t) << std::endl;
                    // execute the time step: update the right hand side and solve the system
                    heatcoarse.run( t, heat.solution() );

                    // save solution at current time
                    heatcoarse.postProcess();

                    // non blocking async comm to receive fine integrator communication
                    // 
                }
                // update work flag to know if we stop or continue
                // work = normL2(_range=elements(mesh),_expr=idv(heatcoarse.solution())-idv(heatcoarse.oldSolution()) < 1e-6;
                work=false;
            }
        }
        else // we are on the other P processors with the fine integrators
        {
            double dt_fine = specs["/Time/dt_fine"_json_pointer].get<double>();
            double t0_fine = t0 + (w->globalRank() - 1) * (T - t0) / (time_partitions + 1);
            double T_fine = t0 + (w->globalRank()) * (T - t0)/(time_partitions+1);

            // fin heat for each time subdomain, local to each subdomain
            Heat<2,1> heatfine(fmt::format("heat-fine-{}",w->globalRank()), specs, w, dt_fine, t0_fine, T_fine );
            
            bool work = true;
            while( work )
            {
                for (double t = t0_fine; t < T_fine; t += dt_fine) {
                    std::cout << "====================================" << std::endl;
                    std::cout << fmt::format("Time interval {}, t = {}", w->globalRank()-1, t) << std::endl;
                    // execute the time step: update the right hand side and solve the system
                    heatfine.run( t, heat.solution() );

                    // save solution at current time
                    heatfine.postProcess();
                }   
                // send fine solution to coarse integrator
                // communication from fine to coarse integrators
                // non-blocking communication

                // update work flag to know whether we have to continue or not
                // communication from coarse to fine integrators
                // blocking communication
                
                work = false;
            }
        }    

        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
}