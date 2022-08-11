
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/utility.hpp>
#include <feelpp/heat.hpp>
#include <boost/mpi/intercommunicator.hpp>
#include <iostream>

int main(int argc, char** argv) {
    using namespace Feel;
    try{
        // mpi::communicator world;

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
        LOG(INFO) << fmt::format("time partitions: {}, n time communicators: {}", time_partitions, time_partitions+1);

        // split the main worldcomm into P+1 subworldcomms
        auto [color, w, wglob] = Environment::worldCommPtr()->split(time_partitions+1);
        LOG(INFO) << fmt::format("color: {}, wglob rank: {}, w global rank: {}, w local rank: {}", color, wglob->globalRank(), w->globalRank(), w->localRank());

        // coarse time step
        ///double dt = specs["/Time/dt"_json_pointer].get<double>();
        //double t0 = specs["/Time/t0"_json_pointer].get<double>();
        // double T = specs["/Time/final_time"_json_pointer].get<double>();
        // double t0 = 0.;
        
        #if 1
        int nb_proc = wglob->globalComm().size();
        int nb_grp = time_partitions+1;
        int nb_dom = nb_proc/nb_grp;

        if(wglob->globalRank()==0){
            int position;
            
            // size = ( wglob->globalComm().size()/(time_partitions+1) line ; time_partitions+1 col )
            std::vector<std::vector<int>> table(nb_dom, std::vector <int>(nb_grp));

            for(int p=1; p < nb_dom*nb_grp; ++p){
                wglob->globalComm().recv( p, p, position );
                table[position/nb_grp][position%nb_grp]=p;
            }
            
            std::cout << "TABLE : " << std::endl;
            for (int i=0;i<nb_dom;++i){
                for (int j=0;j<nb_grp;++j)
                    std::cout << table[i][j] << " ";
                std::cout << std::endl;
            }

            int proc;
            std::vector<int> vec_send(2);
            for(int c=0; c<nb_grp; c++){
                for(int d=0; d<nb_dom; d++){
                    proc = table[d][c];
                    if(proc != 0){
                        vec_send = {c, d};
                        wglob->globalComm().send(proc, proc, vec_send);
                    }
                }
            }
        }
        else{
            wglob->globalComm().send( 0, wglob->globalRank(), w->globalRank()*nb_grp + color);
            std::vector<int> vec_recv(2);
            wglob->globalComm().recv( 0, wglob->globalRank(), vec_recv);
            std::cout << "real : (" << color << "," << w->globalRank() << ") ; recv : (" << vec_recv[0] << "," << vec_recv[1] << ")" << std::endl;
        }
        #endif

        if(wglob->globalRank()==0){
            std::cout << "INTERCOMMUNICATOR : " << std::endl;
        }

        mpi::communicator world = wglob->globalComm();
        mpi::communicator myComm = w->localComm();

        if(color==0){
            // std::cout << world.size() << std::endl;
            mpi::intercommunicator myFirstComm(myComm,0,world,1);
            std::cout << myFirstComm.local_size() << std::endl;
            mpi::intercommunicator mySecondComm(myComm,0,world,2);

            // std::cout << wglob->globalRank() << " : " << myFirstComm.local_rank() << std::endl;
        }


        #if 0
        // on global rank 0, we have the coarse integrator
        if ( color == 0 )
        {
            double t0_coarse = t0;
            double dt_coarse = T / time_partitions;
            double T_coarse = T;
            auto mesh = loadMesh(_mesh = new Mesh<Simplex<2>>(), _worldcomm=w,
                                 _filename = specs["/Meshes/heat/Import/filename"_json_pointer].get<std::string>());
            wglob->barrier();
            // coarse heat mesh: unfortunately we duplicate on the P processor the coarse integrator
            // need to fix that
            Heat<2, 1> heatcoarse( fmt::format("heat-coarse-{}-{}",color,0), specs, mesh, dt_coarse, t0_coarse, T_coarse );
            
            int iteration = 1;
            bool work = true;
            while ( work )
            {
                heatcoarse.resetExporter(fmt::format("heat-coarse-{}-{}", color, iteration) );
                for (double t = dt_coarse; t < T_coarse+dt_coarse; t += dt_coarse) {
                    if ( w->isMasterRank() )
                    {
                        std::cout << "====================================" << std::endl;
                        std::cout << fmt::format("Coarse integrator t = {}", t) << std::endl;
                    }

                    // execute the time step: update the right hand side and solve the system
                    // compute G(U^k_{j-1}), heatcoarse.solution() == U^k_{j-1}
                    heatcoarse.run(t, heatcoarse.solution());
                    // now heatcoarse.solution() == G(U^k_{j-1})

                    // non blocking async comm to receive fine integrator communication
                    // send U^k_j =  G(U^k_{j-1})+(F(U^{k-1}_{j-1})-G(U^{k-1}_{j-1}))
                    // U^k_j =  G(U^k_{j-1}) + correction[j-1]
                    // save solution at current time
                    heatcoarse.postProcess();

                    
                    
                }
                // get in sync with the fine integrators for each coarse time step
                // we have a collection of size P of solution for each coarse time step 
                bool done = true;
                for (double t = dt_coarse; t < T_coarse + dt_coarse; t += dt_coarse)
                {
                    // we are in coarse iteration j

                    // Receive F(U^{k-1}_{j-1}) from the fine integrator
                    // ...
                    //  Compute correction
                    // correction[j-1] = F(U^{k-1}_{j-1}) - G(U^{k-1}_{j-1})

                    // update work flag to know if we stop or continue
                    // done = done && ( normL2(_range=elements(mesh),_expr=idv(U^{k}_{j})-idv(U^{k-1}_{j}) < 1e-6 );
                }
                
                work = !done;
                // broadcast work flag to all processors
                // mpi::broadcast(wglob->globalComm(), 0, work);
                ++iteration;
            }
            LOG(INFO) << fmt::format("== coarse integrator is finished ==================================") << std::endl;
        }
        else // we are on the other P processors with the fine integrators
        {
            int fine_time_interval = color-1;
            double dt_fine = specs["/Time/dt_fine"_json_pointer].get<double>();
            double t0_fine = t0 + (color-1) * (T - t0) / (time_partitions);
            double T_fine = t0 + (color) * (T - t0)/(time_partitions);

            // barrier to make sure that the mesh is created by global rank 0 process
            wglob->barrier();
            auto mesh = loadMesh(_mesh = new Mesh<Simplex<2>>(), _worldcomm = w,
                                 _filename = specs["/Meshes/heat/Import/filename"_json_pointer].get<std::string>());
            // fin heat for each time subdomain, local to each subdomain
            Heat<2,1> heatfine(fmt::format("heat-fine-{}-{}",color,0), specs, mesh, dt_fine, t0_fine, T_fine );
            
            int iteration = 1;
            bool work = true;
            while( work )
            {
                heatfine.resetExporter(fmt::format("heat-fine-{}-{}", color, iteration));

                 
                if ( color > 1 )
                {
                    // receive initial guess from coarse integrator for time_interval
                    //mpi::irecv( );
                }
                for (double t = t0_fine; t < T_fine+dt_fine; t += dt_fine) {
                    if ( w->isMasterRank() )
                    {
                        LOG(INFO) << "====================================" << std::endl;
                        LOG(INFO) << fmt::format("Fine Integrator Time interval {}, t = {}", fine_time_interval, t) << std::endl;
                    }
                    // execute the time step: update the right hand side and solve the system
                    heatfine.run( t, heatfine.solution() );

                    // save solution at current time
                    heatfine.postProcess();
                   
                }   
                // send fine solution to coarse integrator
                // communication from fine to coarse integrators
                // non-blocking communication
                // send F(U^{k}_{T_fine}) to coarse integrator

                // update work flag to know whether we have to continue or not
                // communication from coarse to fine integrators
                // blocking communication
                // mpi::broadcast(wglob->globalComm(), 0, work);
                work = false;
                ++iteration;
            }
            LOG(INFO) << fmt::format("== fine integrator {} is finished ==================================",fine_time_interval ) << std::endl;
        } 
        #endif

        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
}