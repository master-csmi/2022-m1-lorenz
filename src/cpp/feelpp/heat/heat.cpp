#include <boost/mpi/group.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/range.hpp>
#include <feel/feelcore/utility.hpp>
#include <feelpp/heat.hpp>
#include <boost/mpi/group.hpp>
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
        mpi::communicator world;
        // number of time interval for parareal
        int time_partitions = specs["/Time/partitions"_json_pointer].get<int>();
        size_type space_partitions = static_cast<size_type>(world.size()/(time_partitions+1));
        LOG(INFO) << fmt::format("time partitions: {}, n time communicators: {}", time_partitions, time_partitions+1);

        // split the main worldcomm into P+1 subworldcomms
        auto [color, w, wglob] = Environment::worldCommPtr()->split(time_partitions+1);
        LOG(INFO) << fmt::format("color: {}, wglob rank: {}, w global rank: {}, w local rank: {}", color, wglob->globalRank(), w->globalRank(), w->localRank());

        // coarse time step
        //double dt = specs["/Time/dt"_json_pointer].get<double>();
        //double t0 = specs["/Time/t0"_json_pointer].get<double>();
        double T = specs["/Time/final_time"_json_pointer].get<double>();
        double t0 = 0.;
        LOG(INFO) << fmt::format("create  {} groups", space_partitions);

        // now we want to create K communicators between the coarse group and fine groups of processes
        // for each coarse process i we associate the fine process i in each fine group
        std::vector<int> pids;
        mpi::group world_group = world.group();        
        for ( int rank = 0; rank < world.size(); rank++ )
        {
            if ( rank / (time_partitions+1) == w->localRank() )
            {
                pids.push_back(rank);
            }
        }

        LOG(INFO) << fmt::format("-- global rank {} rank {} pids: {}\n", wglob->globalRank(), w->localRank(), pids);
        mpi::group coarsefine = world_group.include(pids.begin(), pids.end());
        LOG(INFO) << fmt::format("-- group {} created\n", w->localRank());
        mpi::communicator c(world, coarsefine );
        LOG(INFO) << fmt::format("communicator {} created\n", w->localRank());
        if ( c )
        {
            CHECK(c.size()==time_partitions+1) << fmt::format("wrong communicator size {} vs {}", c.size(), time_partitions+1);
            int sum = 1;
            sum = mpi::all_reduce(c, sum, std::plus<int>());
            LOG(INFO) << fmt::format("sum over communicator: {}\n", sum);
            CHECK(sum == time_partitions + 1) << fmt::format("on group {} the sum {} is not the number of time partitions + 1 {}", w->localRank(), sum, time_partitions + 1);
        }
        else
        {
            CHECK(false) << fmt::format("communicator not active rank {}, global rank {} pids {} - THIS SHOULD NOT HAPPEN \n", w->localRank(), wglob->globalRank(), pids);
        }
        LOG(INFO) << fmt::format("done with spatial group/communicator size:{}, rank:{}\n",c.size(),c.rank()) << std::endl;
        // create the worldcomm from new communicator
        // we can now send/receive data from spatial groups
        WorldComm wc(c);
#if 0
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
#endif

#if 1
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
                LOG(INFO) << "############################################" << std::endl;
                LOG(INFO) << fmt::format("########## Parareal Iteration = {}, rank: {}", iteration, w->localRank()) << std::endl;

                // initial condition is 0 : need to change that for different initial condition
                heatcoarse.solution().zero();
                heatcoarse.resetExporter(fmt::format("heat-coarse-{}-{}", color, iteration) );
                int k = 1;
                std::vector<mpi::request> reqs;
                for (double t = dt_coarse; t < T_coarse+dt_coarse; t += dt_coarse, ++k) 
                {
                    if ( w->isMasterRank() )
                    {
                        std::cout << "====================================" << std::endl;
                        std::cout << fmt::format("Coarse integrator t = {}", t) << std::endl;
                    }
                    LOG(INFO) << "====================================" << std::endl;
                    LOG(INFO) << fmt::format("Coarse Integrator Time interval t = {}, rank: {}", t, w->localRank()) << std::endl;

                    // execute the time step: update the right hand side and solve the system
                    // compute G(U^k_{j-1}), heatcoarse.solution() == U^k_{j-1}
                    heatcoarse.run(t, heatcoarse.solution());
                    if ( iteration > 1 )
                    {
                        heatcoarse.correction() = heatcoarse.load(t, iteration - 1, "correction");
                        heatcoarse.solution() += heatcoarse.correction();
                    }
                    // now heatcoarse.solution() == G(U^k_{j-1})
                    // save to disk to be able to reload afterward when we receive the correction from the fine integrators
                    heatcoarse.save(heatcoarse.solution(), t,iteration, "solution");

                    if ( std::abs(t-T_coarse) > 1e-10 )
                    {
                        LOG(INFO) << fmt::format("Coarse Integrator non blocking send initiated to {} tag {} size: {} \n", k+1, k, heatcoarse.solution().size()) << std::endl;
                        // receive initial guess from coarse integrator for time_interval
                        reqs.push_back( c.isend(k+1, k, heatcoarse.solution()) );
                    }

                    // non blocking async comm to receive fine integrator communication
                    // send U^k_j =  G(U^k_{j-1})+(F(U^{k-1}_{j-1})-G(U^{k-1}_{j-1}))
                    // U^k_j =  G(U^k_{j-1}) + correction[j-1]
                    // save solution at current time
                    
                }
                LOG(INFO) << fmt::format("Coarse Integrator non blocking send wait_all\n", reqs.size()) << std::endl;
                // this is blocking until we receive the initial guess
                mpi::wait_all(reqs.begin(), reqs.end());
                reqs.clear();

                // get in sync with the fine integrators for each coarse time step
                // we have a collection of size P of solution for each coarse time step 
                bool done = true;
                k = 1;
                for (double t = dt_coarse; t < T_coarse + dt_coarse; t += dt_coarse, ++k )
                {
                    // we are in coarse iteration j
                    auto sol = heatcoarse.load(t, iteration, "solution");
                    // Receive F(U^{k-1}_{j-1}) from the fine integrator
                    c.recv(k, k-1, heatcoarse.correction());
                    heatcoarse.correction() -= sol;
                    sync(heatcoarse.correction());
                    heatcoarse.save(heatcoarse.correction(),t,iteration,"correction");

                    heatcoarse.postProcess();

                    double err = 0;
                    if ( iteration > 1 )
                    {
                        auto sol_prev = heatcoarse.load(t, iteration-1, "solution");
                        // update work flag to know if we stop or continue
                        // the issue feelpp/feelpp#1489 kicks in, need to solve asap
                        err = 1e-5/std::pow(10,iteration); // normL2(_range=elements(mesh),_expr=idv(sol_prev)-idv(sol));
                        LOG(INFO) << fmt::format("Coarse Integrator error: {}\n", err ) << std::endl;
                        done = done && ( err < 1e-10 );
                    }
                    else
                        done = false;
                    //done = true;
                }
                work = !done;
                // broadcast work flag to all processors
                mpi::broadcast(wglob->globalComm(), work, 0);
                LOG(INFO) << fmt::format("work: {}", work) << std::endl;
                if ( work )
                    ++iteration;
            }
            LOG(INFO) << fmt::format("== coarse integrator is finished in {} iterations ==================================",iteration) << std::endl;
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
                LOG(INFO) << "############################################" << std::endl;
                LOG(INFO) << fmt::format("########## Parareal Iteration = {}, rank: {}", iteration, w->localRank()) << std::endl;


                heatfine.resetExporter(fmt::format("heat-fine-{}-{}", color, iteration));
                mpi::request reqs[1];

                if ( color > 1 )
                {
                    LOG(INFO) << fmt::format("Fine Integrator waiting for coarse integrator initial guess from time interval {}, size {}", color - 1, heatfine.solution().size()) << std::endl;
                    
                    // receive initial guess from coarse integrator for time_interval
                    reqs[0] = c.irecv( 0, fine_time_interval, heatfine.solution() );
                    LOG(INFO) << fmt::format("Fine Integrator non blocking receive initiated tag {} rank : {}\n",fine_time_interval,w->localRank()) << std::endl;
                    // this is blocking until we receive the initial guess
                    mpi::wait_all(reqs, reqs + 1);
                    LOG(INFO) << fmt::format("Fine Integrator non blocking receive completed\n") << std::endl;
                    sync(heatfine.solution());
                }
                else
                {
                    // initial condition is 0 at t=0: need to change that for different initial condition
                    heatcoarse.solution().zero();
                }
                
                for (double t = t0_fine; t < T_fine+dt_fine; t += dt_fine) 
                {
                    LOG(INFO) << "====================================" << std::endl;
                    LOG(INFO) << fmt::format("Fine Integrator Time interval {}, t = {}, rank: {}", fine_time_interval, t, w->localRank()) << std::endl;
                
                    // execute the time step: update the right hand side and solve the system
                    heatfine.run( t, heatfine.solution() );

                    // send the final time to the coarse process
                    if ( std::abs(t-T_fine) < 1e-10 )
                    {
                        LOG(INFO) << fmt::format("Fine Integrator non blocking send initiated tag {} rank : {}\n", fine_time_interval, w->localRank()) << std::endl;
                        reqs[0] = c.isend(0, fine_time_interval, heatfine.solution());
                    }
                    // save solution at current time
                    heatfine.postProcess();
                }  
                LOG(INFO) << fmt::format("done with time stepping\n");
                // send fine solution to coarse integrator
                // communication from fine to coarse integrators
                // non-blocking communication
                // send F(U^{k}_{T_fine}) to coarse integrator
                mpi::wait_all(reqs, reqs + 1);
                LOG(INFO) << fmt::format("Fine Integrator non blocking send completed\n") << std::endl;

                // update work flag to know whether we have to continue or not
                // communication from coarse to fine integrators
                // blocking communication
                LOG(INFO) << fmt::format("check if there is still work to do\n") << std::endl;
                mpi::broadcast(wglob->globalComm(), work, 0);
                LOG(INFO) << fmt::format("work: {}", work ) << std::endl;
                if ( work )
                    ++iteration;
            }
            LOG(INFO) << fmt::format("== fine integrator {} is finished in {} iterations ==================================",fine_time_interval, iteration ) << std::endl;
        }
#endif

        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
}