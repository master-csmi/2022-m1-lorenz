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
        double dt = specs["/Time/dt"_json_pointer].get<double>();
        //double t0 = specs["/Time/t0"_json_pointer].get<double>();
        double T = specs["/Time/final_time"_json_pointer].get<double>();

        Heat<2,1> heat( specs );

        // TODO : implement the parareal method
        for (double t = dt; t < T; t += dt) {
            std::cout << "====================================" << std::endl;
            std::cout << fmt::format("t = {}", t) << std::endl;
            // execute the time step: update the right hand side and solve the system
            heat.run( t, heat.solution() );

            // save solution at current time
            heat.postProcess();
        }
        

        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
}