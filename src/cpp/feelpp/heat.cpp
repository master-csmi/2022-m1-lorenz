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

        Heat<2,1> heat( specs );

        for (double t = 0.0; t < 1.0; t += 0.1) {
            heat.run();
            heat.postProcess();
        }
        

        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
}