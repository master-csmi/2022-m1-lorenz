#include <feelpp/heat.hpp>

int main(int argc, const char** argv) {
    try{
        using namespace Feel;
        Environment env(_argc = argc, _argv = argv,
                        _desc = makeOptions(),
                        _about = about(_name = "heat",
                                    _author = "Feel++ Consortium",
                                    _email = "feelpp@cemosis.fr"));
        auto jsonfile = removeComments(readFromFile(Environment::expand(soption("specs"))));
        std::istringstream istr(jsonfile);
        json specs = json::parse(istr);

        Heat<2,1> heat( specs );
        heat.run();
        heat.postProcess();

        return 0;
    }
    catch(...)
    {
        handleException();
    }
}