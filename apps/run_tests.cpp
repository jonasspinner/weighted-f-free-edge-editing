//
// Created by jonas on 19.08.19.
//


#include <boost/program_options.hpp>

#include "../src/tests/FinderTests.h"
#include "../src/tests/GraphTests.h"
#include "../src/tests/SubgraphTests.h"
#include "../src/tests/EditorTests.h"
#include "../src/tests/PermutationTest.h"


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    int seed = 0;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("seed", po::value<int>(&seed)->default_value(seed), "seed for randomized tests")
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    FinderTests(seed).run();
    //GraphTests(seed).run();
    //SubgraphTests(seed).run();
    //EditorTests(seed).run();
    //PermutationTest(seed).run();


    return 0;
}

