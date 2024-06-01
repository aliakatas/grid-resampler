#include "preambler.h"
// #include "appconfig.h"
// #include "appmanager.h"
// #include "build_version_utils.h"

// #include "boost/program_options/options_description.hpp"
// #include "boost/program_options/variables_map.hpp"
// #include "boost/program_options/parsers.hpp"
// #include "boost/program_options/positional_options.hpp"

#include <iostream> 
// #include <filesystem>
// #include <chrono>
// #include <ctime>
// #include <stdexcept>
// #include <string>
// #include <vector>

// namespace bpo = boost::program_options;

int main(int argc, char** argv)
{
    preamble::Preambler preamble;

    try
    {
        preamble.prepare(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';

        std::cout << preamble.get_version() << std::endl;
        std::cout << preamble.get_help() << std::endl;

        return EXIT_FAILURE;
    }
    
    if (preamble.only_help())
    {
        std::cout << preamble.get_help() << std::endl;
        return EXIT_SUCCESS;
    }

    if (preamble.only_version())
    {
        std::cout << preamble.get_version() << std::endl;
        return EXIT_SUCCESS;
    }
    
    // Before proceeding with messages, sort out any logs 
    const preamble::AppOptions* app_config = preamble.get_options();


    std::cout << app_config->summarise() << std::endl;



    
    // //----------------------------------------
    // // Declare the supported options for the command line
    // bpo::variables_map vm;
    // std::string cli_input_options_desc;
    // if (define_and_gather_cli_options(argc, argv, vm, cli_input_options_desc))
    //     return EXIT_SUCCESS;        // just printed help/version and that's all
    
    // // otherwise, produce the banner
    // std::cout << banner_message << std::endl;

    // // also give a description of what was passed 
    // // as input by the user (cli)
    // std::cout << cli_input_options_desc << std::endl;

    //----------------------------------------
    // Start work...
    bool errors_occured = false;
    // AppManager* appmanager = new AppManager;

    // try
    // {
    //     errors_occured = static_cast<bool>(appmanager->run());
    // }
    // catch(const std::exception& e)
    // {
    //     std::cout << "[*Fatal*] Could not recover from: \n";
    //     std::cout << e.what() << "\n" << std::endl;
    //     delete appmanager;

    //     std::cout << "\nTerminated at " << get_time_now() << std::endl;
    //     return EXIT_FAILURE;
    // }
    
    // if (errors_occured)
    //     std::cout << "\n[*Error*] Encountered errors and terminated at " << get_time_now() << "...\n" << std::endl;
    // else
    //     std::cout << "\nTasks completed at " << get_time_now() << "\n" << std::endl;

    //----------------------------------------
    // Exit success/failure - more info in logs generated at runtime
    return errors_occured ? EXIT_FAILURE : EXIT_SUCCESS;
}


