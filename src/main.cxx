#include "preambler.h"
#include "appmanager.h"
#include "get_time_now.h"

#include "boost/log/core.hpp"
#include "boost/log/trivial.hpp"
#include "boost/log/sources/severity_logger.hpp"

#include <iostream>

namespace logging = boost::log;

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
    preamble.initialise_logging();
    boost::log::sources::severity_logger<logging::trivial::severity_level> lg;
    
    BOOST_LOG_SEV(lg, logging::trivial::info) << preamble.get_banner_message();

    const preamble::AppOptions* app_config = preamble.get_options();
    BOOST_LOG_SEV(lg, logging::trivial::info) << app_config->summarise();

    //----------------------------------------
    // Start work...
    bool errors_occured = false;
    AppManager appmanager;

    try
    {
        errors_occured = static_cast<bool>(appmanager.run());
    }
    catch(const std::exception& e)
    {
        BOOST_LOG_SEV(lg, logging::trivial::fatal) << "Could not recover from: " << e.what();
        appmanager.clean_up();

        BOOST_LOG_SEV(lg, logging::trivial::info) << "Terminated at " << get_time_now() << std::endl;
        return EXIT_FAILURE;
    }
    
    if (errors_occured)
        BOOST_LOG_SEV(lg, logging::trivial::error) << "Encountered errors and terminated at " << get_time_now() << "...";
    else
        BOOST_LOG_SEV(lg, logging::trivial::info) << "Tasks completed at " << get_time_now();

    //----------------------------------------
    // Exit success/failure - more info in logs generated at runtime
    return errors_occured ? EXIT_FAILURE : EXIT_SUCCESS;
}


