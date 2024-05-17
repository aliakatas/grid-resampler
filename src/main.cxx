#include "appconfig.h"
#include "appmanager.h"
#include "build_version_utils.h"

#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/positional_options.hpp"

#include <iostream>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>

namespace bpo = boost::program_options;

std::string create_banner_message(int argc, char** argv);
int define_and_gather_cli_options(int argc, char** argv, bpo::variables_map& vm, std::string& s);
std::string get_time_now();

int main(int argc, char** argv)
{
    //----------------------------------------
    // Start with creating some information...
    std::string banner_message = create_banner_message(argc, argv);
    
    //----------------------------------------
    // Declare the supported options for the command line
    bpo::variables_map vm;
    std::string cli_input_options_desc;
    if (define_and_gather_cli_options(argc, argv, vm, cli_input_options_desc))
        return EXIT_SUCCESS;        // just printed help/version and that's all
    
    // otherwise, produce the banner
    std::cout << banner_message << std::endl;

    // also give a description of what was passed 
    // as input by the user (cli)
    std::cout << cli_input_options_desc << std::endl;

    //----------------------------------------
    // Start work...
    bool errors_occured = false;
    AppManager* appmanager = new AppManager;

    try
    {
        errors_occured = static_cast<bool>(appmanager->run());
    }
    catch(const std::exception& e)
    {
        std::cout << "[*Fatal*] Could not recover from: \n";
        std::cout << e.what() << "\n" << std::endl;
        delete appmanager;

        std::cout << "\nTerminated at " << get_time_now() << std::endl;
        return EXIT_FAILURE;
    }
    
    if (errors_occured)
        std::cout << "\n[*Error*] Encountered errors and terminated at " << get_time_now() << "...\n" << std::endl;
    else
        std::cout << "\nTasks completed at " << get_time_now() << "\n" << std::endl;

    //----------------------------------------
    // Exit success/failure - more info in logs generated at runtime
    return errors_occured ? EXIT_FAILURE : EXIT_SUCCESS;
}

//======================================================
std::string create_banner_message(int argc, char** argv) 
{
    std::string start_time = get_time_now();
    
    std::string version_num = std::to_string(APP_VERSION_MAJOR) + std::string(".") + 
        std::to_string(APP_VERSION_MINOR) + std::string(".") + std::to_string(APP_VERSION_PATCH);
    
    std::string start_message = std::string("Starting ") + std::filesystem::path(argv[0]).stem().string() + 
        std::string(" ") + version_num + std::string(" at ") + start_time + std::string("\n");

    std::string compiler_info = std::string("Compiler: ") + build_info::get_compiler_info();
    std::string build_date = std::string("Build date: ") + build_info::get_build_date_time();
    std::string cuda_support = std::string("CUDA support: ") + build_info::get_nvidia_cuda_version();

    int nstars = std::max({start_message.length(), compiler_info.length(), build_date.length(), cuda_support.length()}) + 1;
    std::string star_buffer;
    for (auto i = 0; i < nstars; ++i)
        star_buffer += std::string("*");
    star_buffer += std::string("\n");
    
    std::string banner_message = std::string("\n") + star_buffer + 
        start_message + std::string("\n") +
        compiler_info + std::string("\n") +
        build_date + std::string("\n") +
        cuda_support + std::string("\n") +
        star_buffer;

    return banner_message;
}

//======================================================
int define_and_gather_cli_options(int argc, char** argv, bpo::variables_map& vm, std::string& s)
{
    int opt;
    bpo::options_description desc("Valid options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "print version string")
        ("compression,c", bpo::value<int>(), "set compression level")
        ("optimisation,o", bpo::value<int>(&opt)->default_value(10), "optimisation level")
        ("include-path,I", bpo::value<std::vector<std::string>>(), "include path")
        ("input-file", bpo::value<std::vector<std::string>>(), "input file")
    ;

    bpo::positional_options_description p;
    p.add("input-file", -1);

    bpo::store(bpo::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    bpo::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm.count("version")) {
        std::string version_num = std::to_string(APP_VERSION_MAJOR) + std::string(".") + 
            std::to_string(APP_VERSION_MINOR) + std::string(".") + std::to_string(APP_VERSION_PATCH);
        std::cout << version_num << "\n";
        return 1;
    }

    if (vm.count("compression")) 
        s += std::string("Compression level was set to ") + std::to_string(vm["compression"].as<int>()) + std::string(".\n");
    else 
        s += std::string("Compression level was not set.\n");

    if (vm.count("include-path"))
    {
        s += std::string("Include paths are: \n"); 
        for (auto& incp : vm["include-path"].as<std::vector<std::string>>())
            s += std::string("  ") + incp + std::string("\n");
    }

    if (vm.count("input-file"))
    {
        s += std::string("Input files are: \n");
        for (auto& inpf : vm["input-file"].as<std::vector<std::string>>())
            s += std::string("  ") + inpf + std::string("\n");
    }

    s += std::string("Optimization level is ") + std::to_string(opt) + std::string("\n"); 
    return 0;
}

//======================================================
std::string get_time_now()
{
    auto now = std::chrono::system_clock::now();
    time_t tt = std::chrono::system_clock::to_time_t(now);
    struct tm *localTime = localtime(&tt);
    char time_buffer_str[80];

    strftime(time_buffer_str, sizeof(time_buffer_str), "%Y-%m-%d %H:%M:%S", localTime);

    return std::string(time_buffer_str);
}
