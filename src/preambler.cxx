#include "preambler.h"
#include "get_time_now.h"
#include "appconfig.h"
#include "build_version_utils.h"

#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/positional_options.hpp"

#include <filesystem>
#include <sstream>

namespace preamble
{
    //####################################################
    Preambler::Preambler(int argc, char** argv)
    {
        _create_banner_message(argc, argv);
        _sort_out_app_cli_options(argc, argv);
    }

    //####################################################
    Preambler::~Preambler()
    {

    }

    //####################################################
    bool Preambler::only_help() const
    {
        return m_only_help;
    }

    //####################################################
    std::string Preambler::get_help() const
    {
        std::ostringstream oss;
        oss << m_cli_options;

        return oss.str();
    }

    //####################################################
    bool Preambler::only_version() const
    {
        return m_only_version;
    }

    //####################################################
    std::string Preambler::get_version() const
    {
        return m_version;
    }

    //####################################################
    std::string Preambler::get_banner_message() const
    {
        return m_banner_message;
    }

    //####################################################
    void Preambler::_create_banner_message(int argc, char** argv)
    {
        std::string start_time = get_time_now();
        
        m_version = std::to_string(APP_VERSION_MAJOR) + std::string(".") + 
            std::to_string(APP_VERSION_MINOR) + std::string(".") + std::to_string(APP_VERSION_PATCH);
        
        std::string start_message = std::string("Starting ") + std::filesystem::path(argv[0]).stem().string() + 
            std::string(" ") + m_version + std::string(" at ") + start_time + std::string("\n");

        std::string compiler_info = std::string("Compiler: ") + build_info::get_compiler_info();
        std::string build_date = std::string("Build date: ") + build_info::get_build_date_time();
        std::string cuda_support = std::string("CUDA support: ") + build_info::get_nvidia_cuda_version();

        int nstars = std::max({start_message.length(), compiler_info.length(), build_date.length(), cuda_support.length()}) + 1;
        std::string star_buffer;
        for (auto i = 0; i < nstars; ++i)
            star_buffer += std::string("*");
        star_buffer += std::string("\n");
        
        m_banner_message = std::string("\n") + star_buffer + 
            start_message + std::string("\n") +
            compiler_info + std::string("\n") +
            build_date + std::string("\n") +
            cuda_support + std::string("\n") +
            star_buffer;
    }

    //####################################################
    void Preambler::_sort_out_app_cli_options(int argc, char** argv)
    {
        std::string input_file;
        std::string output_file;
        std::string log_file;

        m_cli_options.add_options()
            ("help,h", "produce help message and exit")
            ("version,v", "print version string and exit")
            ("compression-method,cm", bpo::value<std::string>(&m_app_options.compression_method)->default_value(RASTER_NO_COMPRESSION), "set compression method")
            ("compression-level,cl", bpo::value<int>(&m_app_options.compression_level)->default_value(5), "set compression level")
            ("output-resolution,r", bpo::value<double>(&m_app_options.resolution_scale)->default_value(1.), "set the resolution factor for output grid")
            ("input-file,I", bpo::value<std::string>(&input_file), "input file path")
            ("output-file,O", bpo::value<std::string>(&output_file), "output file path")
            ("log,l", bpo::value<std::string>(&log_file), "log file path")
        ;

        bpo::positional_options_description p;
        p.add("input-file", -1);

        bpo::store(bpo::command_line_parser(argc, argv).options(m_cli_options).positional(p).run(), m_variables_map);
        bpo::notify(m_variables_map);

        m_only_help = (m_variables_map.count("help") > 0) || argc == 1;
        m_only_version = (m_variables_map.count("version") > 0);

        m_app_options.input_file = std::filesystem::path(input_file);
        m_app_options.output_file = std::filesystem::path(output_file);
    }
}
