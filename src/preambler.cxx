#include "preambler.h"
#include "get_time_now.h"
#include "appconfig.h"
#include "build_version_utils.h"

#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/positional_options.hpp"

#include "boost/log/core.hpp"
#include "boost/log/trivial.hpp"
#include "boost/log/expressions.hpp"
#include "boost/log/sinks/text_file_backend.hpp"
#include "boost/log/utility/setup/file.hpp"
#include "boost/log/utility/setup/common_attributes.hpp"
#include "boost/log/sources/severity_logger.hpp"
#include "boost/log/utility/setup/console.hpp"

#include <filesystem>
#include <sstream>
#include <stdexcept>

namespace preamble
{
    namespace logging = boost::log;
    namespace src = boost::log::sources;
    namespace sinks = boost::log::sinks;
    namespace keywords = boost::log::keywords;

    //####################################################
    void normalise_path(const std::filesystem::path& cwd, std::filesystem::path& path)
    {
        if (path.is_relative())
        {
            path = (cwd / path).lexically_normal();
        }
    }

    //####################################################
    void create_file_parent_directories(const std::filesystem::path& file_path) noexcept(false)
    {
        std::filesystem::path output_parent = file_path.parent_path();

        if (!std::filesystem::exists(output_parent))
        {
            if (!std::filesystem::create_directories(output_parent))
            {
                std::string t = std::string(__PRETTY_FUNCTION__) +
                    std::string(": Failed to create directory(-ies): ") + output_parent.string();
                throw std::runtime_error(t.c_str());
            }
        }
    }

    //####################################################
    void AppOptions::process_paths() noexcept(false)
    {
        std::filesystem::path cwd = std::filesystem::current_path();

        if (input_file.string().length() == 0)
        {
            std::string t = std::string(__PRETTY_FUNCTION__) + std::string(": Input file is missing");
            throw std::runtime_error(t.c_str());
        }

        if (output_file.string().length() == 0)
        {
            output_file = input_file.parent_path();
            std::string new_name = input_file.stem().string() + std::string("_resampled");
            output_file = output_file / std::filesystem::path(new_name);
            output_file.replace_extension(".tif");
        }

        normalise_path(cwd, input_file);
        normalise_path(cwd, output_file);

        if (log_file.string().length() > 0)
        {
            normalise_path(cwd, log_file);
        }
        
        if (!std::filesystem::exists(input_file))
        {
            std::string t = std::string(__PRETTY_FUNCTION__) + 
                std::string(": ") + input_file.string() + std::string(" does not exist");
            throw std::runtime_error(t.c_str());
        }

        create_file_parent_directories(output_file);

        if (log_file.string().length() > 0)
        {
            create_file_parent_directories(log_file);
        }
    }

    //####################################################
    std::string AppOptions::summarise() const
    {
        std::string temp = "App Options summary: \n";
        temp += "   Input file:   " + input_file.string() + "\n";
        temp += "   Output file:  " + output_file.string() + "\n";
        temp += "   Log file:     ";
        if (log_file.string().length() > 0)
        {
            temp += log_file.string() + "\n";
        }
        else 
        {
            temp += "stdout \n";
        }
        temp += "   Compression:  " + compression_method + " (" + std::to_string(compression_level) + ")\n";
        temp += "   Output scale: " + std::to_string(resolution_scale) + "\n";

        return temp;
    }

    //####################################################
    Preambler::Preambler()
    {

    }

    //####################################################
    Preambler::Preambler(int argc, char** argv) noexcept(false)
    {
        prepare(argc, argv);
    }

    //####################################################
    Preambler::~Preambler()
    {

    }

    //####################################################
    void Preambler::prepare(int argc, char** argv) noexcept(false)
    {
        _create_banner_message(argc, argv);
        _sort_out_app_cli_options(argc, argv);
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
    const AppOptions* Preambler::get_options() const
    {
        return &m_app_options;
    }

    //####################################################
    void Preambler::initialise_logging() const
    {
        logging::add_console_log(
            std::cout,
            keywords::format = "[%Severity%]: %Message%"
        );
        logging::core::get()->set_filter
        (
            logging::trivial::severity >= logging::trivial::info
        );

        if (m_app_options.log_file.string().length() > 0)
        {
            logging::add_file_log
            (
                keywords::file_name = m_app_options.log_file.string().c_str(),
                keywords::format = "[%TimeStamp%][%Severity%]: %Message%"
            );

            logging::core::get()->set_filter
            (
                logging::trivial::severity >= logging::trivial::info
            );
        }

        logging::add_common_attributes();
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
    void Preambler::_sort_out_app_cli_options(int argc, char** argv) noexcept(false)
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

        m_app_options.process_paths();
    }
}
