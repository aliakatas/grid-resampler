#ifndef PREAMBLER_H
#define PREAMBLER_H

#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/options_description.hpp"

#include <string>
#include <filesystem>

#define RASTER_NO_COMPRESSION       "NONE"
#define RASTER_LZW_COMPRESSION      "LZW"
#define RASTER_DEFLATE_COMPRESSION  "DEFLATE"

namespace preamble
{
    namespace bpo = boost::program_options;

    //+++++++++++++++++++++++++++++++++++++++++
    void normalise_path(const std::filesystem::path& cwd, std::filesystem::path& path);

    //+++++++++++++++++++++++++++++++++++++++++
    void create_file_parent_directories(const std::filesystem::path& file_path) noexcept(false);

    //+++++++++++++++++++++++++++++++++++++++++
    struct AppOptions
    {
        std::string compression_method = RASTER_NO_COMPRESSION;
        int compression_level = 0;
        double resolution_scale = 1.;
        std::filesystem::path input_file;
        std::filesystem::path output_file;
        std::filesystem::path log_file;

        void process_paths() noexcept(false);

        std::string summarise() const;
    };

    //+++++++++++++++++++++++++++++++++++++++++
    class Preambler
    {
    public:
        Preambler();

        Preambler(int argc, char** argv) noexcept(false);

        ~Preambler();

        void prepare(int argc, char** argv) noexcept(false);

        bool only_help() const;
        std::string get_help() const;

        bool only_version() const;
        std::string get_version() const;
        std::string get_banner_message() const;

        const AppOptions* get_options() const;

        void initialise_logging() const;

    private:
        std::string m_banner_message;
        std::string m_version;
        bpo::variables_map m_variables_map;
        bpo::options_description m_cli_options;

        AppOptions m_app_options;

        bool m_only_help = false;
        bool m_only_version = false;

        void _create_banner_message(int argc, char** argv);

        void _sort_out_app_cli_options(int argc, char** argv) noexcept(false);
    };
}

#endif // PREAMBLER_H