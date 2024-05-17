#ifndef BUILD_VERSION_UTILS_H 
#define BUILD_VERSION_UTILS_H  

#include <string>

/**
 * @brief Implements functions that help with general information
 * about the build of the application.
 * 
 */
namespace build_info {

    /**
     * @brief Get the compiler information used to
     * build the module where the function is called
     * from.
     * 
     * @return std::string 
     */
    std::string get_compiler_info();

    /**
     * @brief Get the build date & time of the module
     * where the function is called from.
     * 
     * @return std::string 
     */
    std::string get_build_date_time();

    /**
     * @brief Get the Nvidia CUDA version used to build 
     * the module where the function is called from.
     * 
     * @return std::string 
     */
    std::string get_nvidia_cuda_version();

}

#endif // BUILD_VERSION_UTILS_H
