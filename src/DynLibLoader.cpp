#include "DynLibLoader.h"
#include "get_os_error.h"

#include <stdexcept>
#include <string>

//=====================================================================
explicit DynLibLoader::DynLibLoader(const char* libname, const int load_options) noexcept(false)
{
    p_lib_handle = Dlopen(libname, load_options);
    if (p_lib_handle == nullptr)
    {   // Error handling for loading the library
        std::string message = std::string("Failed to open library: ") + std::string(libname) + std::string(" (");

        std::string what_happened = get_os_error();

        message += what_happened + std::string(")");

        throw std::runtime_error(message);
    }
}

//=====================================================================
explicit DynLibLoader::DynLibLoader(const std::string& libname, const int load_options) noexcept(false) : 
    DynLibLoader(libname.c_str(), load_options)
{}

//=====================================================================
DynLibLoader::~DynLibLoader() noexcept
{   
    Dlclose(p_lib_handle);
}

//=====================================================================
void* DynLibLoader::Handle() const
{
    return p_lib_handle;
}