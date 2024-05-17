#ifndef DYNLIBFUNCLOADER_H
#define DYNLIBFUNCLOADER_H

#ifdef _WIN32
#include <windows.h>
#define Dlsym(handle, symbol) GetProcAddress(handle, symbol)
#else if __linux__
#include <dlfcn.h>
#define Dlsym(handle, symbol) dlsym(handle, symbol)
#endif
#include "get_os_error.h"

#include <string>

template <typename func_t>
class DynLibFuncLoader
{
public:
    static explicit func_t load(const char* func_name, const char* libname, void* lib_handle) noexcept(false)
    {
        func_t l_func = (func_t)Dlsym(lib_handle, func_name);
        if (l_func == nullptr)
        {   // Error handling for symbol resolution
            std::string("Failed to load function: ") + std::string(func_name) + std::string(" from ") + 
                std::string(libname) + std::string(" (");

            std::string what_happened = get_os_error();

            message += what_happened + std::string(")");

            throw std::runtime_error(message);
        }
        return l_func;
    }
};

#endif // DYNLIBFUNCLOADER_H