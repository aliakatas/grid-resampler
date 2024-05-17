#ifndef DYNLIBLOADER_H
#define DYNLIBLOADER_H

#ifdef _WIN32
#include <windows.h>
#define Dlopen(name, flags) LoadLibraryA(name)
#define Dlclose FreeLibrary
#define Dlerror GetLastError
#else if __linux__
#include <dlfcn.h>
#define Dlopen(name, flags) dlopen(name, flags)
#define Dlclose dlclose
#define Dlerror dlerror
#endif

#include <string>

/**
 * @brief Cross-platform loader for dynamic libraries.
*/
class DynLibLoader
{
private:
    void* p_lib_handle;
public:

    /**
     * @brief Loads the dynamic library in memory.
     * 
     * @param libname [in]
     * @param load_options [in] Ignored for Windows.
     * 
     * @throws Runtime error if the load fails.
    */
    explicit DynLibLoader(const char* libname, const int load_options) noexcept(false);

    /**
     * @brief Loads the dynamic library in memory.
     * 
     * @param libname [in]
     * @param load_options [in] Ignored for Windows.
     * 
     * @throws Runtime error if the load fails.
    */
    explicit DynLibLoader(const std::string& libname, const int load_options) noexcept(false);
    
    /**
     * @brief Releases the handle to the library.
    */
    ~DynLibLoader() noexcept;

    /**
     * @brief Returns the handle to the dynamic library.
    */
    void* Handle() const;
};

#endif // DYNLIBLOADER_H
