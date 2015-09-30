#include "sunny_implementation.h"


#ifdef HAVE_DL


SunnyDynamicLoaderImplementation::SunnyDynamicLoaderImplementation(const std::string& a_library_name) :
    super(a_library_name), handle_(nullptr)
{}


void
SunnyDynamicLoaderImplementation::open()
{
    if (handle_)
    {
        throw super::error("already open");
    }
    handle_ = dlopen(library_name().c_str(), RTLD_LAZY);
    if (!handle_)
    {
        throw super::error(dlerror());
    }
}


void
SunnyDynamicLoaderImplementation::close()
{
    if (!handle_)
    {
        throw super::error("not open");
    }
    if (dlclose(handle_) != 0)
    {
        throw super::error(dlerror());
    }
}


void*
SunnyDynamicLoaderImplementation::resolve(const std::string& symbol_name) const
{
    if (!handle_)
    {
        throw super::error("not open");
    }

    dlerror();
    void* symbol = dlsym(handle_, symbol_name.c_str());
    char* message = dlerror();
    if (message != nullptr)
    {
        throw super::error(message);
    }
    if (symbol == nullptr)
    {
        throw super::error("symbol points nowhere");
    }

    return symbol;
}


#endif // HAVE_DL
