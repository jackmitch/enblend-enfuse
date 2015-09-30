#include "gmodule_implementation.h"


#ifdef HAVE_GMODULE


// ANTICIPATED CHANGE: Class GLibDynamicLoader is completely untested,
// but the code looks pretty kosher.  To test use e.g.
//     EXTRACPPFLAGS="-I/usr/include/glib-2.0 -I/usr/include/glib-2.0/include -I/usr/lib/glib-2.0/include"
//     EXTRALDFLAGS="-lgmodule-2.0"
// undefine all relevant prior HAVE_* and define HAVE_GMODULE.

GLibDynamicLoaderImplementation::GLibDynamicLoaderImplementation(const std::string& a_library_name) :
    super(a_library_name), module_(nullptr)
{}


void
GLibDynamicLoaderImplementation::open()
{
    module_ = g_module_open(library_name().c_str(), G_MODULE_BIND_LAZY);
    if (!module_)
    {
        super::error(g_module_error());
    }
}


void
GLibDynamicLoaderImplementation::close()
{
    if (!g_module_close(module_))
    {
        super::error(g_module_error());
    }
}


void*
GLibDynamicLoaderImplementation::resolve(const std::string& symbol_name) const
{
    void* symbol;

    if (!g_module_symbol(module_, symbol_name.c_str(), &symbol))
    {
        super::error(g_module_error());
    }
    if (symbol == nullptr)
    {
        throw super::error("symbol points nowhere");
    }

    return symbol;
}


#endif // HAVE_GMODULE
