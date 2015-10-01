/*
 * Copyright (C) 2013, 2015 Dr. Christoph L. Spiel
 *
 * This file is part of Enblend.
 *
 * Enblend is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Enblend is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef DYNAMIC_LOADER_H_INCLUDED
#define DYNAMIC_LOADER_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#if defined(GMODULE_DL)
#define HAVE_DYNAMICLOADER_IMPL
#include "gmodule_implementation.h"        // g_module_open(3)
#elif defined(POSIX_DL)
#define HAVE_DYNAMICLOADER_IMPL
#include "posix_implementation.h"          // dlopen(3)
#elif defined(WIN32_DL)
#define HAVE_DYNAMICLOADER_IMPL
#include "win32_implementation.h"
#else
// This is the fallback class if no implementation class has been selected yet.
#include "null_implementation.h"
#endif


class DynamicLoader
{
public:
    DynamicLoader() = delete;

    explicit DynamicLoader(const std::string& a_library_name);
    DynamicLoader(const DynamicLoader& another_dynamic_loader);

    DynamicLoader& operator=(const DynamicLoader& another_dynamic_loader);

    ~DynamicLoader();

    // Access symbols that do not require a teardown function to be
    // called on un-linking.
    void* resolve0(const std::string& a_symbol_name) const;

    template <typename T>
    T resolve(const std::string& a_symbol_name) const
    {
        return static_cast<T>(resolve0(a_symbol_name));
    }

    class Teardown
    {
    public:
        Teardown() = delete;

        virtual void teardown(DynamicLoader*) = 0;
        virtual ~Teardown() {}
    };

    // Gain access to a symbol and simultaneously register a clean-up
    // object, which can e.g. run a clean-up function for the symbol.
    void* resolve0(const std::string& a_symbol_name, Teardown* a_teardown_object);

    template <typename T>
    T resolve(const std::string& a_symbol_name, Teardown* a_teardown_object)
    {
        return static_cast<T>(resolve0(a_symbol_name, a_teardown_object));
    }

private:
    void finalize();

    typedef std::vector<Teardown*> observer_list;

    ActualDynamicLoaderImplementation* implementation_;
    observer_list observers_;
}; // class DynamicLoader


#endif // DYNAMIC_LOADER_H_INCLUDED


// Local Variables:
// mode: c++
// End:
