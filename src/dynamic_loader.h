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


#include <stdexcept>
#include <string>
#include <vector>
#include "global.h"  // for enblend::trim

class DynamicLoaderImplementation
{
public:
    DynamicLoaderImplementation() = delete;

    explicit DynamicLoaderImplementation(const std::string& a_library_name) :
        name_(a_library_name)
    {}

    virtual void open() = 0;
    virtual void close() = 0;
    virtual void* resolve(const std::string& a_symbol_name) const = 0;

    const std::string& library_name() const {return name_;}

    virtual ~DynamicLoaderImplementation() {}

    struct error : public std::runtime_error
    {
        error(const std::string& a_message) : std::runtime_error(a_message) {}
    };

private:
    const std::string name_;
};


#if defined(HAVE_DL)

#define HAVE_DYNAMICLOADER_IMPL
#include <dlfcn.h>

class SunnyDynamicLoaderImplementation : public DynamicLoaderImplementation
{
    typedef DynamicLoaderImplementation super;

public:
    explicit SunnyDynamicLoaderImplementation(const std::string& a_library_name) :
        super(a_library_name), handle_(nullptr)
    {}

    void open()
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

    void close()
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

    void* resolve(const std::string& symbol_name) const
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

private:
    void* handle_;
}; // class SunnyDynamicLoaderImplementation

typedef SunnyDynamicLoaderImplementation ActualDynamicLoaderImplementation;

#elif defined(HAVE_GMODULE)

// ANTICIPATED CHANGE: Class GLibDynamicLoader is completely untested,
// but the code looks pretty kosher.  To test use e.g.
//     EXTRACPPFLAGS="-I/usr/include/glib-2.0 -I/usr/include/glib-2.0/include -I/usr/lib/glib-2.0/include"
//     EXTRALDFLAGS="-lgmodule-2.0"
// undefine all relevant prior HAVE_* and define HAVE_GMODULE.

#define HAVE_DYNAMICLOADER_IMPL
#include <gmodule.h>

class GLibDynamicLoaderImplementation : public DynamicLoaderImplementation
{
    typedef DynamicLoaderImplementation super;

public:
    explicit GLibDynamicLoaderImplementation(const std::string& a_library_name) :
        super(a_library_name), module_(nullptr)
    {}

    void open()
    {
        module_ = g_module_open(library_name().c_str(), G_MODULE_BIND_LAZY);
        if (!module_)
        {
            super::error(g_module_error());
        }
    }

    void close()
    {
        if (!g_module_close(module_))
        {
            super::error(g_module_error());
        }
    }

    void* resolve(const std::string& symbol_name) const
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

private:
    GModule* module_;
}; // class GLibDynamicLoaderImplementation

typedef GLibDynamicLoaderImplementation ActualDynamicLoaderImplementation;

#elif defined(WIN32)

#define HAVE_DYNAMICLOADER_IMPL
#include <Windows.h>
#include <boost/lexical_cast.hpp>

class WinDynamicLoaderImplementation : public DynamicLoaderImplementation
{
    typedef DynamicLoaderImplementation super;

public:
    explicit WinDynamicLoaderImplementation(const std::string& a_library_name) :
        super(a_library_name), handle_(nullptr)
    {}

    void open()
    {
        if (handle_)
        {
            throw super::error("already open");
        }
        handle_ = LoadLibrary(library_name().c_str());
        if (!handle_)
        {
            throw super::error(GetLastErrorString());
        }
    }

    void close()
    {
        if (!handle_)
        {
            throw super::error("not open");
        }
        if (!FreeLibrary(handle_))
        {
            throw super::error(GetLastErrorString());
        }
    }

    void* resolve(const std::string& symbol_name) const
    {
        if (!handle_)
        {
            throw super::error("not open");

        }
        void* symbol = (void*) GetProcAddress(handle_, symbol_name.c_str());
        if (symbol == nullptr)
        {
            throw super::error(GetLastErrorString());
        }

        return symbol;
    }

private:
    static const std::string GetLastErrorString()
    {
        LPTSTR lpMsgBuf;
        DWORD lastError = GetLastError();
        std::string errorMsg;
        if (FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                          nullptr, lastError, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),  (LPTSTR) &lpMsgBuf, 0, nullptr) > 0)
        {
            errorMsg=lpMsgBuf;
            // remove leading and trailing white spaces
            enblend::trim(errorMsg);
            LocalFree(lpMsgBuf);
        }
        else
        {
            errorMsg = "Unknown error";
        }
        // add numeric error code
        errorMsg.append(" (Code: ");
        errorMsg.append(boost::lexical_cast<std::string>(lastError));
        errorMsg.append(")");
        return errorMsg;
    }

    HINSTANCE handle_;
}; // class WinDynamicLoaderImplementation

typedef WinDynamicLoaderImplementation ActualDynamicLoaderImplementation;

#else

class NullDynamicLoaderImplementation : public DynamicLoaderImplementation
{
    typedef DynamicLoaderImplementation super;

public:
    explicit NullDynamicLoaderImplementation(const std::string& a_library_name) :
        super(a_library_name)
    {}

    void open() {}
    void close() {}

    void* resolve(const std::string&) const
    {
        return nullptr;
    }
}; // class NullDynamicLoaderImplementation

typedef NullDynamicLoaderImplementation ActualDynamicLoaderImplementation;

#endif


class DynamicLoader
{
public:
    DynamicLoader() = delete;

    explicit DynamicLoader(const std::string& a_library_name) :
        implementation_(new ActualDynamicLoaderImplementation(a_library_name))
    {
        implementation_->open();
    }

    DynamicLoader(const DynamicLoader& another_dynamic_loader) :
        implementation_(new ActualDynamicLoaderImplementation(another_dynamic_loader.implementation_->library_name()))
    {
        implementation_->open();
        // Observer list of the new, copied instance is empty.
    }

    DynamicLoader& operator=(const DynamicLoader& another_dynamic_loader)
    {
        if (this != &another_dynamic_loader)
        {
            finalize();

            observers_ = another_dynamic_loader.observers_;
            implementation_ = another_dynamic_loader.implementation_;
        }

        return *this;
    }

    ~DynamicLoader()
    {
        finalize();
    }

    // Access symbols that do not require a teardown function to be
    // called on un-linking.
    void* resolve0(const std::string& a_symbol_name) const
    {
        return implementation_->resolve(a_symbol_name);
    }

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
    void* resolve0(const std::string& a_symbol_name, Teardown* a_teardown_object)
    {
        observers_.push_back(a_teardown_object);
        return resolve0(a_symbol_name);
    }

    template <typename T>
    T resolve(const std::string& a_symbol_name, Teardown* a_teardown_object)
    {
        return static_cast<T>(resolve0(a_symbol_name, a_teardown_object));
    }

private:
    void finalize()
    {
        for (observer_list::iterator x = observers_.begin(); x != observers_.end(); ++x)
        {
            (*x)->teardown(this);
        }

        implementation_->close();
        delete implementation_;
    }

    typedef std::vector<Teardown*> observer_list;

    ActualDynamicLoaderImplementation* implementation_;
    observer_list observers_;
}; // class DynamicLoader


#endif // DYNAMIC_LOADER_H_INCLUDED


// Local Variables:
// mode: c++
// End:
