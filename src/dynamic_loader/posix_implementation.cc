/*
 * Copyright (C) 2013-2017 Dr. Christoph L. Spiel
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


// See: http://pubs.opengroup.org/onlinepubs/9699919799/functions/dlopen.html


#include "posix_implementation.h"


PosixDynamicLoaderImplementation::PosixDynamicLoaderImplementation(const std::string& a_library_name) :
    super(a_library_name), handle_(nullptr)
{}


void
PosixDynamicLoaderImplementation::open()
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
PosixDynamicLoaderImplementation::close()
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
PosixDynamicLoaderImplementation::resolve(const std::string& symbol_name) const
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
