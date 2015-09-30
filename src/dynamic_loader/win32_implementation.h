/*
 * Copyright (C) 2015 Thomas Modes
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
#ifndef WIN32_IMPLEMENTATION_H_INCLUDED
#define WIN32_IMPLEMENTATION_H_INCLUDED


#include "dynamic_loader_implementation.h"

#include "config.h"


#ifdef WIN32

#define HAVE_DYNAMICLOADER_IMPL

#include <Windows.h>

class WinDynamicLoaderImplementation : public DynamicLoaderImplementation
{
    typedef DynamicLoaderImplementation super;

public:
    explicit WinDynamicLoaderImplementation(const std::string& a_library_name);
    void open();
    void close();
    void* resolve(const std::string& symbol_name) const;

private:
    static const std::string GetLastErrorString();

    HINSTANCE handle_;
}; // class WinDynamicLoaderImplementation

typedef WinDynamicLoaderImplementation ActualDynamicLoaderImplementation;

#endif // WIN32

#endif // WIN32_IMPLEMENTATION_H_INCLUDED


// Local Variables:
// mode: c++
// End:
