/*
 * Copyright (C) 2015-2017 Thomas Modes
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define NOMINMAX

#include <Windows.h>

#include "dynamic_loader_implementation.h"


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


class ActualDynamicLoaderImplementation : public WinDynamicLoaderImplementation
{
public:
    explicit ActualDynamicLoaderImplementation(const std::string& a_library_name) :
        WinDynamicLoaderImplementation(a_library_name)
    {}
}; // class ActualDynamicLoaderImplementation


#endif // WIN32_IMPLEMENTATION_H_INCLUDED


// Local Variables:
// mode: c++
// End:
