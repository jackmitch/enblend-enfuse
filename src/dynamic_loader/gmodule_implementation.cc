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


// See: https://developer.gnome.org/glib/unstable/glib-Dynamic-Loading-of-Modules.html


#include "gmodule_implementation.h"


GLibDynamicLoaderImplementation::GLibDynamicLoaderImplementation(const std::string& a_library_name) :
    super(a_library_name), module_(nullptr)
{}


void
GLibDynamicLoaderImplementation::open()
{
    module_ = g_module_open(library_name().c_str(), G_MODULE_BIND_LAZY);
    if (!module_)
    {
        super::error(static_cast<const char*>(g_module_error()));
    }
}


void
GLibDynamicLoaderImplementation::close()
{
    if (!g_module_close(module_))
    {
        super::error(static_cast<const char*>(g_module_error()));
    }
}


void*
GLibDynamicLoaderImplementation::resolve(const std::string& symbol_name) const
{
    void* symbol;

    if (!g_module_symbol(module_, symbol_name.c_str(), &symbol))
    {
        super::error(static_cast<const char*>(g_module_error()));
    }
    if (symbol == nullptr)
    {
        throw super::error("symbol points nowhere");
    }

    return symbol;
}
