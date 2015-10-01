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
#ifndef SUNNY_IMPLEMENTATION_H_INCLUDED
#define SUNNY_IMPLEMENTATION_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <dlfcn.h>

#include "dynamic_loader_implementation.h"


class SunnyDynamicLoaderImplementation : public DynamicLoaderImplementation
{
    typedef DynamicLoaderImplementation super;

public:
    explicit SunnyDynamicLoaderImplementation(const std::string& a_library_name);
    void open();
    void close();
    void* resolve(const std::string& symbol_name) const;

private:
    void* handle_;
}; // class SunnyDynamicLoaderImplementation

typedef SunnyDynamicLoaderImplementation ActualDynamicLoaderImplementation;


#endif // SUNNY_IMPLEMENTATION_H_INCLUDED


// Local Variables:
// mode: c++
// End:
