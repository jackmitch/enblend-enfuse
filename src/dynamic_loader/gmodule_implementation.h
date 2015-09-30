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
#ifndef GMODULE_IMPLEMENTATION_H_INCLUDED
#define GMODULE_IMPLEMENTATION_H_INCLUDED


#include "dynamic_loader_implementation.h"

#include "config.h"


#ifdef HAVE_GMODULE

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
    explicit GLibDynamicLoaderImplementation(const std::string& a_library_name);
    void open();
    void close();
    void* resolve(const std::string& symbol_name) const;

private:
    GModule* module_;
}; // class GLibDynamicLoaderImplementation

typedef GLibDynamicLoaderImplementation ActualDynamicLoaderImplementation;

#endif // HAVE_GMODULE

#endif // GMODULE_IMPLEMENTATION_H_INCLUDED


// Local Variables:
// mode: c++
// End:
