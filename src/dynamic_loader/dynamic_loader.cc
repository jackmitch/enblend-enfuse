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

#include "dynamic_loader.h"


DynamicLoader::DynamicLoader(const std::string& a_library_name) :
    implementation_(new ActualDynamicLoaderImplementation(a_library_name))
{
    implementation_->open();
}


DynamicLoader::DynamicLoader(const DynamicLoader& another_dynamic_loader) :
    implementation_(new ActualDynamicLoaderImplementation(another_dynamic_loader.implementation_->library_name()))
{
    implementation_->open();
    // Observer list of the new, copied instance is empty.
}


DynamicLoader&
DynamicLoader::operator= (const DynamicLoader& another_dynamic_loader)
{
    if (this != &another_dynamic_loader)
    {
        finalize();

        observers_ = another_dynamic_loader.observers_;
        implementation_ = another_dynamic_loader.implementation_;
    }

    return *this;
}


DynamicLoader::~DynamicLoader()
{
    finalize();
}


// Access symbols that do not require a teardown function to be
// called on un-linking.
void*
DynamicLoader::resolve0(const std::string& a_symbol_name) const
{
    return implementation_->resolve(a_symbol_name);
}


// Gain access to a symbol and simultaneously register a clean-up
// object, which can e.g. run a clean-up function for the symbol.
void*
DynamicLoader::resolve0(const std::string& a_symbol_name, Teardown* a_teardown_object)
{
    observers_.push_back(a_teardown_object);
    return resolve0(a_symbol_name);
}


void
DynamicLoader::finalize()
{
    for (observer_list::iterator x = observers_.begin(); x != observers_.end(); ++x)
    {
        (*x)->teardown(this);
    }

    implementation_->close();
    delete implementation_;
}
