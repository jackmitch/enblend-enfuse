/*
 * Copyright (C) 2016, 2017 Dr. Christoph L. Spiel
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

#ifndef METADATA_H_INCLUDED
#define METADATA_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdexcept>
#include <string>
#include <vector>


#ifdef HAVE_EXIV2
#include <exiv2/image.hpp>
#endif


namespace metadata
{
#ifdef HAVE_EXIV2
    struct Warning : public std::runtime_error
    {
        Warning() = delete;
        explicit Warning(const char* a_message) : std::runtime_error(a_message) {}
        explicit Warning(const std::string& a_message) : std::runtime_error(a_message) {}
    }; // struct Warning


    class Named
    {
        typedef Exiv2::Image::AutoPtr::element_type* meta_pointer;

    public:
        Named() = delete;
        Named(const std::string& a_filename, const meta_pointer a_meta, bool is_desired) :
            filename_(a_filename), meta_(a_meta), is_desired_(is_desired) {}

        const std::string& filename() const {return filename_;}
        meta_pointer meta() const {return meta_;}
        bool is_desired() const {return is_desired_;}

    private:
        const std::string& filename_;
        const meta_pointer meta_;
        bool is_desired_;
    }; // class Named


    typedef std::vector<Named> named_meta_array;


    Exiv2::Image::AutoPtr read(const std::string& an_image_filename);

    named_meta_array::const_iterator write(const std::string& an_image_filename,
                                           named_meta_array::const_iterator some_named_meta_begin,
                                           named_meta_array::const_iterator some_named_meta_end);
#endif
} // namespace metadata


#endif // METADATA_H_INCLUDED

// Local Variables:
// mode: c++
// End:
