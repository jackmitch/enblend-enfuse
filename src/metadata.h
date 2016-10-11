/*
 * Copyright (C) 2016 Dr. Christoph L. Spiel
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


#ifdef HAVE_EXIV2
#include <exiv2/image.hpp>
#endif



namespace metadata
{
#ifdef HAVE_EXIV2
    Exiv2::Image::AutoPtr read(const std::string& an_image_filename);
    void write(const std::string& an_image_filename, Exiv2::Image::AutoPtr some_input_metadata);
#endif
} // namespace metadata


#endif // METADATA_H_INCLUDED

// Local Variables:
// mode: c++
// End:
