/*
 * Copyright (C) 2012-2016 Christoph L. Spiel
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


#include <memory>

#include "metadata.h"


namespace metadata
{
#ifdef HAVE_EXIV2
    static void
    copy_exif(Exiv2::ExifData& output_exif, const Exiv2::ExifData& input_exif)
    {
        for (auto x : input_exif)
        {
            if (output_exif.findKey(Exiv2::ExifKey(x.key())) == output_exif.end())
            {
                output_exif.add(x);
            }
        }
    }


    static void
    copy_iptc(Exiv2::IptcData& output_iptc, const Exiv2::IptcData& input_iptc)
    {
        for (auto x : input_iptc)
        {
            if (output_iptc.findKey(Exiv2::IptcKey(x.key())) == output_iptc.end())
            {
                output_iptc.add(x);
            }
        }
    }


    static void
    copy_xmp(Exiv2::XmpData& output_xmp, const Exiv2::XmpData& input_xmp)
    {
        for (auto x : input_xmp)
        {
            if (output_xmp.findKey(Exiv2::XmpKey(x.key())) == output_xmp.end())
            {
                output_xmp.add(x);
            }
        }
    }


    static void
    copy(Exiv2::Image* some_output_metadata, const Exiv2::Image* some_input_metadata)
    {
        copy_exif(some_output_metadata->exifData(), some_input_metadata->exifData());
        copy_iptc(some_output_metadata->iptcData(), some_input_metadata->iptcData());
        copy_xmp(some_output_metadata->xmpData(), some_input_metadata->xmpData());
    }


#if defined(ENBLEND_SOURCE)
#define PROCESSING_SOFTWARE "enblend"
#elif defined(ENFUSE_SOURCE)
#define PROCESSING_SOFTWARE "enfuse"
#else
#error "Neither ENBLEND_SOURCE nor ENFUSE_SOURCE are defined."
#endif


    static void
    augment(Exiv2::Image* some_output_metadata)
    {
        Exiv2::ExifData& output_exif {some_output_metadata->exifData()};

        std::unique_ptr<Exiv2::AsciiValue>
            processing_software {new Exiv2::AsciiValue(PROCESSING_SOFTWARE " " VERSION)};

        output_exif.add(Exiv2::ExifKey("Exif.Image.ProcessingSoftware"), processing_software.get());
    }


    Exiv2::Image::AutoPtr
    read(const std::string& an_image_filename)
    {
        Exiv2::Image::AutoPtr meta {Exiv2::ImageFactory::open(an_image_filename)};

        if (meta.get() && meta->good())
        {
            meta->readMetadata();
        }

        return meta;
    }


    void
    write(const std::string& an_image_filename, const Exiv2::Image::AutoPtr some_input_metadata)
    {
        Exiv2::Image::AutoPtr output_meta {Exiv2::ImageFactory::open(an_image_filename)};

        if (some_input_metadata.get() && some_input_metadata->good() &&
            output_meta.get() && output_meta->good())
        {
            copy(output_meta.get(), some_input_metadata.get());
            augment(output_meta.get());
            output_meta->writeMetadata();
        }
    }
#endif // HAVE_EXIV2
} // namespace metadata
