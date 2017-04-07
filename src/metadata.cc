/*
 * Copyright (C) 2016, 2017 Christoph L. Spiel
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


#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <regex>
#include <sstream>
#include <vector>

#include "parameter.h"

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


    static const std::array<const char*, 1> static_blacklist_keys
    {
        R"(Exif\.Image\.[XY]Resolution)"
    };


    using regexp_list = std::vector<std::regex>;


    static void
    compile_blacklist(std::back_insert_iterator<regexp_list> a_blacklist)
    {
        const std::regex::flag_type regex_flags {std::regex::icase | std::regex::nosubs};

        // Static part.
        for (auto x : static_blacklist_keys)
        {
            *a_blacklist++ = std::regex(x, regex_flags);
        }

        // Dynamic part.
        const std::regex delimiter {R"(\s*[,;:]\s*)"};
        const std::string dynamic_keys {parameter::as_string("blacklist-exif-keys", "")};
        std::sregex_token_iterator
            dynamic_key_iterator(dynamic_keys.begin(), dynamic_keys.end(), delimiter, -1);
        while (dynamic_key_iterator != std::sregex_token_iterator())
        {
            const std::string key {dynamic_key_iterator->str()};
            if (!key.empty())
            {
                *a_blacklist++ = std::regex(key, regex_flags);
            }
            ++dynamic_key_iterator;
        }
    }


    static void
    enforce_blacklist(Exiv2::Image* some_output_metadata,
                      regexp_list::const_iterator a_blacklist_begin,
                      regexp_list::const_iterator a_blacklist_end)
    {
        Exiv2::ExifData& output_exif {some_output_metadata->exifData()};

        auto exif_datum = output_exif.begin();
        while (exif_datum != output_exif.end())
        {
            for (auto x = a_blacklist_begin; x != a_blacklist_end; ++x)
            {
                if (std::regex_match(exif_datum->key(), *x))
                {
#ifdef DEBUG_METADATA
                    std::cout << "+ EXIF key `" << exif_datum->key() << "' is blacklisted -- not copied\n";
#endif
                    exif_datum = output_exif.erase(exif_datum);
                }
                else
                {
                    ++exif_datum;
                }
            }
        }
    }


    Exiv2::Image::AutoPtr
    read(const std::string& an_image_filename)
    {
        Exiv2::Image::AutoPtr meta {Exiv2::ImageFactory::open(an_image_filename)};

        if (meta.get() && meta->good())
        {
            meta->readMetadata();
            return meta;
        }
        else
        {
            return Exiv2::Image::AutoPtr(nullptr);
        }
    }


    named_meta_array::const_iterator
    write(const std::string& an_image_filename,
          named_meta_array::const_iterator some_named_meta_begin,
          named_meta_array::const_iterator some_named_meta_end)
    {
        Exiv2::Image::AutoPtr output_meta {Exiv2::ImageFactory::open(an_image_filename)};

        if (output_meta.get() && output_meta->good())
        {
            assert(some_named_meta_begin != some_named_meta_end);

            named_meta_array::const_iterator
                desired_meta(std::find_if(some_named_meta_begin, some_named_meta_end,
                                          [](const Named& named){return named.is_desired();}));
            named_meta_array::const_iterator safe_meta(desired_meta == some_named_meta_end ?
                                                       some_named_meta_begin :
                                                       desired_meta);

            copy(output_meta.get(), safe_meta->meta());
            augment(output_meta.get());
            {
                regexp_list blacklist;
                compile_blacklist(std::back_inserter(blacklist));
                enforce_blacklist(output_meta.get(), blacklist.begin(), blacklist.end());
            }
            output_meta->writeMetadata();

            if (desired_meta == some_named_meta_end)
            {
                std::stringstream message;
                message <<
                    "could not use metadata of selected image \"" << desired_meta->filename() <<
                    "\" - have used \"" << safe_meta->filename() << "\" instead";
                throw Warning(message.str());
            }

            return safe_meta;
        }
        else
        {
            throw Warning("initial output metadata is not valid");
        }
    }
#endif // HAVE_EXIV2
} // namespace metadata
