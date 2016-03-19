/*
 * Copyright (C) 2009-2016 Dr. Christoph L. Spiel
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
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

// Here we define macros and types that we already need in the
// definitions of global variables.

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cerrno>
#include <sstream>
#include <stdexcept>
#include <locale>

#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
#else
#include <map>
#endif

#include <vigra/numerictraits.hxx>


// Defines to control how many -v flags are required for each type
// of message to be produced on stdout.
#define VERBOSE_ASSEMBLE_MESSAGES           1
#define VERBOSE_CHECKPOINTING_MESSAGES      1
#define VERBOSE_OPENCL_MESSAGES             1 // ANTICIPATED CHANGE: raise level when we are done with debugging

#define VERBOSE_BLEND_MESSAGES              2
#define VERBOSE_MASK_MESSAGES               2
#define VERBOSE_NFT_MESSAGES                2
#define VERBOSE_PYRAMID_MESSAGES            2
#define VERBOSE_SIGNATURE_REPORTING         2
#define VERBOSE_TIFF_MESSAGES               2
#define VERBOSE_VERSION_REPORTING           2

#define VERBOSE_COLOR_CONVERSION_MESSAGES   3
#define VERBOSE_DIFFERENCE_STATISTICS       3
#define VERBOSE_LAYER_SELECTION             3
#define VERBOSE_RESPONSE_FILES              3

#define VERBOSE_ABB_MESSAGES                4
#define VERBOSE_IBB_MESSAGES                4
#define VERBOSE_INPUT_IMAGE_INFO_MESSAGES   4
#define VERBOSE_INPUT_UNION_SIZE_MESSAGES   4
#define VERBOSE_ROIBB_SIZE_MESSAGES         4
#define VERBOSE_UBB_MESSAGES                4

#define VERBOSE_CFI_MESSAGES                5
#define VERBOSE_GDA_MESSAGES                5

#define VERBOSE_MEMORY_ESTIMATION_MESSAGES  6


//< default-output-filename a.tif
#define DEFAULT_OUTPUT_FILENAME "a.tif"


namespace enblend
{
    inline static void
    to_lower(std::string& a_string)
    {
        for (auto& c : a_string)
        {
            c = std::tolower(c, std::locale());
        }
    }


    inline static void
    to_upper(std::string& a_string)
    {
        for (auto& c : a_string)
        {
            c = std::toupper(c, std::locale());
        }
    }


    inline static std::string
    to_lower_copy(const std::string& a_string)
    {
        std::string result;
        std::for_each(a_string.begin(), a_string.end(),
                      [&](char c) {result.push_back(std::tolower(c, std::locale()));});
        return result;
    }


    inline static std::string
    to_upper_copy(const std::string& a_string)
    {
        std::string result;
        std::for_each(a_string.begin(), a_string.end(),
                      [&](char c) {result.push_back(std::toupper(c, std::locale()));});
        return result;
    }


    inline static void
    trim(std::string& a_string)
    {
        const size_t start = a_string.find_first_not_of(" \t\n\r");
        const size_t end = a_string.find_last_not_of(" \t\n\r");

        a_string.assign(a_string, start, 1U + end - start);
    }
} // end namespace enblend


/** The different kinds of boundary conditions we can impose upon an
 *  image. */
typedef enum BoundaryKind
{
    UnknownWrapAround,          // unknown kind
    OpenBoundaries,             // contractible
    HorizontalStrip,            // contractible along 2nd axis
    VerticalStrip,              // contractible along 1st axis
    DoubleStrip                 // non-contractible
} boundary_t;

enum MainAlgo {
    NFT, GraphCut
};


// Colorspaces available for pyramidal blending operations
typedef enum
{
    UndeterminedColorspace,    // explicit `not-a-colorspace' value
    IdentitySpace,             // (1d) luminance interval for grayscale images and (3d) RGB-cube for RGB-images
    CIELAB,
    CIELUV,
    CIECAM
} blend_colorspace_t;


//< default-tiff-resolution 300
#define DEFAULT_TIFF_RESOLUTION 300.0f // units are dots-per-inch ("DPI")


struct TiffResolution {
    TiffResolution() : x(0.0f), y(0.0f) {}

    TiffResolution(float anXresolution, float aYresolution) :
        x(anXresolution), y(aYresolution) {}

    bool operator==(const TiffResolution& anOther) const {
        return this->x == anOther.x && this->y == anOther.y;
    }

    bool operator!=(const TiffResolution& anOther) const {
        return !operator==(anOther);
    }

    float x;
    float y;
};


#endif /* __GLOBAL_H__ */

// Local Variables:
// mode: c++
// End:
