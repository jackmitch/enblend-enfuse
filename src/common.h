/*
 * Copyright (C) 2004-2007 Andrew Mihal
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
#ifndef __COMMON_H__
#define __COMMON_H__

#include <map>
#include <stdexcept>

#include <boost/assign/list_inserter.hpp>
#include <boost/assign/list_of.hpp>

#include "vigra/numerictraits.hxx"


// Defines to control how many -v flags are required for each type
// of message to be produced on stdout.
#define VERBOSE_ASSEMBLE_MESSAGES           0
#define VERBOSE_ABB_MESSAGES                1
#define VERBOSE_UBB_MESSAGES                1
#define VERBOSE_IBB_MESSAGES                1
#define VERBOSE_BLEND_MESSAGES              0
#define VERBOSE_NUMLEVELS_MESSAGES          0
#define VERBOSE_ROIBB_SIZE_MESSAGES         1
#define VERBOSE_MEMORY_ESTIMATION_MESSAGES  1
#define VERBOSE_CHECKPOINTING_MESSAGES      0
#define VERBOSE_INPUT_IMAGE_INFO_MESSAGES   1
#define VERBOSE_INPUT_UNION_SIZE_MESSAGES   1
#define VERBOSE_COLOR_CONVERSION_MESSAGES   0
#define VERBOSE_NFT_MESSAGES                0
#define VERBOSE_MASK_MESSAGES               0
#define VERBOSE_PYRAMID_MESSAGES            0
#define VERBOSE_CFI_MESSAGES                2
#define VERBOSE_GDA_MESSAGES                2
#define VERBOSE_VERSION_REPORTING           1

// Select our preferred type of image depending on what ./configure
// tells us.
#ifdef ENBLEND_CACHE_IMAGES
#define IMAGETYPE CachedFileImage
#else
#define IMAGETYPE BasicImage
#endif

namespace enblend {

/** The different image overlap classifications. */
enum Overlap {NoOverlap, PartialOverlap, CompleteOverlap};


/** Answer aString converted to uppercase letters. */
std::string
toUppercase(const std::string& aString)
{
    std::string result(aString);
    std::transform(aString.begin(), aString.end(), result.begin(), toupper);
    return result;
}


/** Answer aString converted to lowercase letters. */
std::string
toLowercase(const std::string& aString)
{
    std::string result(aString);
    std::transform(aString.begin(), aString.end(), result.begin(), tolower);
    return result;
}


/** Answer the VIGRA file type as determined by the extension of
 *  aFileName. */
std::string
getFileType(const std::string& aFileName)
{
    const std::string ext =
        toUppercase(aFileName.substr(aFileName.rfind(".") + 1));

    if (ext == "JPG") return "JPEG";
    else if (ext == "TIF") return "TIFF";
    else if (ext == "VIF") return "VIFF";
    else if (ext == "PBM" || ext == "PGM" || ext == "PPM") return "PNM";
    else return ext;
}


/** Convert an anOutputDepth to a "pixel type" string understood by
 * VIGRA. */
std::string
outputPixelTypeOfString(const char* anOutputDepth)
{
    typedef std::map<std::string, std::string> Str2StrMapType;
    Str2StrMapType depthMap;

    boost::assign::insert(depthMap)
        ("INT16", "INT16")
        ("INT32", "INT32")

        ("8", "UINT8")
        ("16", "UINT16")
        ("32", "UINT32")
        ("UINT8", "UINT8")
        ("UINT16", "UINT16")
        ("UINT32", "UINT32")

        ("DOUBLE", "DOUBLE")
        ("FLOAT", "FLOAT")
        ("R32", "FLOAT")
        ("R64", "DOUBLE")
        ("REAL32", "FLOAT")
        ("REAL64", "DOUBLE");

    Str2StrMapType::const_iterator p = depthMap.find(toUppercase(anOutputDepth));
    if (p == depthMap.end())
    {
        throw std::invalid_argument(std::string("unknown output depth \"") +
                                    anOutputDepth + "\"");
    }
    else
    {
        return p->second;
    }
}


/** Answer the best pixel type of an image given aFileType with
 * respect to aPixelType.  This is the type with the largest range. */
std::string
bestPixelType(const std::string& aFileType, const std::string& aPixelType)
{
    if (aFileType == "BMP" || aFileType == "JPEG" || aFileType == "RAS")
        return "UINT8";
    else if (aFileType == "PNG" &&
             (aPixelType == "INT32" || aPixelType == "UINT32" ||
              aPixelType == "FLOAT" || aPixelType == "DOUBLE"))
        return "UINT16";
    else if (aFileType == "EXR")
        return "FLOAT";
    else
        return aPixelType;
}


typedef std::pair<double, double> range_t;


/** Answer the maximum range of values aPixelType can represent. */
range_t
rangeOfPixelType(const std::string& aPixelType)
{
    typedef std::map<std::string, range_t> Str2PairMapType;
    Str2PairMapType rangeMap;

    boost::assign::insert(rangeMap)
        ("INT8", std::make_pair(vigra::NumericTraits<vigra::Int8>::min(),
                                vigra::NumericTraits<vigra::Int8>::max()))
        ("INT16", std::make_pair(vigra::NumericTraits<vigra::Int16>::min(),
                                 vigra::NumericTraits<vigra::Int16>::max()))
        ("INT32", std::make_pair(vigra::NumericTraits<vigra::Int32>::min(),
                                 vigra::NumericTraits<vigra::Int32>::max()))

        ("UINT8", std::make_pair(0.0, vigra::NumericTraits<vigra::UInt8>::max()))
        ("UINT16", std::make_pair(0.0, vigra::NumericTraits<vigra::UInt16>::max()))
        ("UINT32", std::make_pair(0.0, vigra::NumericTraits<vigra::UInt32>::max()))

        ("FLOAT", std::make_pair(0.0, 1.0))
        ("DOUBLE", std::make_pair(0.0, 1.0));

    assert(!aPixelType.empty());
    Str2PairMapType::const_iterator r = rangeMap.find(aPixelType);
    if (r == rangeMap.end())
    {
        throw std::invalid_argument(std::string("unknown pixel type \"") + aPixelType + "\"");
    }
    else
    {
        return r->second;
    }
}


/** Answer whether aPixelType defines a range that is so larges that
 *  it includes both aRange and anotherRange. */
bool
includesBothRanges(const std::string& aPixelType,
                   const range_t& aRange,
                   const range_t& anotherRange)
{
    const range_t range = rangeOfPixelType(aPixelType);

    return (aRange.first >= range.first && aRange.second <= range.second &&
            anotherRange.first >= range.first && anotherRange.second <= range.second);
}


/** Answer the smallest pixel type that is larger or equal to both
 *  aPixelType and anotherPixelType. */
std::string
maxPixelType(const std::string& aPixelType, const std::string& anotherPixelType)
{
    const range_t range1 = rangeOfPixelType(aPixelType);
    const range_t range2 = rangeOfPixelType(anotherPixelType);

    if (aPixelType == "DOUBLE" || anotherPixelType == "DOUBLE") {
        return "DOUBLE";
    } else if (aPixelType == "FLOAT" || anotherPixelType == "FLOAT") {
        return "FLOAT";
    } else if (range1.first <= range2.first && range1.second >= range2.second) {
        return aPixelType;      // first includes second
    } else if (range2.first <= range1.first && range2.second >= range1.second) {
        return anotherPixelType; // second includes first
    } else {
        // Types are different: look for the smallest containing type
        using namespace boost::assign;

        typedef std::vector<std::string> string_array;
        typedef string_array::const_iterator string_array_ci;

        if (range1.first < 0 || range2.first < 0) {
            const string_array types = list_of("INT8")("INT16")("INT32");
            for (string_array_ci i = types.begin(); i != types.end(); ++i) {
                if (includesBothRanges(*i, range1, range2)) {
                    return *i;
                }
            }
            return "INT32";
        } else {
            const string_array types = list_of("UINT8")("UINT16")("UINT32");
            for (string_array_ci i = types.begin(); i != types.end(); ++i) {
                if (includesBothRanges(*i, range1, range2)) {
                    return *i;
                }
            }
            return "UINT32";
        }
    }
}


} // namespace enblend

#endif /* __COMMON_H__ */

// Local Variables:
// mode: c++
// End:
