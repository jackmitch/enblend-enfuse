/*
 * Copyright (C) 2004-2009 Andrew Mihal
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

#include <iomanip>
#include <map>
#include <stdexcept>

#include <boost/assign/list_inserter.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/scoped_ptr.hpp>

#include "vigra/numerictraits.hxx"

#include "filenameparse.h"

#define NUMERIC_OPTION_DELIMITERS ";:/"
#define PATH_OPTION_DELIMITERS ",;:"

#define MASK_COMPRESSION "DEFLATE"

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


// Colors used in the optimizer visualization
#define VISUALIZE_RGB_COLOR_BLUE1    RGBValue<vigra::UInt8>(  0,   0, 255)
#define VISUALIZE_RGB_COLOR_BLUE2    RGBValue<vigra::UInt8>(  0,   0, 238)
#define VISUALIZE_RGB_COLOR_BLUE3    RGBValue<vigra::UInt8>(  0,   0, 205)
#define VISUALIZE_RGB_COLOR_BLUE4    RGBValue<vigra::UInt8>(  0,   0, 139)
#define VISUALIZE_RGB_COLOR_CYAN1    RGBValue<vigra::UInt8>(  0, 255, 255)
#define VISUALIZE_RGB_COLOR_CYAN2    RGBValue<vigra::UInt8>(  0, 238, 238)
#define VISUALIZE_RGB_COLOR_CYAN3    RGBValue<vigra::UInt8>(  0, 205, 205)
#define VISUALIZE_RGB_COLOR_CYAN4    RGBValue<vigra::UInt8>(  0, 139, 139)
#define VISUALIZE_RGB_COLOR_GREEN1   RGBValue<vigra::UInt8>(  0, 255,   0)
#define VISUALIZE_RGB_COLOR_GREEN2   RGBValue<vigra::UInt8>(  0, 238,   0)
#define VISUALIZE_RGB_COLOR_GREEN3   RGBValue<vigra::UInt8>(  0, 205,   0)
#define VISUALIZE_RGB_COLOR_GREEN4   RGBValue<vigra::UInt8>(  0, 139,   0)
#define VISUALIZE_RGB_COLOR_MAGENTA1 RGBValue<vigra::UInt8>(255,   0, 255)
#define VISUALIZE_RGB_COLOR_MAGENTA2 RGBValue<vigra::UInt8>(238,   0, 238)
#define VISUALIZE_RGB_COLOR_MAGENTA3 RGBValue<vigra::UInt8>(205,   0, 205)
#define VISUALIZE_RGB_COLOR_MAGENTA4 RGBValue<vigra::UInt8>(139,   0, 139)
#define VISUALIZE_RGB_COLOR_RED1     RGBValue<vigra::UInt8>(255,   0,   0)
#define VISUALIZE_RGB_COLOR_RED2     RGBValue<vigra::UInt8>(238,   0,   0)
#define VISUALIZE_RGB_COLOR_RED3     RGBValue<vigra::UInt8>(205,   0,   0)
#define VISUALIZE_RGB_COLOR_RED4     RGBValue<vigra::UInt8>(139,   0,   0)
#define VISUALIZE_RGB_COLOR_YELLOW1  RGBValue<vigra::UInt8>(255, 255,   0)
#define VISUALIZE_RGB_COLOR_YELLOW2  RGBValue<vigra::UInt8>(238, 238,   0)
#define VISUALIZE_RGB_COLOR_YELLOW3  RGBValue<vigra::UInt8>(205, 205,   0)
#define VISUALIZE_RGB_COLOR_YELLOW4  RGBValue<vigra::UInt8>(139, 139,   0)

#define VISUALIZE_SHORT_PATH_VALUE   VISUALIZE_RGB_COLOR_YELLOW1
#define VISUALIZE_FIRST_VERTEX_VALUE VISUALIZE_RGB_COLOR_GREEN3
#define VISUALIZE_NEXT_VERTEX_VALUE  VISUALIZE_RGB_COLOR_GREEN2
#define VISUALIZE_NO_OVERLAP_VALUE   VISUALIZE_RGB_COLOR_RED4
#define VISUALIZE_STATE_SPACE        VISUALIZE_RGB_COLOR_BLUE3
#define VISUALIZE_STATE_SPACE_INSIDE VISUALIZE_RGB_COLOR_CYAN1
#define VISUALIZE_STATE_SPACE_UNCONVERGED VISUALIZE_RGB_COLOR_MAGENTA1


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


/** String tokenizer similar to strtok_r().
 *  In contrast to strtok_r this function returns an empty string for
 *  each pair of successive delimiters.  Function strtok_r skips them.
 */
char*
strtoken_r(char *str, const char *delim, char **save_ptr)
{
    char *s = str == NULL ? *save_ptr : str;
    if (s == NULL) return NULL;
    else
    {
        char *token = s;
        while (*s != 0 && strchr(delim, (int) *s) == NULL) s++;
        *save_ptr = *s == 0 ? NULL : s + 1;
        *s = 0;
        return token;
    }
}


/** Answer the error message associated with anErrorNumber.
 */
std::string
errorMessage(int anErrorNumber)
{
#if HAVE_STRERROR_R || HAVE_STRERROR
    const size_t size = 256;
    boost::scoped_ptr<char> message(new char[size]);
#if HAVE_STRERROR_R
    strerror_r(anErrorNumber, message.get(), size);
#elif HAVE_STRERROR
    strncpy(message.get(), strerror(anErrorNumber), size);
#endif
    return std::string(message.get());
#else
    std::ostringstream oss;
    oss << "no message available, error: " << anErrorNumber;
    return oss.str();
#endif
}


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


/** Compute the integral logarithm of n to the base 10.  We do not
 *  need to take special care of the case n == 0 for our purposes. */
unsigned
ilog10(unsigned n)
{
    return n <= 9 ? 0 : 1 + ilog10(n / 10);
}


/** Expand aTemplate filling the variable parts with anInputFilename,
 *  anOutputFilename, and aNumber depending on the conversion
 *  specifiers in aTemplate.
 *
 *  Conversion Specifiers - lowercase characters refer to
 *  anInputFilename whereas uppercase ones refer to anOutputFilename:
 *      %%        A single '%'-sign
 *      %i        aNumber unaltered
 *      %n        successor of aNumber
 *      %p        aFilename unaltered
 *      %d        directory part of aFilename
 *      %b        non-directory part (aka basename) of aFilename
 *      %f        basename of aFilename without extension
 *      %e        extension of aFilename (including the leading dot)
 *  All other characters in aTemplate are passed through literally.
 *
 *  The "%i" and "%n" conversions honor a flag which is either
 *      '0'       pad with zeros (default) or
 *      PUNCT     i.e. any punctuation character to pad with
 *  and a width specification.  If no width is requested, the
 *  width is computed based on aNumberOfImages.
 *
 *  For example
 *          expandFilenameTemplate("mask-%04n.tif", 2, "foobar.jpg", 9)
 *  evaluates to
 *          mask-0009.tif
 */
std::string
expandFilenameTemplate(const std::string& aTemplate,
                       unsigned aNumberOfImages,
                       const std::string& anInputFilename,
                       const std::string& anOutputFilename,
                       unsigned aNumber)
{
    std::string result;

    for (std::string::const_iterator c = aTemplate.begin();
         c != aTemplate.end();
         ++c)
    {
        if (*c == '%')
        {
            ++c;
            if (c == aTemplate.end())
            {
                result.push_back(*c);
            }
            else
            {
                char pad = 0;
                while (c != aTemplate.end() && (*c == '0' || ispunct(*c)))
                {
                    pad = *c;
                    ++c;
                }

                std::string width;
                while (c != aTemplate.end() && isdigit(*c))
                {
                    width.push_back(*c);
                    ++c;
                }

                if (c != aTemplate.end())
                {
                    switch (*c)
                    {
                    case '%':
                        result.push_back(*c);
                        break;
                    case 'n':
                        ++aNumber;
                        ++aNumberOfImages;
                        // FALLTHROUGH
                    case 'i':
                    {
                        std::ostringstream oss;
                        oss <<
                            std::setw(width.empty() ? 1 + ilog10(aNumberOfImages - 1) : atoi(width.c_str())) <<
                            std::setfill(pad == 0 ? '0' : pad) <<
                            aNumber;
                        result.append(oss.str());
                        break;
                    }
                    case 'P':
                        result.append(anOutputFilename);
                        break;
                    case 'p':
                        result.append(anInputFilename);
                        break;
                    case 'D':
                        result.append(extractDirname(anOutputFilename));
                        break;
                    case 'd':
                        result.append(extractDirname(anInputFilename));
                        break;
                    case 'B':
                        result.append(extractBasename(anOutputFilename));
                        break;
                    case 'b':
                        result.append(extractBasename(anInputFilename));
                        break;
                    case 'F':
                        result.append(extractFilename(anOutputFilename));
                        break;
                    case 'f':
                        result.append(extractFilename(anInputFilename));
                        break;
                    case 'E':
                        result.append(extractExtension(anOutputFilename));
                        break;
                    case 'e':
                        result.append(extractExtension(anInputFilename));
                        break;
                    default:
                        std::cerr <<
                            command <<
                            ": warning: ignoring unknown variable character ";
                        if (isprint(*c))
                        {
                            std::cerr << "'" << *c << "'";
                        }
                        else
                        {
                            std::cerr << "0x" << std::hex << *c;
                        }
                        std::cerr << " in\n"
                                  << command
                                  << ": warning:     filename template \""
                                  << aTemplate
                                  << "\""
                                  << std::endl;
                    } // switch (*c)
                }
            }
        }
        else
        {
            result.push_back(*c);
        }
    }

    return result;
}


} // namespace enblend

#endif /* __COMMON_H__ */

// Local Variables:
// mode: c++
// End:
