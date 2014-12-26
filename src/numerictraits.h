/*
 * Copyright (C) 2004-2014 Andrew Mihal
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
#ifndef NUMERICTRAITS_H_INCLUDED_
#define NUMERICTRAITS_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vigra/basicimage.hxx>
#include <vigra/rgbvalue.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/utilities.hxx>

#include "common.h"


namespace enblend
{
    struct Error_EnblendNumericTraits_not_specialized_for_this_case {};


    template <class A>
    struct EnblendNumericTraits
    {
        // Types related to input images
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImagePixelComponentType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImagePixelType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImageType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImageIsScalar;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case AlphaPixelType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case AlphaType;

        // Types related to the mask
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskPixelType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskType;

        // Types related to image pyramids
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImagePyramidPixelComponentType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImagePyramidPixelType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case ImagePyramidType;
        enum {ImagePyramidIntegerBits = 0};
        enum {ImagePyramidFractionBits = 0};

        // Pixel type used by SKIPSM algorithm for intermediate image pixel calculations
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMImagePixelType;

        // Pixel type used by SKIPSM algorithm for intermediate alpha pixel calculations
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMAlphaPixelType;

        // Types related to mask pyramid
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskPyramidPixelType;
        typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskPyramidType;
        enum {MaskPyramidIntegerBits = 0};
        enum {MaskPyramidFractionBits = 0};

        typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMMaskPixelType;
    };


#define ENBLEND_NUMERICTRAITS(IMAGE, IMAGECOMPONENT, ALPHA, MASK, \
                              PYRAMIDCOMPONENT, PYRAMIDINTEGER, PYRAMIDFRACTION, \
                              SKIPSMIMAGE, SKIPSMALPHA, \
                              MASKPYRAMID, MASKPYRAMIDINTEGER, MASKPYRAMIDFRACTION, SKIPSMMASK) \
    template <>                                                         \
    struct EnblendNumericTraits<IMAGECOMPONENT>                         \
    {                                                                   \
        typedef IMAGECOMPONENT ImagePixelComponentType;                 \
        typedef IMAGECOMPONENT ImagePixelType;                          \
        typedef IMAGE<IMAGECOMPONENT> ImageType;                        \
        typedef vigra::VigraTrueType ImageIsScalar;                     \
        typedef ALPHA AlphaPixelType;                                   \
        typedef IMAGE<ALPHA> AlphaType;                                 \
        typedef MASK MaskPixelType;                                     \
        typedef IMAGE<MASK> MaskType;                                   \
        typedef PYRAMIDCOMPONENT ImagePyramidPixelComponentType;        \
        typedef PYRAMIDCOMPONENT ImagePyramidPixelType;                 \
        typedef IMAGE<PYRAMIDCOMPONENT> ImagePyramidType;               \
        enum {ImagePyramidIntegerBits = PYRAMIDINTEGER};                \
        enum {ImagePyramidFractionBits = PYRAMIDFRACTION};              \
        typedef SKIPSMIMAGE SKIPSMImagePixelComponentType;              \
        typedef SKIPSMIMAGE SKIPSMImagePixelType;                       \
        static_assert(ImagePyramidIntegerBits >= 0, "ImagePyramidIntegerBits must be non-negative"); \
        static_assert(ImagePyramidFractionBits >= 0, "ImagePyramidFractionBits must be non-negative"); \
        static_assert(ImagePyramidIntegerBits + ImagePyramidFractionBits + 6 <= 8 * sizeof(SKIPSMImagePixelType) - 1, \
                      "ImagePyramidIntegerBits + ImagePyramidFractionBits do not fit into SKIPSMImagePixelType"); \
        typedef SKIPSMALPHA SKIPSMAlphaPixelType;                       \
        typedef MASKPYRAMID MaskPyramidPixelType;                       \
        typedef IMAGE<MASKPYRAMID> MaskPyramidType;                     \
        enum {MaskPyramidIntegerBits = MASKPYRAMIDINTEGER};             \
        enum {MaskPyramidFractionBits = MASKPYRAMIDFRACTION};           \
        typedef SKIPSMMASK SKIPSMMaskPixelType;                         \
        static_assert(MaskPyramidIntegerBits >= 0, "MaskPyramidIntegerBits must be non-negative"); \
        static_assert(MaskPyramidFractionBits >= 0, "MaskPyramidFractionBits must be non-negative"); \
        static_assert(MaskPyramidIntegerBits + MaskPyramidFractionBits + 6 <= 8 * sizeof(SKIPSMMaskPixelType) - 1, \
                      "MaskPyramidIntegerBits + MaskPyramidFractionBits do not fit into SKIPSMMaskPixelType"); \
    };                                                                  \
                                                                        \
    template <>                                                         \
    struct EnblendNumericTraits<vigra::RGBValue<IMAGECOMPONENT, 0, 1, 2> > \
    {                                                                   \
        typedef IMAGECOMPONENT ImagePixelComponentType;                 \
        typedef vigra::RGBValue<IMAGECOMPONENT, 0, 1, 2> ImagePixelType; \
        typedef IMAGE<vigra::RGBValue<IMAGECOMPONENT, 0, 1, 2> > ImageType; \
        typedef vigra::VigraFalseType ImageIsScalar;                    \
        typedef ALPHA AlphaPixelType;                                   \
        typedef IMAGE<ALPHA> AlphaType;                                 \
        typedef MASK MaskPixelType;                                     \
        typedef IMAGE<MASK> MaskType;                                   \
        typedef PYRAMIDCOMPONENT ImagePyramidPixelComponentType;        \
        typedef vigra::RGBValue<PYRAMIDCOMPONENT, 0, 1, 2> ImagePyramidPixelType; \
        typedef IMAGE<vigra::RGBValue<PYRAMIDCOMPONENT, 0, 1, 2> > ImagePyramidType; \
        enum {ImagePyramidIntegerBits = PYRAMIDINTEGER};                \
        enum {ImagePyramidFractionBits = PYRAMIDFRACTION};              \
        typedef SKIPSMIMAGE SKIPSMImagePixelComponentType;              \
        typedef vigra::RGBValue<SKIPSMIMAGE, 0, 1, 2> SKIPSMImagePixelType; \
        typedef SKIPSMALPHA SKIPSMAlphaPixelType;                       \
        typedef MASKPYRAMID MaskPyramidPixelType;                       \
        typedef IMAGE<MASKPYRAMID> MaskPyramidType;                     \
        enum {MaskPyramidIntegerBits = MASKPYRAMIDINTEGER};             \
        enum {MaskPyramidFractionBits = MASKPYRAMIDFRACTION};           \
        typedef SKIPSMMASK SKIPSMMaskPixelType;                         \
    }


    // Traits for converting between image pixel types and pyramid
    // pixel types Pyramids require one more bit of precision than the
    // regular image type.
    //
    // Notes
    //   * SKIPSM math requires 6 bits on top of the pyramid type.
    //   * The maximum of an N-bit signed integer on a
    //     twos-complement machine is 2^(N - 1) - 1, which is why we
    //     dutifully subtract one of each inequalities' RHS.
    //   * Image Pyramids
    //     (ImagePyramidIntegerBits + 1) + ImagePyramidFractionBits + 6  <=  sizeof(SKIPSMImagePixelType) - 1
    //          (8 + 1) +  7 + 6  =  22  <=  32 - 1
    //         (16 + 1) +  7 + 6  =  30  <=  32 - 1
    //         (32 + 1) + 23 + 6  =  62  <=  64 - 1
    //             8    +  0 + 6  =  14      # fake value for floating-point scales
    //   * Mask Pyramids
    //     (MaskPyramidIntegerBits + 1) + MaskPyramidFractionBits + 6  <=  sizeof(SKIPSMMaskPixelType) - 1
    //          (8 + 1) +  7 + 6  =  22  <=  32 - 1
    //          (8 + 1) + 15 + 6  =  30  <=  32 - 1
    //
    //
    //                                    IMAGE-            ALPHA          MASK         PYRAMID-      IMG-PYR.      SKIPSM-         SKIPSM-          MASK-       MASK-PYR.     SKIPSM-
    //                                   COMPONENT                                      COMPONENT    INT   FRAC      IMAGE           ALPHA          PYRAMID     INT   FRAC      MASK
    //                                 ==============   =============  =============  =============  ====  ====  =============   =============   =============  ====  ====  ============
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::Int8,     vigra::UInt8,  vigra::UInt8,  vigra::Int16,    9,    7,  vigra::Int32,   vigra::Int16,   vigra::Int16,    9,    7,  vigra::Int32);
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::UInt8,    vigra::UInt8,  vigra::UInt8,  vigra::Int16,    9,    7,  vigra::Int32,   vigra::Int16,   vigra::Int16,    9,    7,  vigra::Int32);

    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::Int16,    vigra::UInt8,  vigra::UInt8,  vigra::Int32,   17,    7,  vigra::Int32,   vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::UInt16,   vigra::UInt8,  vigra::UInt8,  vigra::Int32,   17,    7,  vigra::Int32,   vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);

#ifdef PREFER_DOUBLE_TO_INT64_AS_SKIPSM_IMAGE_TYPE
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::Int32,    vigra::UInt8,  vigra::UInt8,  double,          8,    0,  double,         vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::UInt32,   vigra::UInt8,  vigra::UInt8,  double,          8,    0,  double,         vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
#else
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::Int32,    vigra::UInt8,  vigra::UInt8,  vigra::Int64,   33,   23,  vigra::Int64,   vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::UInt32,   vigra::UInt8,  vigra::UInt8,  vigra::Int64,   33,   23,  vigra::Int64,   vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
#endif

    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::Int64,    vigra::UInt8,  vigra::UInt8,  double,          8,    0,  double,         vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   vigra::UInt64,   vigra::UInt8,  vigra::UInt8,  double,          8,    0,  double,         vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);

    ENBLEND_NUMERICTRAITS(IMAGETYPE,   float,           vigra::UInt8,  vigra::UInt8,  double,          8,    0,  double,         vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
#ifdef PREFER_LONG_DOUBLE_TO_DOUBLE_AS_SKIPSM_IMAGE_TYPE
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   double,          vigra::UInt8,  vigra::UInt8,  long double,     8,    0,  long double,    vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
#else
    ENBLEND_NUMERICTRAITS(IMAGETYPE,   double,          vigra::UInt8,  vigra::UInt8,  double,          8,    0,  double,         vigra::Int16,   vigra::Int32,    9,   15,  vigra::Int32);
#endif

#undef ENBLEND_NUMERICTRAITS


    // Traits for correctly handling alpha/mask values in floating point files.
    // This differs from NumericTraits for floating point values.
#define ALPHA_TRAITS(T1, S)                     \
    template <>                                 \
    struct AlphaTraits<T1>                      \
    {                                           \
        static T1 max()                         \
        {                                       \
            return S;                           \
        }                                       \
                                                \
        static T1 zero()                        \
        {                                       \
            return 0;                           \
        }                                       \
    };                                          \
                                                \
    template <>                                 \
    struct AlphaTraits<vigra::RGBValue<T1> >    \
    {                                           \
        static T1 max()                         \
        {                                       \
            return S;                           \
        }                                       \
                                                \
        static T1 zero()                        \
        {                                       \
            return 0;                           \
        }                                       \
    }


    template <class T1> struct AlphaTraits;


    // Integral: 8 bits
    ALPHA_TRAITS(unsigned char, std::numeric_limits<unsigned char>::max());
    ALPHA_TRAITS(signed char, std::numeric_limits<signed char>::max());

    // Integral: 16 bits
    ALPHA_TRAITS(unsigned short, std::numeric_limits<unsigned short>::max());
    ALPHA_TRAITS(signed short, std::numeric_limits<signed short>::max());

    // Integral: 32 bits
    ALPHA_TRAITS(unsigned int, std::numeric_limits<unsigned int>::max());
    ALPHA_TRAITS(signed int, std::numeric_limits<signed int>::max());

    // Integral: 64 bits
    ALPHA_TRAITS(unsigned long long int, std::numeric_limits<unsigned long long int>::max());
    ALPHA_TRAITS(signed long long int, std::numeric_limits<signed long long int>::max());

    // Floating Point: 32/64/80 bits
    ALPHA_TRAITS(float, 1.0);
    ALPHA_TRAITS(double, 1.0);
    ALPHA_TRAITS(long double, 1.0);

#undef ALPHA_TRAITS
} // namespace enblend


#endif // NUMERICTRAITS_H_INCLUDED_


// Local Variables:
// mode: c++
// End:
