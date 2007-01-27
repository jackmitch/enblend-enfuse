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
#ifndef __NUMERICTRAITS_H__
#define __NUMERICTRAITS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vigra/basicimage.hxx"
#include "vigra/cachedfileimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/utilities.hxx"

using vigra::BasicImage;
using vigra::CachedFileImage;
using vigra::NumericTraits;
using vigra::RGBValue;
using vigra::VigraFalseType;
using vigra::VigraTrueType;

using vigra::Int8;
using vigra::Int16;
using vigra::Int32;
using vigra::Int64;
using vigra::UInt8;
using vigra::UInt16;
using vigra::UInt32;
using vigra::UInt64;

#ifdef ENBLEND_CACHE_IMAGES
    #define IMAGETYPE CachedFileImage
#else
    #define IMAGETYPE BasicImage
#endif

namespace enblend {

struct Error_EnblendNumericTraits_not_specialized_for_this_case { };

template<class A>
struct EnblendNumericTraits {
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
    enum { ImagePyramidIntegerBits = 0 };
    enum { ImagePyramidFractionBits = 0 };

    // Pixel type used by SKIPSM algorithm for intermediate image pixel calculations
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMImagePixelType;

    // Pixel type used by SKIPSM algorithm for intermediate alpha pixel calculations
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMAlphaPixelType;

    // Types related to mask pyramid
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskPyramidPixelType;
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskPyramidType;
    enum { MaskPyramidIntegerBits = 0 };
    enum { MaskPyramidFractionBits = 0 };

    typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMMaskPixelType;
};

#define DEFINE_ENBLENDNUMERICTRAITS(IMAGE, IMAGECOMPONENT, ALPHA, MASK, PYRAMIDCOMPONENT, PYRAMIDINTEGER, PYRAMIDFRACTION, SKIPSMIMAGE, SKIPSMALPHA, MASKPYRAMID, MASKPYRAMIDINTEGER, MASKPYRAMIDFRACTION, SKIPSMMASK) \
template<> \
struct EnblendNumericTraits<IMAGECOMPONENT> { \
    typedef IMAGECOMPONENT ImagePixelComponentType; \
    typedef IMAGECOMPONENT ImagePixelType; \
    typedef IMAGE<IMAGECOMPONENT> ImageType; \
    typedef VigraTrueType ImageIsScalar; \
    typedef ALPHA AlphaPixelType; \
    typedef IMAGE<ALPHA> AlphaType; \
    typedef MASK MaskPixelType; \
    typedef IMAGE<MASK> MaskType; \
    typedef PYRAMIDCOMPONENT ImagePyramidPixelComponentType; \
    typedef PYRAMIDCOMPONENT ImagePyramidPixelType; \
    typedef IMAGE<PYRAMIDCOMPONENT> ImagePyramidType; \
    enum {ImagePyramidIntegerBits = PYRAMIDINTEGER}; \
    enum {ImagePyramidFractionBits = PYRAMIDFRACTION}; \
    typedef SKIPSMIMAGE SKIPSMImagePixelType; \
    typedef SKIPSMALPHA SKIPSMAlphaPixelType; \
    typedef MASKPYRAMID MaskPyramidPixelType; \
    typedef IMAGE<MASKPYRAMID> MaskPyramidType; \
    enum {MaskPyramidIntegerBits = MASKPYRAMIDINTEGER}; \
    enum {MaskPyramidFractionBits = MASKPYRAMIDFRACTION}; \
    typedef SKIPSMMASK SKIPSMMaskPixelType; \
};\
template<> \
struct EnblendNumericTraits<RGBValue<IMAGECOMPONENT,0,1,2> > { \
    typedef IMAGECOMPONENT ImagePixelComponentType; \
    typedef RGBValue<IMAGECOMPONENT,0,1,2> ImagePixelType; \
    typedef IMAGE<RGBValue<IMAGECOMPONENT,0,1,2> > ImageType; \
    typedef VigraFalseType ImageIsScalar; \
    typedef ALPHA AlphaPixelType; \
    typedef IMAGE<ALPHA> AlphaType; \
    typedef MASK MaskPixelType; \
    typedef IMAGE<MASK> MaskType; \
    typedef PYRAMIDCOMPONENT ImagePyramidPixelComponentType; \
    typedef RGBValue<PYRAMIDCOMPONENT,0,1,2> ImagePyramidPixelType; \
    typedef IMAGE<RGBValue<PYRAMIDCOMPONENT,0,1,2> > ImagePyramidType; \
    enum {ImagePyramidIntegerBits = PYRAMIDINTEGER}; \
    enum {ImagePyramidFractionBits = PYRAMIDFRACTION}; \
    typedef RGBValue<SKIPSMIMAGE,0,1,2> SKIPSMImagePixelType; \
    typedef SKIPSMALPHA SKIPSMAlphaPixelType; \
    typedef MASKPYRAMID MaskPyramidPixelType; \
    typedef IMAGE<MASKPYRAMID> MaskPyramidType; \
    enum {MaskPyramidIntegerBits = MASKPYRAMIDINTEGER}; \
    enum {MaskPyramidFractionBits = MASKPYRAMIDFRACTION}; \
    typedef SKIPSMMASK SKIPSMMaskPixelType; \
};

// Traits for converting between image pixel types and pyramid pixel types
// Pyramids require one more bit of precision than the regular image type.
// SKIPSM math requires 6 bits on top of the pyramid type.
//                                     IMCOMP  ALPHA  MASK   IMPYR   I   F  SKIPI   SKIPA  MASKPYR I  F  SKIPM
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, Int8,   UInt8, UInt8, Int16,  9,  7, Int32,  Int16, Int16, 9,  7, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, UInt8,  UInt8, UInt8, Int16,  9,  7, Int32,  Int16, Int16, 9,  7, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, Int16,  UInt8, UInt8, Int32, 17,  9, Int32,  Int16, Int32, 9, 17, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, UInt16, UInt8, UInt8, Int32, 17,  7, Int32,  Int16, Int32, 9, 15, Int32)
//DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, Int32,  UInt8, UInt8, Int64, 33, 25, Int64,  Int16, Int32, 9, 17, Int32)
//DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, UInt32, UInt8, UInt8, Int64, 33, 25, Int64,  Int16, Int32, 9, 17, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, Int32,  UInt8, UInt8, double, 8,  0, double, Int16, Int32, 9, 17, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, UInt32, UInt8, UInt8, double, 8,  0, double, Int16, Int32, 9, 17, Int32)
//DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, Int64,  UInt8, UInt8, double, 8,  0, double, Int16, Int32, 9, 17, Int32)
//DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, UInt64, UInt8, UInt8, double, 8,  0, double, Int16, Int32, 9, 17, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, float,  UInt8, UInt8, double, 8,  0, double, Int16, Int32, 9, 17, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGETYPE, double, UInt8, UInt8, double, 8,  0, double, Int16, Int32, 9, 17, Int32)

} // namespace enblend

#endif /* __FIXMATH_H__ */
