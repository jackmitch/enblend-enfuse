/*
 * Copyright (C) 2004-2005 Andrew Mihal
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
#ifndef __FIXMATH_H__
#define __FIXMATH_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "vigra/basicimage.hxx"
#include "vigra/cachedfileimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/utilities.hxx"

using std::cout;
using std::pair;

using vigra::BasicImage;
using vigra::CachedFileImage;
using vigra::NumericTraits;
using vigra::RGBValue;
using vigra::triple;
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
    #define IMAGE CachedFileImage
#else
    #define IMAGE BasicImage
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
    enum { ImagePyramidIntegerBits = 0 };
    enum { ImagePyramidFractionBits = 0 };

    // Pixel type used by SKIPSM algorithm for intermediate image pixel calculations
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMImagePixelType;

    // Pixel type used by SKIPSM algorithm for intermediate alpha pixel calculations
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case SKIPSMAlphaPixelType;

    // Types related to mask pyramid
    typedef Error_EnblendNumericTraits_not_specialized_for_this_case MaskPyramidPixelType;
    enum { MaskPyramidIntegerBits = 0 };
    enum { MaskPyramidFractionBits = 0 };
};

#define DEFINE_ENBLENDNUMERICTRAITS(I, A, B, C, D, E) \
template<> \
struct EnblendNumericTraits<A> { \
    typedef A ImagePixelComponentType; \
    typedef A ImagePixelType; \
    typedef I<A> ImageType; \
    typedef B PyramidPixelComponentType; \
    typedef B PyramidPixelType; \
    enum {PyramidIntegerBits = C}; \
    enum {PyramidFractionBits = D}; \
    typedef I<B> PyramidType; \
    typedef E SKIPSMPixelComponentType; \
    typedef E SKIPSMPixelType; \
};\
template<> \
struct EnblendNumericTraits<RGBValue<A,0,1,2> > { \
    typedef A ImagePixelComponentType; \
    typedef RGBValue<A,0,1,2> ImagePixelType; \
    typedef I<RGBValue<A,0,1,2> > ImageType; \
    typedef B PyramidPixelComponentType; \
    typedef RGBValue<B,0,1,2> PyramidPixelType; \
    enum {PyramidIntegerBits = C}; \
    enum {PyramidFractionBits = D}; \
    typedef I<RGBValue<B,0,1,2> > PyramidType; \
    typedef E SKIPSMPixelComponentType; \
    typedef RGBValue<E,0,1,2> SKIPSMPixelType; \
};

// Traits for converting between image pixel types and pyramid pixel types
// Pyramids require one more bit of precision than the regular image type.
// SKIPSM math requires 6 bits on top of the pyramid type.
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, Int8, Int16, 9, 7, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, UInt8, Int16, 9, 7, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, Int16, Int32, 17, 9, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, UInt16, Int32, 17, 9, Int32)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, Int32, Int64, 33, 25, Int64)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, UInt32, Int64, 33, 25, Int64)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, Int64, double, 0, 0, double)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, UInt64, double, 0, 0, double)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, float, double, 0, 0, double)
DEFINE_ENBLENDNUMERICTRAITS(IMAGE, double, double, 0, 0, double)

/** A functor for converting scalar pixel values to the number representation used
 *  for pyramids. These are either fixed-point integers or floating-point numbers.
 */
template <typename SrcPixelType, typename PyramidPixelType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertScalarToPyramidFunctor {

public:
    ConvertScalarToPyramidFunctor() { }

    inline PyramidPixelType operator()(const SrcPixelType &v) const {
        return doConvert(v, SrcIsIntegral(), PyramidIsIntegral());
    }

protected:

    typedef typename NumericTraits<SrcPixelType>::isIntegral SrcIsIntegral;
    typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

    // Convert an integral pixel type to an integral pyramid value type.
    inline PyramidPixelType doConvert(const SrcPixelType &v, VigraTrueType, VigraTrueType) const {
        return convertIntegerToFixedPoint(v);
    }

    // Convert an integral pixel type to a real pyramid value type.
    inline PyramidPixelType doConvert(const SrcPixelType &v, VigraTrueType, VigraFalseType) const {
        return NumericTraits<SrcPixelType>::toRealPromote(v);
    }

    // Convert a real pixel type to an integral pyramid value type.
    inline PyramidPixelType doConvert(const SrcPixelType &v, VigraFalseType, VigraTrueType) const {
        return convertDoubleToFixedPoint(v);
    }

    // Convert a real pixel type to a real pyramid value type.
    inline PyramidPixelType doConvert(const SrcPixelType &v, VigraFalseType, VigraFalseType) const {
        return v;
    }

    inline PyramidPixelType convertDoubleToFixedPoint(const double &v) const {
        // Shift v to get the appropriate number of fraction bits into the integer part,
        // then fromRealPromote this value into the fixed-point type.
        return NumericTraits<PyramidPixelType>::fromRealPromote(v * (double)(1 << PyramidFractionBits));
    };

    inline PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType &v) const {
        // Shift v left to move the decimal point and set the fraction bits to zero.
        return (PyramidPixelType)v << PyramidFractionBits;
    };

};

/** A functor for converting numbers stored in the pyramid number representation back
 *  into normal pixel values.
 */
template <typename DestPixelType, typename PyramidPixelType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertPyramidToScalarFunctor {

public:
    ConvertPyramidToScalarFunctor() { }

    inline DestPixelType operator()(const PyramidPixelType &v) const {
        return doConvert(v, DestIsIntegral(), PyramidIsIntegral());
    }

protected:

    typedef typename NumericTraits<DestPixelType>::isIntegral DestIsIntegral;
    typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

    // test time with floating-point dithering: 100.01 sec
    // test time with integer dithering: 94.89 sec
    // Convert an integral pyramid pixel to an integral image pixel.
    inline DestPixelType doConvert(const PyramidPixelType &v, VigraTrueType, VigraTrueType) const {
        // Integer Dithering
        PyramidPixelType half = 1 << (PyramidFractionBits-1);
        PyramidPixelType quarter = 1 << (PyramidFractionBits-2);
        PyramidPixelType threeQuarter = 3 << (PyramidFractionBits-2);

        PyramidPixelType vFraction = v & ((1 << PyramidFractionBits) - 1);

        if ((vFraction >= quarter) && (vFraction < threeQuarter)) {
            PyramidPixelType random = (PyramidPixelType(::Twister()) & (half - 1)) + quarter;
            if (random <= vFraction) {
                return DestPixelType(NumericTraits<DestPixelType>::fromPromote((v >> PyramidFractionBits) + 1));
            } else {
                return DestPixelType(NumericTraits<DestPixelType>::fromPromote(v >> PyramidFractionBits));
            }
        }
        else if (vFraction >= quarter) {
            return DestPixelType(NumericTraits<DestPixelType>::fromPromote((v >> PyramidFractionBits) + 1));
        }
        else {
            return DestPixelType(NumericTraits<DestPixelType>::fromPromote(v >> PyramidFractionBits));
        }

    }

    // Convert a real pyramid pixel to an integral image pixel.
    inline DestPixelType doConvert(const PyramidPixelType &v, VigraTrueType, VigraFalseType) const {
        double d = dither(v);
        return NumericTraits<DestPixelType>::fromRealPromote(d);
    }

    // Convert an integral pyramid pixel to a real image pixel.
    inline DestPixelType doConvert(const PyramidPixelType &v, VigraFalseType, VigraTrueType) const {
        return convertFixedPointToDouble(v);
    }

    // Convert a real pyramid pixel to a real image pixel.
    inline DestPixelType doConvert(const PyramidPixelType &v, VigraFalseType, VigraFalseType) const {
        return v;
    }

    // Dithering is used to fool the eye into seeing gradients that are finer
    // than the precision of the pixel type.
    // This prevents the occurence of cleanly-bordered regions in the output where
    // the pixel values suddenly change from N to N+1.
    // Such regions are especially objectionable in the green channel of 8-bit images.
    inline double dither(const double &v) const {
        double vFraction = v - floor(v);
        // Only dither values within a certain range of the rounding cutoff point.
        if (vFraction > 0.25 && vFraction <= 0.75) {
            // Generate a random number between 0 and 0.5.
            double random = 0.5 * (double)::Twister() / UINT_MAX;
            if ((vFraction - 0.25) >= random) {
                return ceil(v);
            } else {
                return floor(v);
            }
        } else {
            return v;
        }
    }

    inline double convertFixedPointToDouble(const PyramidPixelType &v) const {
        return NumericTraits<PyramidPixelType>::toRealPromote(v) / (double)(1 << PyramidFractionBits);
    };

};

/** Wrapper for vector pixel types. */
template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToPyramidFunctor {

    typedef typename SrcVectorType::value_type SrcComponentType;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    //typedef typename EnblendNumericTraits<SrcVectorType>::PyramidPixelType PyramidVectorType;
    //typedef typename EnblendNumericTraits<SrcVectorType>::PyramidPixelComponentType PyramidPixelComponentType;
    //typedef typename EnblendNumericTraits<SrcVectorType>::ImagePixelComponentType SrcComponentType;
    //enum {PyramidIntegerBits = EnblendNumericTraits<SrcVectorType>::PyramidIntegerBits};
    //enum {PyramidFractionBits = EnblendNumericTraits<SrcVectorType>::PyramidFractionBits};
    typedef ConvertScalarToPyramidFunctor<SrcComponentType, PyramidComponentType,
            PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertVectorToPyramidFunctor() : cf() {}

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        return PyramidVectorType(cf(v.red()), cf(v.green()), cf(v.blue()));
    }

protected:
    ConvertFunctorType cf;
};

/** Wrapper for vector pixel types. */
template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertPyramidToVectorFunctor {

    typedef typename DestVectorType::value_type DestComponentType;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    //typedef typename EnblendNumericTraits<DestVectorType>::PyramidPixelType PyramidVectorType;
    //typedef typename EnblendNumericTraits<DestVectorType>::PyramidPixelComponentType PyramidPixelComponentType;
    //typedef typename EnblendNumericTraits<DestVectorType>::ImagePixelComponentType DestComponentType;
    //enum {PyramidIntegerBits = EnblendNumericTraits<DestVectorType>::PyramidIntegerBits};
    //enum {PyramidFractionBits = EnblendNumericTraits<DestVectorType>::PyramidFractionBits};
    typedef ConvertPyramidToScalarFunctor<DestComponentType, PyramidComponentType,
            PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertPyramidToVectorFunctor() : cf() {}

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        return DestVectorType(cf(v.red()), cf(v.green()), cf(v.blue()));
    }

protected:
    ConvertFunctorType cf;
};

/** Copy a scalar image into a scalar pyramid image. */
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
void copyToPyramidImage(
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename PyramidImageType::traverser dest_upperleft,
        typename PyramidImageType::Accessor da,
        VigraTrueType) {

    typedef typename SrcImageType::value_type SrcPixelType;
    typedef typename PyramidImageType::value_type PyramidPixelType;

    //typedef typename EnblendNumericTraits<SrcPixelType>::PyramidPixelType PyramidPixelType;
    //enum {PyramidIntegerBits = EnblendNumericTraits<SrcPixelType>::PyramidIntegerBits};
    //enum {PyramidFractionBits = EnblendNumericTraits<SrcPixelType>::PyramidFractionBits};

    transformImage(src_upperleft, src_lowerright, sa,
            dest_upperleft, da,
            ConvertScalarToPyramidFunctor<SrcPixelType, PyramidPixelType, PyramidIntegerBits, PyramidFractionBits>());
};

/** Copy a vector image into a vector pyramid image.
 *  Uses an optional color space conversion.
 */
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
void copyToPyramidImage(
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename PyramidImageType::traverser dest_upperleft,
        typename PyramidImageType::Accessor da,
        VigraFalseType) {

    typedef typename SrcImageType::value_type SrcVectorType;
    typedef typename PyramidImageType::value_type PyramidVectorType;

    transformImage(src_upperleft, src_lowerright, sa,
            dest_upperleft, da,
            ConvertVectorToPyramidFunctor<SrcVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());

};

// Compile-time switch based on scalar or vector image type.
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyToPyramidImage(
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename PyramidImageType::traverser dest_upperleft,
        typename PyramidImageType::Accessor da) {

    typedef typename NumericTraits<typename SrcImageType::value_type>::isScalar src_is_scalar;

    copyToPyramidImage<SrcImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits>(
            src_upperleft,
            src_lowerright,
            sa,
            dest_upperleft,
            da,
            src_is_scalar());
};

// Version using argument object factories.
template <typename SrcImageType, typename PyramidImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyToPyramidImage(
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        pair<typename PyramidImageType::traverser, typename PyramidImageType::Accessor> dest) {
    copyToPyramidImage<SrcImageType, PyramidImageType, PyramidIntegerBits, PyramidFractionBits>(
            src.first,
            src.second,
            src.third,
            dest.first,
            dest.second);
};

/** Copy a scalar pyramid image into a scalar image. */
template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        VigraTrueType) {

    typedef typename DestImageType::value_type DestPixelType;
    typedef typename PyramidImageType::value_type PyramidPixelType;
    //typedef typename EnblendNumericTraits<DestPixelType>::PyramidPixelType PyramidPixelType;
    //enum {PyramidIntegerBits = EnblendNumericTraits<DestPixelType>::PyramidIntegerBits};
    //enum {PyramidFractionBits = EnblendNumericTraits<DestPixelType>::PyramidFractionBits};

    transformImageIf(src_upperleft, src_lowerright, sa,
            mask_upperleft, ma,
            dest_upperleft, da,
            ConvertPyramidToScalarFunctor<DestPixelType, PyramidPixelType, PyramidIntegerBits, PyramidFractionBits>());

};

/** Copy a vector pyramid image into a vector image.
 *  Uses an optional color space conversion.
 */
template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        VigraFalseType) {

    typedef typename DestImageType::value_type DestVectorType;
    typedef typename PyramidImageType::value_type PyramidVectorType;

    transformImageIf(src_upperleft, src_lowerright, sa,
            mask_upperleft, ma,
            dest_upperleft, da,
            ConvertPyramidToVectorFunctor<DestVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());

};

// Compile-time switch based on scalar or vector image type.
template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da) {

    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar src_is_scalar;

    copyFromPyramidImageIf<PyramidImageType, MaskImageType, DestImageType, PyramidIntegerBits, PyramidFractionBits>(
            src_upperleft,
            src_lowerright,
            sa,
            mask_upperleft,
            ma,
            dest_upperleft,
            da,
            src_is_scalar());
};

// Version using argument object factories.
template <typename PyramidImageType, typename MaskImageType, typename DestImageType, int PyramidIntegerBits, int PyramidFractionBits>
inline void copyFromPyramidImageIf(
        triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
        pair<typename MaskImageType::const_traverser, typename MaskImageType::ConstAccessor> mask,
        pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest) {

    copyFromPyramidImageIf<PyramidImageType, MaskImageType, DestImageType, PyramidIntegerBits, PyramidFractionBits>(
            src.first,
            src.second,
            src.third,
            mask.first,
            mask.second,
            dest.first,
            dest.second);
};

} // namespace enblend

#endif /* __FIXMATH_H__ */
