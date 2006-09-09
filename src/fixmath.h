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

#include "vigra/colorconversions.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/utilities.hxx"

using std::cout;
using std::pair;

using vigra::Lab2RGBPrimeFunctor;
using vigra::NumericTraits;
using vigra::RGB2RGBPrimeFunctor;
using vigra::RGBPrime2LabFunctor;
using vigra::RGBPrime2RGBFunctor;
using vigra::RGBPrime2YPrimePbPrFunctor;
using vigra::triple;
using vigra::UInt32;
using vigra::VigraFalseType;
using vigra::VigraTrueType;
using vigra::YPrimePbPr2RGBPrimeFunctor;

namespace enblend {

/** A functor for converting scalar pixel values to the number representation used
 *  for pyramids. These are either fixed-point integers or floating-point numberss.
 *  For fixed-point integers, use the SrcIntegerBits to indicate how many bits are
 *  necessary to represent the integer part of source pixel data.
 *  One additional bit of range above the maximum range of the source pixels is
 *  necessary for Laplacian pyramid calculations.
 */
template <typename SrcPixelType, typename PyramidPixelType, int SrcIntegerBits=1+8*sizeof(SrcPixelType), int FractionBits=8*sizeof(PyramidPixelType)-SrcIntegerBits>
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
        return NumericTraits<PyramidPixelType>::fromRealPromote(v * (double)(1 << FractionBits));
    };

    inline PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType &v) const {
        // Shift v left to move the decimal point and set the fraction bits to zero.
        return (PyramidPixelType)v << FractionBits;
    };

};

/** A functor for converting numbers stored in the pyramid number representation back
 *  into normal pixel values.
 *  The IntegerBits parameter is supposed to match the corresponding parameter used in
 *  the ConvertScalarToPyramidFunctor that originally coverted the pixels to the
 *  pyramid representation.
 */
template <typename DestPixelType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(DestPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertPyramidToScalarFunctor {
public:
    ConvertPyramidToScalarFunctor() /*: overflows(0), underflows(0)*/ { }

    //~ConvertPyramidToScalarFunctor() {
    //    cout << "overflows=" << overflows << endl;
    //    cout << "underflows=" << underflows << endl;
    //}

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
        PyramidPixelType half = 1 << (FractionBits-1);
        PyramidPixelType quarter = 1 << (FractionBits-2);
        PyramidPixelType threeQuarter = 3 << (FractionBits-2);

        PyramidPixelType vFraction = v & ((1 << FractionBits) - 1);

        if ((vFraction >= quarter) && (vFraction < threeQuarter)) {
            PyramidPixelType random = (PyramidPixelType(::Twister()) & (half - 1)) + quarter;
            if (random <= vFraction) {
                return DestPixelType(v >> FractionBits) + 1;
            } else {
                return DestPixelType(v >> FractionBits);
            }
        }
        else if (vFraction >= quarter) {
            return DestPixelType(v >> FractionBits) + 1;
        }
        else {
            return DestPixelType(v >> FractionBits);
        }

        // Floating-point dithering
        //double d = convertFixedPointToDouble(v);
        //d = dither(d);
        ////if (d > NumericTraits<DestPixelType>::max()) overflows++;
        ////if (d < NumericTraits<DestPixelType>::min()) underflows++;
        //return NumericTraits<DestPixelType>::fromRealPromote(d);
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
        return NumericTraits<PyramidPixelType>::toRealPromote(v) / (double)(1 << FractionBits);
    };

    //mutable int overflows;
    //mutable int underflows;

};

/** Wrapper for vector pixel types. */
template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(SrcPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertVectorToPyramidFunctor {
public:
    ConvertVectorToPyramidFunctor() : cf() {}

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        //PyramidVectorType pv;
        //SrcVectorIterator svi = v.begin();
        //PyramidVectorIterator pvi = pv.begin();
        //for (; svi != v.end(); ++svi, ++pvi) {
        //    *pvi = cf(*svi);
        //}
        //return pv;
        return PyramidVectorType(cf(v.red()), cf(v.green()), cf(v.blue()));
    }

protected:
    typedef ConvertScalarToPyramidFunctor<SrcPixelType, PyramidPixelType, IntegerBits, FractionBits> ConvertFunctorType;
    typedef typename PyramidVectorType::iterator PyramidVectorIterator;
    typedef typename SrcVectorType::const_iterator SrcVectorIterator;

    ConvertFunctorType cf;
};

/** Wrapper for vector pixel types. */
template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(DestPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertPyramidToVectorFunctor {
public:
    ConvertPyramidToVectorFunctor() : cf() {}

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        //DestVectorType dv;
        //PyramidVectorIterator pvi = v.begin();
        //DestVectorIterator dvi = dv.begin();
        //for (; pvi != v.end(); ++pvi, ++dvi) {
        //    *dvi = cf(*pvi);
        //}
        //return dv;
        return DestVectorType(cf(v.red()), cf(v.green()), cf(v.blue()));
    }

protected:
    typedef ConvertPyramidToScalarFunctor<DestPixelType, PyramidPixelType, IntegerBits, FractionBits> ConvertFunctorType;
    typedef typename PyramidVectorType::const_iterator PyramidVectorIterator;
    typedef typename DestVectorType::iterator DestVectorIterator;

    ConvertFunctorType cf;
};

/** Copy a scalar image into a scalar pyramid image. */
template <typename SrcImageType, typename PyramidImageType>
void copyToPyramidImage(
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename PyramidImageType::traverser dest_upperleft,
        typename PyramidImageType::Accessor da,
        VigraTrueType) {

    typedef typename SrcImageType::value_type SrcPixelType;
    typedef typename PyramidImageType::value_type PyramidPixelType;

    transformImage(src_upperleft, src_lowerright, sa,
            dest_upperleft, da,
            ConvertScalarToPyramidFunctor<SrcPixelType, PyramidPixelType>());
};

/** Copy a vector image into a vector pyramid image.
 *  Uses an optional color space conversion.
 */
template <typename SrcImageType, typename PyramidImageType>
void copyToPyramidImage(
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename PyramidImageType::traverser dest_upperleft,
        typename PyramidImageType::Accessor da,
        VigraFalseType) {

    typedef typename SrcImageType::value_type SrcVectorType;
    typedef typename SrcVectorType::value_type SrcPixelType;
    typedef typename PyramidImageType::value_type PyramidVectorType;
    typedef typename PyramidVectorType::value_type PyramidPixelType;

    //if (UseLabColor) {
    //    transformImage(src_upperleft, src_lowerright, sa,
    //            dest_upperleft, da,
    //            ConvertToLabPyramidFunctor<SrcVectorType, SrcPixelType, PyramidVectorType, PyramidPixelType>());
    //} else {
        transformImage(src_upperleft, src_lowerright, sa,
                dest_upperleft, da,
                ConvertVectorToPyramidFunctor<SrcVectorType, SrcPixelType, PyramidVectorType, PyramidPixelType>());
    //}

};

// Compile-time switch based on scalar or vector image type.
template <typename SrcImageType, typename PyramidImageType>
inline void copyToPyramidImage(
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename PyramidImageType::traverser dest_upperleft,
        typename PyramidImageType::Accessor da) {

    typedef typename NumericTraits<typename SrcImageType::value_type>::isScalar src_is_scalar;

    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src_upperleft,
            src_lowerright,
            sa,
            dest_upperleft,
            da,
            src_is_scalar());
};

// Version using argument object factories.
template <typename SrcImageType, typename PyramidImageType>
inline void copyToPyramidImage(
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        pair<typename PyramidImageType::traverser, typename PyramidImageType::Accessor> dest) {
    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src.first,
            src.second,
            src.third,
            dest.first,
            dest.second);
};

/** Copy a scalar pyramid image into a scalar image. */
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        VigraTrueType) {

    typedef typename PyramidImageType::value_type PyramidPixelType;
    typedef typename DestImageType::value_type DestPixelType;

    transformImageIf(src_upperleft, src_lowerright, sa,
            mask_upperleft, ma,
            dest_upperleft, da,
            ConvertPyramidToScalarFunctor<DestPixelType, PyramidPixelType>());

};

/** Copy a vector pyramid image into a vector image.
 *  Uses an optional color space conversion.
 */
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        VigraFalseType) {

    typedef typename PyramidImageType::value_type PyramidVectorType;
    typedef typename PyramidVectorType::value_type PyramidPixelType;
    typedef typename DestImageType::value_type DestVectorType;
    typedef typename DestVectorType::value_type DestPixelType;

    //if (UseLabColor) {
    //    transformImageIf(src_upperleft, src_lowerright, sa,
    //            mask_upperleft, ma,
    //            dest_upperleft, da,
    //            ConvertFromLabPyramidFunctor<DestVectorType, DestPixelType, PyramidVectorType, PyramidPixelType>());
    //} else {
        transformImageIf(src_upperleft, src_lowerright, sa,
                mask_upperleft, ma,
                dest_upperleft, da,
                ConvertPyramidToVectorFunctor<DestVectorType, DestPixelType, PyramidVectorType, PyramidPixelType>());
    //}

};

// Compile-time switch based on scalar or vector image type.
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da) {

    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar src_is_scalar;

    copyFromPyramidImageIf<DestImageType, PyramidImageType, MaskImageType>(
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
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
        pair<typename MaskImageType::const_traverser, typename MaskImageType::ConstAccessor> mask,
        pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest) {
    copyFromPyramidImageIf<DestImageType, PyramidImageType, MaskImageType>(
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
