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
#include "vigra/mathutil.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/utilities.hxx"

using std::pair;

using vigra::NumericTraits;
using vigra::triple;

namespace enblend {

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

/** Fixed point converter that uses ICC profile transformation */
template <typename SrcVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertVectorToJCHPyramidFunctor {

    typedef typename SrcVectorType::value_type SrcComponentType;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertScalarToPyramidFunctor<double, PyramidComponentType,
            PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertVectorToJCHPyramidFunctor() : cf() {
        scale = 1.0 / NumericTraits<SrcComponentType>::toRealPromote(
                            NumericTraits<SrcComponentType>::max());
    }

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        // rgb values must be in range [0,1]
        double rgb[3];
        rgb[0] = scale * NumericTraits<SrcComponentType>::toRealPromote(v.red());
        rgb[1] = scale * NumericTraits<SrcComponentType>::toRealPromote(v.green());
        rgb[2] = scale * NumericTraits<SrcComponentType>::toRealPromote(v.blue());

        double xyz[3];
        cmsDoTransform(InputToXYZTransform, rgb, xyz, 1);
        // xyz values are in range [0,1]

        // relative xyz values must be in range [0,100]
        cmsCIEXYZ cmsxyz;
        cmsxyz.X = xyz[0] * 100.0;
        cmsxyz.Y = xyz[1] * 100.0;
        cmsxyz.Z = xyz[2] * 100.0;

        cmsJCh jch;
        cmsCIECAM02Forward(CIECAMTransform, &cmsxyz, &jch);
        // J in range [0,100], C in range [0,120], h in range [0,360]

        // convert cylindrical to cartesian
        double theta = jch.h * M_PI / 180.0;
        jch.h = jch.C * cos(theta);
        jch.C = jch.C * sin(theta);

        // Scale to maximize usage of fixed-point type
        jch.J *= exp2(PyramidIntegerBits - 1 - 7);
        jch.C *= exp2(PyramidIntegerBits - 1 - 7);
        jch.h *= exp2(PyramidIntegerBits - 1 - 7);

        return PyramidVectorType(cf(jch.J), cf(jch.C), cf(jch.h));
    }

protected:
    ConvertFunctorType cf;
    double scale;
};

/** Fixed point converter that uses ICC profile transformation */
template <typename DestVectorType, typename PyramidVectorType, int PyramidIntegerBits, int PyramidFractionBits>
class ConvertJCHPyramidToVectorFunctor {

    typedef typename DestVectorType::value_type DestComponentType;
    typedef typename PyramidVectorType::value_type PyramidComponentType;
    typedef ConvertPyramidToScalarFunctor<double, PyramidComponentType,
            PyramidIntegerBits, PyramidFractionBits> ConvertFunctorType;

public:
    ConvertJCHPyramidToVectorFunctor() : cf() {
        scale = NumericTraits<DestComponentType>::toRealPromote(
                        NumericTraits<DestComponentType>::max());
    }

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        cmsJCh jch;
        jch.J = cf(v.red());
        jch.C = cf(v.green());
        jch.h = cf(v.blue());

        // Scale back to range J[0,100], C[0,120], h[0,120]
        jch.J /= exp2(PyramidIntegerBits - 1 - 7);
        jch.C /= exp2(PyramidIntegerBits - 1 - 7);
        jch.h /= exp2(PyramidIntegerBits - 1 - 7);

        // convert cartesian to cylindrical
        double r = sqrt(jch.C * jch.C + jch.h * jch.h);
        jch.h = (180.0 / M_PI) * atan2(jch.C, jch.h);
        if (jch.h < 0.0) jch.h += 360.0;
        jch.C = r;

        cmsCIEXYZ cmsxyz;
        cmsCIECAM02Reverse(CIECAMTransform, &jch, &cmsxyz);
        // xyz values in range [0,100]

        // scale xyz values to range [0,1]
        double xyz[3];
        xyz[0] = cmsxyz.X / 100.0;
        xyz[1] = cmsxyz.Y / 100.0;
        xyz[2] = cmsxyz.Z / 100.0;

        double rgb[3];
        cmsDoTransform(XYZToInputTransform, xyz, rgb, 1);
        // rgb values in range [0,1]

        return DestVectorType(NumericTraits<DestComponentType>::fromRealPromote(rgb[0] * scale),
                              NumericTraits<DestComponentType>::fromRealPromote(rgb[1] * scale),
                              NumericTraits<DestComponentType>::fromRealPromote(rgb[2] * scale));
    }

protected:
    ConvertFunctorType cf;
    double scale;
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

    if (UseCIECAM) {
        if (Verbose > VERBOSE_COLOR_CONVERSION_MESSAGES) {
            cout << "CIECAM02 color conversion..." << endl;
        }
        transformImage(src_upperleft, src_lowerright, sa,
                dest_upperleft, da,
                ConvertVectorToJCHPyramidFunctor<SrcVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());
    } else {
        transformImage(src_upperleft, src_lowerright, sa,
                dest_upperleft, da,
                ConvertVectorToPyramidFunctor<SrcVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());
    }

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

    if (UseCIECAM) {
        if (Verbose > VERBOSE_COLOR_CONVERSION_MESSAGES) {
            cout << "CIECAM02 color conversion..." << endl;
        }
        transformImageIf(src_upperleft, src_lowerright, sa,
                mask_upperleft, ma,
                dest_upperleft, da,
                ConvertJCHPyramidToVectorFunctor<DestVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());
    } else {
        transformImageIf(src_upperleft, src_lowerright, sa,
                mask_upperleft, ma,
                dest_upperleft, da,
                ConvertPyramidToVectorFunctor<DestVectorType, PyramidVectorType, PyramidIntegerBits, PyramidFractionBits>());
    }

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
