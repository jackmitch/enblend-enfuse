/*
 * Copyright (C) 2004 Andrew Mihal
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
using vigra::VigraFalseType;
using vigra::VigraTrueType;
using vigra::YPrimePbPr2RGBPrimeFunctor;

namespace enblend {

template <typename SrcPixelType, typename PyramidPixelType, int SrcIntegerBits=1+8*sizeof(SrcPixelType), int FractionBits=8*sizeof(PyramidPixelType)-SrcIntegerBits>
class ConvertScalarToPyramidFunctor {
public:
    typedef typename NumericTraits<SrcPixelType>::isIntegral SrcIsIntegral;
    typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

    ConvertScalarToPyramidFunctor() { }

    inline PyramidPixelType operator()(const SrcPixelType &v) const {
        return doConvert(v, SrcIsIntegral(), PyramidIsIntegral());
    }

protected:

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
        return NumericTraits<PyramidPixelType>::fromRealPromote(v * (double)(1 << FractionBits));
    };

    inline PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType &v) const {
        return (PyramidPixelType)v << FractionBits;
    };

};

template <typename DestPixelType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(DestPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertPyramidToScalarFunctor {
public:
    typedef typename NumericTraits<DestPixelType>::isIntegral DestIsIntegral;
    typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

    ConvertPyramidToScalarFunctor() { }

    inline DestPixelType operator()(const PyramidPixelType &v) const {
        return doConvert(v, DestIsIntegral(), PyramidIsIntegral());
    }

protected:

    // Convert an integral pyramid pixel to an integral image pixel.
    inline DestPixelType doConvert(const PyramidPixelType &v, VigraTrueType, VigraTrueType) const {
        double d = convertFixedPointToDouble(v);
        d = dither(d);
        return NumericTraits<DestPixelType>::fromRealPromote(d);
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

    inline double dither(const double &v) const {
        double vFraction = v - floor(v);
        // Only dither values within a certain range of the rounding cutoff point.
        if (vFraction > 0.25 && vFraction <= 0.75) {
            // Generate a random number between 0 and 0.5.
            double random = 0.5 * (double)Twister() / UINT_MAX;
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

};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToLabPyramidFunctor {
public:
    typedef RGBPrime2LabFunctor<SrcPixelType> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // L*a*b* components:
    //      0      <= L <= 100
    //    -86.1813 <= a <=  98.2352
    //   -107.862  <= b <=  94.4758
    // Needs 8 integer bits
    // Needs 9 integer bits for pyramid math
    typedef ConvertScalarToPyramidFunctor<ColorFunctorResultComponent, PyramidPixelType, 9> ConvertFunctorType;

    ConvertToLabPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {
        if (Verbose > 0) {
            cout << "R'G'B' to L*a*b* color space conversion..." << endl;
        }
    }

    PyramidVectorType operator()(const SrcVectorType &v) const {
        ColorFunctorResultType labVector = colorFunctor(v);
        PyramidPixelType l = convertFunctor(labVector[0]);
        PyramidPixelType a = convertFunctor(labVector[1]);
        PyramidPixelType b = convertFunctor(labVector[2]);
        return PyramidVectorType(l, a, b);
    }

protected:
    ColorFunctorType colorFunctor;
    ConvertFunctorType convertFunctor;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertFromLabPyramidFunctor {
public:
    typedef Lab2RGBPrimeFunctor<double> ColorFunctorType;
    typedef typename ColorFunctorType::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // L*a*b* fixed-point pyramid uses 9 integer bits.
    typedef ConvertPyramidToScalarFunctor<double, PyramidPixelType, 9> DoubleFunctorType;
    typedef ConvertPyramidToScalarFunctor<DestPixelType, double, 9> DestFunctorType;

    ConvertFromLabPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {
        if (Verbose > 0) {
            cout << "L*a*b* to R'G'B' color space conversion..." << endl;
        }
    }

    DestVectorType operator()(const PyramidVectorType &v) const {
        double l = doubleFunctor(v.red());
        double a = doubleFunctor(v.green());
        double b = doubleFunctor(v.blue());

        ColorFunctorArgumentType labVector(l, a, b);
        ColorFunctorResultType rgbpVector = colorFunctor(labVector);

        DestPixelType red = destFunctor(rgbpVector.red());
        DestPixelType green = destFunctor(rgbpVector.green());
        DestPixelType blue = destFunctor(rgbpVector.blue());

        return DestVectorType(red, green, blue);
    }

protected:
    ColorFunctorType colorFunctor;
    DoubleFunctorType doubleFunctor;
    DestFunctorType destFunctor;
};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToYPrimePbPrPyramidFunctor {
public:
    typedef RGBPrime2YPrimePbPrFunctor<SrcPixelType> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // Y'PbPr components:
    //      0 <= Y' <= 1
    //   -0.5 <= Pb <= 0.5
    //   -0.5 <= Pr <= 0.5
    // Needs 2 integer bits for pyramid math
    typedef ConvertScalarToPyramidFunctor<ColorFunctorResultComponent, PyramidPixelType, 2> ConvertFunctorType;

    ConvertToYPrimePbPrPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {
        if (Verbose > 0) {
            cout << "R'G'B' to Y'PbPr color space conversion..." << endl;
        }
    }

    PyramidVectorType operator()(const SrcVectorType &v) const {
        ColorFunctorResultType ypbprVector = colorFunctor(v);
        PyramidPixelType y  = convertFunctor(ypbprVector[0]);
        PyramidPixelType pb = convertFunctor(ypbprVector[1]);
        PyramidPixelType pr = convertFunctor(ypbprVector[2]);
        return PyramidVectorType(y, pb, pr);
    }

protected:
    ColorFunctorType colorFunctor;
    ConvertFunctorType convertFunctor;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertFromYPrimePbPrPyramidFunctor {
public:
    typedef YPrimePbPr2RGBPrimeFunctor<double> ColorFunctorType;
    typedef typename ColorFunctorType::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // Y'PbPr fixed-point pyramid uses 2 integer bits.
    typedef ConvertPyramidToScalarFunctor<double, PyramidPixelType, 2> DoubleFunctorType;
    typedef ConvertPyramidToScalarFunctor<DestPixelType, double, 2> DestFunctorType;

    ConvertFromYPrimePbPrPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {
        if (Verbose > 0) {
            cout << "Y'PbPr to R'G'B' color space conversion..." << endl;
        }
    }

    DestVectorType operator()(const PyramidVectorType &v) const {
        double y  = doubleFunctor(v.red());
        double pb = doubleFunctor(v.green());
        double pr = doubleFunctor(v.blue());

        ColorFunctorArgumentType ypbprVector(y, pb, pr);
        ColorFunctorResultType rgbpVector = colorFunctor(ypbprVector);

        DestPixelType red = destFunctor(rgbpVector.red());
        DestPixelType green = destFunctor(rgbpVector.green());
        DestPixelType blue = destFunctor(rgbpVector.blue());

        return DestVectorType(red, green, blue);
    }

protected:
    ColorFunctorType colorFunctor;
    DoubleFunctorType doubleFunctor;
    DestFunctorType destFunctor;
};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToRGBPyramidFunctor {
public:
    typedef RGBPrime2RGBFunctor<SrcPixelType, double> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // RGB components are the same range as R'G'B' components
    // Use default number of integer bits for pyramid math
    typedef ConvertScalarToPyramidFunctor<ColorFunctorResultComponent, PyramidPixelType, 1+8*sizeof(SrcPixelType)> ConvertFunctorType;

    ConvertToRGBPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {}

    PyramidVectorType operator()(const SrcVectorType &v) const {
        ColorFunctorResultType rgbVector = colorFunctor(v);
        PyramidPixelType r = convertFunctor(rgbVector[0]);
        PyramidPixelType g = convertFunctor(rgbVector[1]);
        PyramidPixelType b = convertFunctor(rgbVector[2]);
        return PyramidVectorType(r, g, b);
    }

protected:
    ColorFunctorType colorFunctor;
    ConvertFunctorType convertFunctor;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertFromRGBPyramidFunctor {
public:
    typedef RGB2RGBPrimeFunctor<double, double> ColorFunctorType;
    typedef typename ColorFunctorType::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // RGB components are the same range as R'G'B' components
    // Use default number of integer bits for pyramid math
    typedef ConvertPyramidToScalarFunctor<double, PyramidPixelType, 1+8*sizeof(DestPixelType)> DoubleFunctorType;
    typedef ConvertPyramidToScalarFunctor<DestPixelType, double, 1+8*sizeof(DestPixelType)> DestFunctorType;

    ConvertFromRGBPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {}

    DestVectorType operator()(const PyramidVectorType &v) const {
        double r  = doubleFunctor(v.red());
        double g = doubleFunctor(v.green());
        double b = doubleFunctor(v.blue());

        ColorFunctorArgumentType rgbVector(r, g, b);
        ColorFunctorResultType rgbpVector = colorFunctor(rgbVector);

        DestPixelType red = destFunctor(rgbpVector.red());
        DestPixelType green = destFunctor(rgbpVector.green());
        DestPixelType blue = destFunctor(rgbpVector.blue());

        return DestVectorType(red, green, blue);
    }

protected:
    ColorFunctorType colorFunctor;
    DoubleFunctorType doubleFunctor;
    DestFunctorType destFunctor;
};

// copy scalar image to scalar fixed-point pyramid image.
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

// Copy vector image to vector fixed-point pyramid image.
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

    if (UseLabColor) {
        transformImage(src_upperleft, src_lowerright, sa,
                dest_upperleft, da,
                ConvertToLabPyramidFunctor<SrcVectorType, SrcPixelType, PyramidVectorType, PyramidPixelType>());
    } else {
        transformImage(src_upperleft, src_lowerright, sa,
                dest_upperleft, da,
                ConvertToRGBPyramidFunctor<SrcVectorType, SrcPixelType, PyramidVectorType, PyramidPixelType>());
    }

};

// Switch based on vector or scalar image types.
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

// Using argument object factory
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

// copy scalar fixed-point pyramid image to scalar image.
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

// copy vector fixed-point pyramid image to vector image.
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

    if (UseLabColor) {
        transformImageIf(src_upperleft, src_lowerright, sa,
                mask_upperleft, ma,
                dest_upperleft, da,
                ConvertFromLabPyramidFunctor<DestVectorType, DestPixelType, PyramidVectorType, PyramidPixelType>());
    } else {
        transformImageIf(src_upperleft, src_lowerright, sa,
                mask_upperleft, ma,
                dest_upperleft, da,
                ConvertFromRGBPyramidFunctor<DestVectorType, DestPixelType, PyramidVectorType, PyramidPixelType>());
    }

};

// Switch based on vector or scalar image types.
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

// Using argument object factory
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
