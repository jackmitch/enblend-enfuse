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
        return NumericTraits<PyramidPixelType>::fromRealPromote(v * (double)(1 << FractionBits));
    };

    inline PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType &v) const {
        return (PyramidPixelType)v << FractionBits;
    };

};

template <typename DestPixelType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(DestPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertPyramidToScalarFunctor {
public:
    ConvertPyramidToScalarFunctor() { }

    inline DestPixelType operator()(const PyramidPixelType &v) const {
        return doConvert(v, DestIsIntegral(), PyramidIsIntegral());
    }

protected:
    typedef typename NumericTraits<DestPixelType>::isIntegral DestIsIntegral;
    typedef typename NumericTraits<PyramidPixelType>::isIntegral PyramidIsIntegral;

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

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(SrcPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertVectorToPyramidFunctor {
public:
    ConvertVectorToPyramidFunctor() : cf() {}

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        PyramidVectorType pv;
        SrcVectorIterator svi = v.begin();
        PyramidVectorIterator pvi = pv.begin();
        for (; svi != v.end(); ++svi, ++pvi) {
            *pvi = cf(*svi);
        }
        return pv;
    }

protected:
    typedef ConvertScalarToPyramidFunctor<SrcPixelType, PyramidPixelType, IntegerBits, FractionBits> ConvertFunctorType;
    typedef typename PyramidVectorType::iterator PyramidVectorIterator;
    typedef typename SrcVectorType::const_iterator SrcVectorIterator;

    ConvertFunctorType cf;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType, int IntegerBits=1+8*sizeof(DestPixelType), int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
class ConvertPyramidToVectorFunctor {
public:
    ConvertPyramidToVectorFunctor() : cf() {}

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        DestVectorType dv;
        PyramidVectorIterator pvi = v.begin();
        DestVectorIterator dvi = dv.begin();
        for (; pvi != v.end(); ++pvi, ++dvi) {
            *dvi = cf(*pvi);
        }
        return dv;
    }

protected:
    typedef ConvertPyramidToScalarFunctor<DestPixelType, PyramidPixelType, IntegerBits, FractionBits> ConvertFunctorType;
    typedef typename PyramidVectorType::const_iterator PyramidVectorIterator;
    typedef typename DestVectorType::iterator DestVectorIterator;

    ConvertFunctorType cf;
};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToLabPyramidFunctor {
public:
    ConvertToLabPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {
        if (Verbose > VERBOSE_COLOR_CONVERSION_MESSAGES) {
            cout << "R'G'B' to L*a*b* color space conversion..." << endl;
        }
    }

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        return convertFunctor(colorFunctor(v));
    }

protected:
    typedef RGBPrime2LabFunctor<SrcPixelType> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // L*a*b* components:
    //      0      <= L <= 100
    //    -86.1813 <= a <=  98.2352
    //   -107.862  <= b <=  94.4758
    // Needs 8 integer bits
    // Needs 9 integer bits for pyramid math
    typedef ConvertVectorToPyramidFunctor<ColorFunctorResultType, ColorFunctorResultComponent, PyramidVectorType, PyramidPixelType, 9> ConvertFunctorType;

    ColorFunctorType colorFunctor;
    ConvertFunctorType convertFunctor;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertFromLabPyramidFunctor {
public:
    ConvertFromLabPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {
        if (Verbose > VERBOSE_COLOR_CONVERSION_MESSAGES) {
            cout << "L*a*b* to R'G'B' color space conversion..." << endl;
        }
    }

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        return destFunctor(colorFunctor(doubleFunctor(v)));
    }

protected:
    typedef Lab2RGBPrimeFunctor<double> ColorFunctorType;
    typedef typename ColorFunctorType::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // L*a*b* fixed-point pyramid uses 9 integer bits.
    typedef ConvertPyramidToVectorFunctor<ColorFunctorArgumentType, double, PyramidVectorType, PyramidPixelType, 9> DoubleFunctorType;
    typedef ConvertPyramidToVectorFunctor<DestVectorType, DestPixelType, ColorFunctorResultType, ColorFunctorResultComponent, 9> DestFunctorType;

    ColorFunctorType colorFunctor;
    DoubleFunctorType doubleFunctor;
    DestFunctorType destFunctor;
};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToYPrimePbPrPyramidFunctor {
public:
    ConvertToYPrimePbPrPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {
        if (Verbose > VERBOSE_COLOR_CONVERSION_MESSAGES) {
            cout << "R'G'B' to Y'PbPr color space conversion..." << endl;
        }
    }

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        return convertFunctor(colorFunctor(v));
    }

protected:
    typedef RGBPrime2YPrimePbPrFunctor<SrcPixelType> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // Y'PbPr components:
    //      0 <= Y' <= 1
    //   -0.5 <= Pb <= 0.5
    //   -0.5 <= Pr <= 0.5
    // Needs 2 integer bits for pyramid math
    typedef ConvertVectorToPyramidFunctor<ColorFunctorResultType, ColorFunctorResultComponent, PyramidVectorType, PyramidPixelType, 2> ConvertFunctorType;

    ColorFunctorType colorFunctor;
    ConvertFunctorType convertFunctor;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertFromYPrimePbPrPyramidFunctor {
public:
    ConvertFromYPrimePbPrPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {
        if (Verbose > VERBOSE_COLOR_CONVERSION_MESSAGES) {
            cout << "Y'PbPr to R'G'B' color space conversion..." << endl;
        }
    }

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        return destFunctor(colorFunctor(doubleFunctor(v)));
    }

protected:
    typedef YPrimePbPr2RGBPrimeFunctor<double> ColorFunctorType;
    typedef typename ColorFunctorType::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // Y'PbPr fixed-point pyramid uses 2 integer bits.
    typedef ConvertPyramidToVectorFunctor<ColorFunctorArgumentType, double, PyramidVectorType, PyramidPixelType, 2> DoubleFunctorType;
    typedef ConvertPyramidToVectorFunctor<DestVectorType, DestPixelType, ColorFunctorResultType, ColorFunctorResultComponent, 2> DestFunctorType;

    ColorFunctorType colorFunctor;
    DoubleFunctorType doubleFunctor;
    DestFunctorType destFunctor;
};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToRGBPyramidFunctor {
public:
    ConvertToRGBPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {}

    inline PyramidVectorType operator()(const SrcVectorType &v) const {
        return convertFunctor(colorFunctor(v));
    }

protected:
    typedef RGBPrime2RGBFunctor<SrcPixelType, double> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // RGB components are the same range as R'G'B' components
    // Use default number of integer bits for pyramid math
    typedef ConvertVectorToPyramidFunctor<ColorFunctorResultType, ColorFunctorResultComponent, PyramidVectorType, PyramidPixelType, 1+8*sizeof(SrcPixelType)> ConvertFunctorType;

    ColorFunctorType colorFunctor;
    ConvertFunctorType convertFunctor;
};

template <typename DestVectorType, typename DestPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertFromRGBPyramidFunctor {
public:
    ConvertFromRGBPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {}

    inline DestVectorType operator()(const PyramidVectorType &v) const {
        return destFunctor(colorFunctor(doubleFunctor(v)));
    }

protected:
    typedef RGB2RGBPrimeFunctor<double, double> ColorFunctorType;
    typedef typename ColorFunctorType::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // RGB components are the same range as R'G'B' components
    // Use default number of integer bits for pyramid math
    typedef ConvertPyramidToVectorFunctor<ColorFunctorArgumentType, double, PyramidVectorType, PyramidPixelType, 1+8*sizeof(DestPixelType)> DoubleFunctorType;
    typedef ConvertPyramidToVectorFunctor<DestVectorType, DestPixelType, ColorFunctorResultType, ColorFunctorResultComponent, 1+8*sizeof(DestPixelType)> DestFunctorType;

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
                ConvertVectorToPyramidFunctor<SrcVectorType, SrcPixelType, PyramidVectorType, PyramidPixelType>());
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
                ConvertPyramidToVectorFunctor<DestVectorType, DestPixelType, PyramidVectorType, PyramidPixelType>());
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
