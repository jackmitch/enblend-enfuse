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

#include <boost/random/mersenne_twister.hpp>
#include <boost/static_assert.hpp>
#include "vigra/colorconversions.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/utilities.hxx"

using std::pair;

using vigra::NumericTraits;
//using vigra::RGBPrime2YPrimeCbCrFunctor;
using vigra::RGBPrime2LabFunctor;
using vigra::Lab2RGBPrimeFunctor;
//using vigra::RGBPrime2RGBFunctor;
//using vigra::RGB2RGBPrimeFunctor;
using vigra::triple;
using vigra::VigraFalseType;
using vigra::VigraTrueType;
//using vigra::YPrimeCbCr2RGBPrimeFunctor;

namespace enblend {

template <typename SrcPixelType, typename PyramidPixelType, int SrcIntegerBits=1+8*sizeof(SrcPixelType)>
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

    template <int FractionBits=8*sizeof(PyramidPixelType)-SrcIntegerBits>
    inline PyramidPixelType convertDoubleToFixedPoint(const double &v) {
        return NumericTraits<PyramidPixelType>::fromRealPromote(v * (double)(1 << FractionBits));
    };

    template <int FractionBits=8*sizeof(PyramidPixelType)-SrcIntegerBits>
    inline PyramidPixelType convertIntegerToFixedPoint(const SrcPixelType &v) {
        return (PyramidPixelType)v << FractionBits;
    };

};

template <typename DestPixelType, typename PyramidPixelType, int IntegerBits>
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

    template <int FractionBits=8*sizeof(PyramidPixelType)-IntegerBits>
    inline double convertFixedPointToDouble(const PyramidPixelType &v) {
        return NumericTraits<PyramidPixelType>::toRealPromote(v) / (double)(1 << FractionBits);
    };

};

template <typename SrcVectorType, typename SrcPixelType, typename PyramidVectorType, typename PyramidPixelType>
class ConvertToLabPyramidFunctor {
public:
    typedef RBGPrime2LabFunctor<SrcPixelType> ColorFunctorType;
    typedef typename ColorFunctorType::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // L*a*b* components:
    //      0      <= L <= 100
    //    -86.1813 <= a <=  98.2352
    //   -107.862  <= b <=  94.4758
    // Needs 8 integer bits
    // Needs 9 integer bits for pyramid math
    typedef ConvertScalarToPyramidFunctor<ColorFunctorResultComponent, PyramidPixelType, 9> ConvertFunctorType;

    ConvertToLabPyramidFunctor() : colorFunctor(NumericTraits<SrcPixelType>::max()), convertFunctor() {}

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

    ConvertFromLabPyramidFunctor() : colorFunctor(NumericTraits<DestPixelType>::max()), doubleFunctor(), destFunctor() {}

    DestVectorType operator()(const PyramidVectorType &v) const {
        double l = doubleFunctor(v.red());
        double a = doubleFunctor(v.green());
        double b = doubleFunctor(v.blue());

        ColorFunctorArgumentType labVector(l, a, b);
        ColorFunctorResultType rgbpVector = colorFunctor(labVector);

        DestPixelType r = destFunctor(rgbpVector.red());
        DestPixelType g = destFunctor(rgbpVector.green());
        DestPixelType b = destFunctor(rgbpVector.blue());

        return DestVectorType(r, g, b);
    }

protected:
    ColorFunctorType colorFunctor;
    DoubleFunctorType doubleFunctor;
    DestFunctorType destFunctor;
};









// Convert type T1 to fixed-point value stored in type T2.
template <typename T1, typename T2>
inline T2 convertToPyramidMath(T1 v) {
    // Fixed-point math conversion not defined for this type pair.
    BOOST_STATIC_ASSERT(false);
};

template <typename T1>
inline double convertFromPyramidMath(T1 v) {
    // Fixed-point math conversion not defined for this type pair.
    BOOST_STATIC_ASSERT(false);
};

/* convertToPyramidMath functions expect v to be in the range -256 to 255.
 * They are encoded as a signed fixed-point type with 9 integer bits and
 * as many fraction bits as possible.
 */
template <>
inline short convertToPyramidMath<double, short>(double v) {
    return NumericTraits<short>::fromRealPromote(v * 128.0);
};

template <>
inline short convertToPyramidMath<float, short>(float v) {
    return NumericTraits<short>::fromRealPromote(v * 128.0);
};

template <>
inline double convertFromPyramidMath<short>(short v) {
    return NumericTraits<short>::toRealPromote(v) / 128.0;
};

template <>
inline int convertToPyramidMath<double, int>(double v) {
    return NumericTraits<int>::fromRealPromote(v * 8388608.0);
};

template <>
inline int convertToPyramidMath<float, int>(float v) {
    return NumericTraits<int>::fromRealPromote(v * 8388608.0);
};

template <>
inline double convertFromPyramidMath<int>(int v) {
    return NumericTraits<int>::toRealPromote(v) / 8388608.0;
};

template <>
inline double convertToPyramidMath<double, double>(double v) {
    return v;
};

template <>
inline double convertToPyramidMath<float, double>(float v) {
    return v;
};

template <>
inline double convertFromPyramidMath<double>(double v) {
    return v;
};

template <>
inline short convertToPyramidMath<unsigned char, short>(unsigned char v) {
    // cast to short does not sign-extend
    // Shift left 7 bits
    return (short)v << 7;
};

template <>
inline int convertToPyramidMath<unsigned char, int>(unsigned char v) {
    // Shift left 23 bits
    return (int)v << 23;
};

template <>
inline double convertToPyramidMath<unsigned char, double>(unsigned char v) {
    return (double)v / NumericTraits<unsigned char>::max();
};

template <>
inline int convertToPyramidMath<short, int>(short v) {
    // cast to int sign-extends
    // Shift left 15 bits
    return (int)v << 15;
};

template <>
inline int convertToPyramidMath<unsigned short, int>(unsigned short v) {
    // cast to int does not sign-extend
    // Shift left 15 bits
    return (int)v << 15;
};

template <>
inline double convertToPyramidMath<int, double>(int v) {
    return NumericTraits<int>::toRealPromote(v);
};

template <>
inline double convertToPyramidMath<unsigned int, double>(unsigned int v) {
    return NumericTraits<unsigned int>::toRealPromote(v);
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

    typedef typename SrcImageType::const_traverser SrcTraverser;
    typedef typename PyramidImageType::traverser DestTraverser;

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;

        for (; sx.x != send.x; ++sx.x, ++dx.x) {
            da.set(convertToPyramidMath<SrcPixelType, PyramidPixelType>(sa(sx)), dx);
        }
    }

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

    typedef typename SrcVectorType::const_iterator SrcVectorIterator;
    typedef typename PyramidVectorType::iterator PyramidVectorIterator;

    typedef typename SrcImageType::const_traverser SrcTraverser;
    typedef typename PyramidImageType::traverser DestTraverser;

    typedef RGBPrime2LabFunctor<SrcPixelType> ColorFunctor;
    typedef typename ColorFunctor::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    // Functor for color space conversion.
    // L*a*b* components need a 8-bit signed integer part and arbitrary fractional bits.
    ColorFunctor cf(NumericTraits<SrcPixelType>::max());

    if (Verbose > 0) {
        cout << "R'G'B' to L*a*b* color space conversion..." << endl;
    }

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;

        for (; sx.x != send.x; ++sx.x, ++dx.x) {
            ColorFunctorResultType labVector = cf(sa(sx));

            PyramidPixelType l = convertToPyramidMath<ColorFunctorResultComponent, PyramidPixelType>(labVector[0]);
            PyramidPixelType a = convertToPyramidMath<ColorFunctorResultComponent, PyramidPixelType>(labVector[1]);
            PyramidPixelType b = convertToPyramidMath<ColorFunctorResultComponent, PyramidPixelType>(labVector[2]);

            PyramidVectorType r(l, a, b);
            da.set(r, dx);

            //r = convertToPyramidMath<double, PyramidPixelType>(f(sa(sx)));
            //ColorFunctorResult converted = f(sa(sx));
            //unsigned char yp = NumericTraits<unsigned char>::fromRealPromote(converted[0]);
            //unsigned char cb = NumericTraits<unsigned char>::fromRealPromote(converted[1]);
            //unsigned char cr = NumericTraits<unsigned char>::fromRealPromote(converted[2]);
            //r.setRed(convertToPyramidMath<unsigned char, PyramidPixelType>(yp));
            //r.setGreen(convertToPyramidMath<unsigned char, PyramidPixelType>(cb));
            //r.setBlue(convertToPyramidMath<unsigned char, PyramidPixelType>(cr));

            //SrcVectorIterator svi = sa(sx).begin();
            //PyramidVectorIterator pvi = r.begin();
            //for (; pvi != r.end(); ++pvi, ++svi) {
            //    *pvi = convertToPyramidMath<SrcPixelType, PyramidPixelType>(*svi);
            //}

            //da.set(r, dx);
        }
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
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        VigraTrueType) {

    typedef typename PyramidImageType::value_type PyramidPixelType;
    typedef typename DestImageType::value_type DestPixelType;
    typedef typename NumericTraits<DestPixelType>::isIntegral DestPixelIsIntegral;

    typedef typename PyramidImageType::const_traverser SrcTraverser;
    typedef typename DestImageType::traverser DestTraverser;
    typedef typename MaskImageType::const_traverser MaskTraverser;

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;
    MaskTraverser my = mask_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y, ++my.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;
        MaskTraverser mx = my;

        for (; sx.x != send.x; ++sx.x, ++dx.x, ++mx.x) {
            if (ma(mx)) {
                double p = convertFromPyramidMath<PyramidPixelType>(sa(sx));
                p = dither(p, DestPixelIsIntegral());
                da.set(NumericTraits<DestPixelType>::fromRealPromote(p), dx);
            }
        }
    }

};

// copy vector fixed-point pyramid image to vector image.
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma,
        VigraFalseType) {

    typedef typename PyramidImageType::value_type PyramidVectorType;
    typedef typename PyramidVectorType::value_type PyramidPixelType;
    typedef typename DestImageType::value_type DestVectorType;
    typedef typename DestVectorType::value_type DestPixelType;
    typedef typename NumericTraits<DestPixelType>::isIntegral DestPixelIsIntegral;

    typedef typename PyramidVectorType::const_iterator PyramidVectorIterator;
    typedef typename DestVectorType::iterator DestVectorIterator;

    typedef typename PyramidImageType::const_traverser SrcTraverser;
    typedef typename DestImageType::traverser DestTraverser;
    typedef typename MaskImageType::const_traverser MaskTraverser;

    typedef Lab2RGBPrimeFunctor<double> ColorFunctor;
    typedef typename ColorFunctor::argument_type ColorFunctorArgumentType;
    typedef typename ColorFunctor::result_type ColorFunctorResultType;
    typedef typename ColorFunctorResultType::value_type ColorFunctorResultComponent;

    ColorFunctor cf(NumericTraits<DestPixelType>::max());

    if (Verbose > 0) {
        cout << "L*a*b* to R'G'B' color space conversion..." << endl;
    }

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;
    MaskTraverser my = mask_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y, ++my.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;
        MaskTraverser mx = my;

        for (; sx.x != send.x; ++sx.x, ++dx.x, ++mx.x) {
            if (ma(mx)) {
                // Convert from fixed point to floating point
                PyramidVectorType p = sa(sx);
                double l = convertFromPyramidMath<PyramidPixelType>(p.red());
                double a = convertFromPyramidMath<PyramidPixelType>(p.green());
                double b = convertFromPyramidMath<PyramidPixelType>(p.blue());
                ColorFunctorArgumentType labVector(l, a, b);

                // Convert from L*a*b* color space back to gamma-corrected RGB.
                ColorFunctorResultType rgbpVector = cf(labVector);

                // Convert from floating point to DestVectorType.
                // If DestPixelType is integral, use dithering.
                rgbpVector.setRed(dither(rgbpVector.red(), DestPixelIsIntegral()));
                rgbpVector.setGreen(dither(rgbpVector.green(), DestPixelIsIntegral()));
                rgbpVector.setBlue(dither(rgbpVector.blue(), DestPixelIsIntegral()));
                da.set(NumericTraits<DestVectorType>::fromRealPromote(rgbpVector), dx);

                //r = f(sa(sx));
                //PyramidVectorType p = sa(sx);
                //ColorArgType c;
                //c[0] = convertFromPyramidMath<DestPixelType, PyramidPixelType>(p.red());
                //c[1] = convertFromPyramidMath<DestPixelType, PyramidPixelType>(p.green());
                //c[2] = convertFromPyramidMath<DestPixelType, PyramidPixelType>(p.blue());
                //r = f(c);
                //PyramidVectorIterator pvi = sa(sx).begin();
                //DestVectorIterator dvi = r.begin();
                //for (; dvi != r.end(); ++dvi, ++pvi) {
                //    *dvi = convertFromPyramidMath<DestPixelType, PyramidPixelType>(*pvi);
                //}

                //da.set(r, dx);
            }
        }
    }

};

// Switch based on vector or scalar image types.
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        typename MaskImageType::const_traverser mask_upperleft,
        typename MaskImageType::ConstAccessor ma) {

    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar src_is_scalar;

    copyFromPyramidImageIf<DestImageType, PyramidImageType, MaskImageType>(
            src_upperleft,
            src_lowerright,
            sa,
            dest_upperleft,
            da,
            mask_upperleft,
            ma,
            src_is_scalar());
};

// Using argument object factory
template <typename DestImageType, typename PyramidImageType, typename MaskImageType>
inline void copyFromPyramidImageIf(
        triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
        pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest,
        pair<typename MaskImageType::const_traverser, typename MaskImageType::ConstAccessor> mask) {
    copyFromPyramidImageIf<DestImageType, PyramidImageType, MaskImageType>(
            src.first,
            src.second,
            src.third,
            dest.first,
            dest.second,
            mask.first,
            mask.second);
};

} // namespace enblend

#endif /* __FIXMATH_H__ */
