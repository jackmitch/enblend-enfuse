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
                PyramidVectorType p = sa(sx);
                double l = convertFromPyramidMath<PyramidPixelType>(p.red());
                double a = convertFromPyramidMath<PyramidPixelType>(p.green());
                double b = convertFromPyramidMath<PyramidPixelType>(p.blue());
                ColorFunctorArgumentType labVector(l, a, b);
                ColorFunctorResultType rgbpVector = cf(labVector);
                DestVectorType r = NumericTraits<DestVectorType>::fromRealPromote(rgbpVector);
                da.set(r, dx);

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

template <typename DestImageType, typename PyramidImageType>
inline void copyLabPyramid(
        triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
        pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest) {

};

} // namespace enblend

#endif /* __FIXMATH_H__ */
