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
#include "vigra/numerictraits.hxx"
#include "vigra/utilities.hxx"

using std::pair;

using vigra::NumericTraits;
using vigra::triple;
using vigra::VigraFalseType;
using vigra::VigraTrueType;

namespace enblend {

// Convert type T1 to fixed-point value stored in type T2.
template <typename T1, typename T2>
inline T2 convertToPyramidMath(T1 v) {
    // Fixed-point math conversion not defined for this type pair.
    BOOST_STATIC_ASSERT(false);
};

// Convert to type T1 from fixed-point value stored in type T2.
template <typename T1, typename T2>
inline T1 convertFromPyramidMath(T2 v) {
    // Fixed-point math conversion not defined for this type pair.
    BOOST_STATIC_ASSERT(false);
};

template <>
inline short convertToPyramidMath<unsigned char, short>(unsigned char v) {
    // cast to short does not sign-extend
    // Shift left 7 bits
    return (short)v << 7;
};

template <>
inline unsigned char convertFromPyramidMath<unsigned char, short>(short v) {
    double r = (double)v / 128.0;
    return NumericTraits<unsigned char>::fromRealPromote(r);
};

template <>
inline int convertToPyramidMath<short, int>(short v) {
    // cast to int sign-extends
    // Shift left 15 bits
    return (int)v << 15;
};

template <>
inline short convertFromPyramidMath<short, int>(int v) {
    double r = (double)v / 32768.0;
    return NumericTraits<short>::fromRealPromote(r);
};

template <>
inline int convertToPyramidMath<unsigned short, int>(unsigned short v) {
    // cast to int does not sign-extend
    // Shift left 15 bits
    return (int)v << 15;
};

template <>
inline unsigned short convertFromPyramidMath<unsigned short, int>(int v) {
    double r = (double)v / 32768.0;
    return NumericTraits<unsigned short>::fromRealPromote(r);
};

template <>
inline double convertToPyramidMath<int, double>(int v) {
    return NumericTraits<int>::toRealPromote(v);
};

template <>
inline int convertFromPyramidMath<int, double>(double v) {
    return NumericTraits<int>::fromRealPromote(v);
};

template <>
inline double convertToPyramidMath<unsigned int, double>(unsigned int v) {
    return NumericTraits<unsigned int>::toRealPromote(v);
};

template <>
inline unsigned int convertFromPyramidMath<unsigned int, double>(double v) {
    return NumericTraits<unsigned int>::fromRealPromote(v);
};

template <>
inline double convertToPyramidMath<float, double>(float v) {
    return (double)v;
};

template <>
inline float convertFromPyramidMath<float, double>(double v) {
    return (float)v;
};

template <>
inline double convertToPyramidMath<double, double>(double v) {
    return v;
};

template <>
inline double convertFromPyramidMath<double, double>(double v) {
    return v;
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

    typedef typename SrcVectorType::iterator SrcVectorIterator;
    typedef typename PyramidVectorType::iterator PyramidVectorIterator;

    typedef typename SrcImageType::const_traverser SrcTraverser;
    typedef typename PyramidImageType::traverser DestTraverser;

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;

        for (; sx.x != send.x; ++sx.x, ++dx.x) {
            PyramidVectorType r;

            SrcVectorIterator svi = sa(sx).begin();
            PyramidVectorIterator pvi = r.begin();
            for (; pvi != r.end(); ++pvi, ++svi) {
                *pvi = convertToPyramidMath<SrcPixelType, PyramidPixelType>(*svi);
            }

            da.set(r, dx);
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
template <typename DestImageType, typename PyramidImageType>
inline void copyFromPyramidImage(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        VigraTrueType) {

    typedef typename PyramidImageType::value_type PyramidPixelType;
    typedef typename DestImageType::value_type DestPixelType;

    typedef typename PyramidImageType::const_traverser SrcTraverser;
    typedef typename DestImageType::traverser DestTraverser;

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;

        for (; sx.x != send.x; ++sx.x, ++dx.x) {
            da.set(convertFromPyramidMath<DestPixelType, PyramidPixelType>(sa(sx)), dx);
        }
    }

};

// copy vector fixed-point pyramid image to vector image.
template <typename DestImageType, typename PyramidImageType>
inline void copyFromPyramidImage(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da,
        VigraFalseType) {

    typedef typename PyramidImageType::value_type PyramidVectorType;
    typedef typename PyramidVectorType::value_type PyramidPixelType;
    typedef typename DestImageType::value_type DestVectorType;
    typedef typename DestVectorType::value_type DestPixelType;

    typedef typename PyramidVectorType::iterator PyramidVectorIterator;
    typedef typename DestVectorType::iterator DestVectorIterator;

    typedef typename PyramidImageType::const_traverser SrcTraverser;
    typedef typename DestImageType::traverser DestTraverser;

    SrcTraverser sy = src_upperleft;
    SrcTraverser send = src_lowerright;
    DestTraverser dy = dest_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y) {
        SrcTraverser sx = sy;
        DestTraverser dx = dy;

        for (; sx.x != send.x; ++sx.x, ++dx.x) {
            DestVectorType r;

            PyramidVectorIterator pvi = sa(sx).begin();
            DestVectorIterator dvi = r.begin();
            for (; dvi != r.end(); ++dvi, ++pvi) {
                *dvi = convertFromPyramidMath<DestPixelType, PyramidPixelType>(*pvi);
            }

            da.set(r, dx);
        }
    }

};

// Switch based on vector or scalar image types.
template <typename DestImageType, typename PyramidImageType>
inline void copyFromPyramidImage(
        typename PyramidImageType::const_traverser src_upperleft,
        typename PyramidImageType::const_traverser src_lowerright,
        typename PyramidImageType::ConstAccessor sa,
        typename DestImageType::traverser dest_upperleft,
        typename DestImageType::Accessor da) {

    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar src_is_scalar;

    copyFromPyramidImage<DestImageType, PyramidImageType>(
            src_upperleft,
            src_lowerright,
            sa,
            dest_upperleft,
            da,
            src_is_scalar());
};

// Using argument object factory
template <typename DestImageType, typename PyramidImageType>
inline void copyFromPyramidImage(
        triple<typename PyramidImageType::const_traverser, typename PyramidImageType::const_traverser, typename PyramidImageType::ConstAccessor> src,
        pair<typename DestImageType::traverser, typename DestImageType::Accessor> dest) {
    copyFromPyramidImage<DestImageType, PyramidImageType>(
            src.first,
            src.second,
            src.third,
            dest.first,
            dest.second);
};

} // namespace enblend

#endif /* __FIXMATH_H__ */
