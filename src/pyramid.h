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
#ifndef __PYRAMID_H__
#define __PYRAMID_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <boost/static_assert.hpp>

#include "vigra/numerictraits.hxx"

using std::vector;
using vigra::DRGBImage;
using vigra::NumericTraits;

namespace enblend {

static const double A = 0.4;
static const double W[] = {0.25 - A / 2.0, 0.25, A, 0.25, 0.25 - A / 2.0};
static const unsigned int A100 = 40;
static const unsigned int W100[] = {25 - A100 / 2, 25, A100, 25, 25 - A100 / 2};

/** Calculate the half width of a level n filter, taking into account
 *  pixel precision and rounding method.
 */
template <typename PixelType>
unsigned int filterHalfWidth(const unsigned int level) {

    // This is the arithmetic half width (true for level > 0).
    unsigned int length = 1 + (1 << level);

    PixelType *f = new PixelType[length];
    for (unsigned int i = 0; i < length; i++) {
        f[i] = NumericTraits<PixelType>::zero();
    }

    // input f(x) is the step function u(-x)
    f[0] = NumericTraits<PixelType>::max();
    double maxPixelValue = NumericTraits<PixelType>::toRealPromote(f[0]);

    for (unsigned int l = 1; l <= level; l++) {
        // sample 0 from level l-1
        double pZero = NumericTraits<PixelType>::toRealPromote(f[0]);

        // Sample 1 from level l-1
        double pOne = NumericTraits<PixelType>::toRealPromote(f[1 << (l-1)]);

        // Sample 0 on level l
        double nZero = (pZero * W[2]) + (pOne * W[3])
                + (maxPixelValue * W[0])
                + (maxPixelValue * W[1]);
        f[0] = NumericTraits<PixelType>::fromRealPromote(nZero);

        // Sample 1 on level l
        double nOne = (pZero * W[0]) + (pOne * W[1]);
        f[1 << l] = NumericTraits<PixelType>::fromRealPromote(nOne);

        // Remaining samples on level l are zero.

        // If sample 1 was rounded down to zero, then sample 1 on
        // level l-1 is the rightmost nonzero value.
        if (f[1 << l] == NumericTraits<PixelType>::zero()) {
            delete[] f;
            // return the index of the rightmost nonzero value.
            return (1 << (l-1));
        }
    }

    // Else there is no round-to-zero issue.
    delete[] f;
    return (length - 1);
}

//template <typename SrcImageType, typename PyramidImageType>
//vector<int*> *gaussianPyramid(unsigned int numLevels,
//        SrcImageIterator src_upperleft,
//        SrcImageIterator src_lowerright,
//        SrcAccessor sa) {
//    //FIXME
//    return new vector<int*>();
//}

template <typename SrcImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::const_accessor sa,
        VigraFalseType) {

    // ImageTypes are both vectors.
}

template <typename SrcImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::const_accessor sa,
        VigraTrueType) {

    // ImageTypes are both scalars.

}
template <typename SrcImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::const_accessor sa) {

    // This indicates if the source image is a vector or a scalar.
    typedef typename NumericTraits<typename SrcImageType::value_type>::isScalar src_is_scalar;
    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar pyramid_is_scalar;

    // The source image and the pyramid image must both be vectors or both be scalars.
    BOOST_STATIC_ASSERT(src_is_scalar() == pyramid_is_scalar());

    return gaussianPyramid<SrcImageType, PyramidImageType>(numLevels,
            src_upperleft,
            src_lowerright,
            sa,
            pyramid_is_scalar());
}

} // namespace enblend

#endif /* __PYRAMID_H__ */
