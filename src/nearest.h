/*
 * Copyright (C) 2009 Christoph L. Spiel
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
#ifndef __NEAREST_H__
#define __NEAREST_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#ifdef _WIN32
#include <cmath>
#else
#include <math.h>
#endif
#include <stdlib.h>
#include <utility>

#include "vigra/distancetransform.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdcachedfileimage.hxx"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;

using vigra::NumericTraits;
using vigra::triple;

namespace enblend {

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void
quadruple_image(SrcImageIterator src_upperleft,
                SrcImageIterator src_lowerright, SrcAccessor sa,
                DestImageIterator dest_upperleft, DestAccessor da,
                boundary_t boundary)
{
    const vigra::Diff2D size_x(src_lowerright.x - src_upperleft.x, 0);
    const vigra::Diff2D size_y(0, src_lowerright.y - src_upperleft.y);

    switch (boundary)
    {
    case OpenBoundaries:
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft, da);
        break;

    case HorizontalStrip:
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft, da); // 11
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft + size_x, da); // 12
        break;

    case VerticalStrip:
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft, da); // 11
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft + size_y, da); // 21
        break;

    case DoubleStrip:
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft, da); // 11
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft + size_x, da); // 12
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft + size_y, da); // 21
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft + size_x + size_y, da); // 22
        break;

    default:
        assert(false);
    }
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void
quadruple_image(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                vigra::pair<DestImageIterator, DestAccessor> dest,
                boundary_t boundary)
{
    quadruple_image(src.first, src.second, src.third,
                    dest.first, dest.second,
                    boundary);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void
quater_image(SrcImageIterator src_upperleft,
             SrcImageIterator src_lowerright, SrcAccessor sa,
             DestImageIterator dest_upperleft, DestAccessor da,
             boundary_t boundary)
{
    const vigra::Diff2D size_x(src_lowerright.x - src_upperleft.x, 0);
    const vigra::Diff2D size_y(0, src_lowerright.y - src_upperleft.y);
    const vigra::Diff2D size_x2(size_x / 2);
    const vigra::Diff2D size_y2(size_y / 2);
    const vigra::Diff2D size_x4(size_x2 / 2);
    const vigra::Diff2D size_y4(size_y2 / 2);

    // destination image
    //  | 11  12 |
    //  |        |
    //  | 21  22 |

    switch (boundary)
    {
    case OpenBoundaries:
        copyImage(src_upperleft, src_lowerright, sa,
                  dest_upperleft, da);
        break;

    case HorizontalStrip:
        copyImage(src_upperleft + size_x2,
                  src_upperleft + size_x2 + size_x4 + size_y,
                  sa,
                  dest_upperleft,
                  da); // 11
        copyImage(src_upperleft + size_x4,
                  src_upperleft + size_x2 + size_y,
                  sa,
                  dest_upperleft + size_x4,
                  da); // 12
        break;

    case VerticalStrip:
        copyImage(src_upperleft + size_y2,
                  src_upperleft + size_y2 + size_y4 + size_x,
                  sa,
                  dest_upperleft,
                  da); // 21
        copyImage(src_upperleft + size_y4,
                  src_upperleft + size_y2 + size_x,
                  sa,
                  dest_upperleft + size_y4,
                  da); // 22
        break;

    case DoubleStrip:
        copyImage(src_upperleft + size_x2 + size_y2,
                  src_upperleft + size_x2 + size_y2 + size_x4 + size_y4,
                  sa,
                  dest_upperleft,
                  da); // 11
        copyImage(src_upperleft + size_x4 + size_y2,
                  src_upperleft + size_x2 + size_y2 + size_y4,
                  sa,
                  dest_upperleft + size_x4,
                  da); // 12
        copyImage(src_upperleft + size_x2 + size_y4,
                  src_upperleft + size_x2 + size_x4 + size_y2,
                  sa,
                  dest_upperleft + size_y4,
                  da); // 21
        copyImage(src_upperleft + size_x4 + size_y4,
                  src_upperleft + size_x2 + size_y2,
                  sa,
                  dest_upperleft + size_x4 + size_y4,
                  da); // 22
        break;

    default:
        assert(false);
    }
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void
quater_image(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
             vigra::pair<DestImageIterator, DestAccessor> dest,
             boundary_t boundary)
{
    quater_image(src.first, src.second, src.third,
                 dest.first, dest.second,
                 boundary);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class ValueType>
void
periodicDistanceTransform(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor sa,
                          DestImageIterator dest_upperleft, DestAccessor da,
                          ValueType background, int norm, boundary_t boundary)
{
    typedef typename SrcImageIterator::value_type SrcValueType;
    typedef typename DestImageIterator::value_type DestValueType;

    const vigra::Diff2D size(src_lowerright.x - src_upperleft.x,
                             src_lowerright.y - src_upperleft.y);
    int size_x;
    int size_y;

    switch (boundary)
    {
    case OpenBoundaries:
        size_x = size.x;
        size_y = size.y;
        break;

    case HorizontalStrip:
        size_x = 2 * size.x;
        size_y = size.y;
        break;

    case VerticalStrip:
        size_x = size.x;
        size_y = 2 * size.y;
        break;

    case DoubleStrip:
        size_x = 2 * size.x;
        size_y = 2 * size.y;
        break;

    default:
        size_x = -1;
        size_y = -1;
        assert(false);
    }

    vigra::BasicImage<SrcValueType> periodic(size_x, size_y);
    vigra::BasicImage<DestValueType> distance(periodic.size());

    quadruple_image(src_upperleft, src_lowerright, sa,
                    periodic.upperLeft(), periodic.accessor(),
                    boundary);
    vigra::distanceTransform(srcImageRange(periodic), destImage(distance),
                             background, norm);
    quater_image(srcImageRange(distance), destIter(dest_upperleft, da), boundary);
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor,
          class ValueType>
inline void
periodicDistanceTransform(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                          vigra::pair<DestImageIterator, DestAccessor> dest,
                          ValueType background, int norm, boundary_t boundary)
{
    periodicDistanceTransform(src.first, src.second, src.third,
                              dest.first, dest.second,
                              background, norm, boundary);
}


template <class ValueType>
struct saturating_subtract
{
    typedef ValueType first_argument_type;
    typedef ValueType second_argument_type;
    typedef ValueType result_type;

    result_type operator()(const first_argument_type& v1,
                           const second_argument_type& v2) const
    {
        typedef vigra::NumericTraits<result_type> traits;

        return v2 < v1 ? v1 - v2 : traits::zero();
    }
};


// Compute a mask (dest) that defines the seam line given the
// blackmask (src1) and the whitemask (src2) of the overlapping
// images.
//
// The idea of the algorithm is from
//     Yalin Xiong, Ken Turkowski
//     "Registration, Calibration and Blending in Creating High Quality Panoramas"
//     Proceedings of the 4th IEEE Workshop on Applications of Computer Vision (WACV'98)
// where we find:
//     "To locate the mask boundary, we perform the grassfire
//      transform on two images individually.  The resulting distance
//      maps represent how far away each pixel is from its nearest
//      boundary.  The pixel values of the blend mask is then set to
//      either 0 or 1 by comparing the distance values at each pixel
//      in the two distance maps."
//
// Though we prefer the Distance Transform to the Grassfire Transform.

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void
nearestFeatureTransform(SrcImageIterator src1_upperleft, SrcImageIterator src1_lowerright, SrcAccessor sa1,
                        SrcImageIterator src2_upperleft, SrcAccessor sa2,
                        DestImageIterator dest_upperleft, DestAccessor da,
                        nearest_neigbor_metric_t norm, boundary_t boundary)
{
    typedef typename SrcAccessor::value_type SrcPixelType;
    typedef vigra::NumericTraits<SrcPixelType> SrcPixelTraits;
    typedef typename SrcPixelTraits::Promote SrcPromoteType;

    typedef typename DestAccessor::value_type DestPixelType;
    typedef vigra::NumericTraits<DestPixelType> DestPixelTraits;

    const SrcPixelType background = SrcPixelTraits::zero();
    const Diff2D size(src1_lowerright.x - src1_upperleft.x,
                      src1_lowerright.y - src1_upperleft.y);

    IMAGETYPE<SrcPromoteType> dist12(size);
    IMAGETYPE<SrcPromoteType> dist21(size);
    if (Verbose > VERBOSE_NFT_MESSAGES)
    {
        cerr << command << ": info: creating blend mask: 1/3";
        cerr.flush();
    }

#ifdef OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef OPENMP
#pragma omp section
#endif
        {
            IMAGETYPE<SrcPixelType> diff12(size);
            combineTwoImages(src1_upperleft, src1_lowerright, sa1,
                             src2_upperleft, sa2,
                             diff12.upperLeft(), diff12.accessor(),
                             saturating_subtract<SrcPixelType>());
            switch (boundary)
            {
            case OpenBoundaries:
                distanceTransform(srcImageRange(diff12), destImage(dist12),
                                  background, norm);
                break;

            case HorizontalStrip: // FALLTHROUGH
            case VerticalStrip:   // FALLTHROUGH
            case DoubleStrip:
                periodicDistanceTransform(srcImageRange(diff12), destImage(dist12),
                                          background, norm, boundary);
                break;

            default:
                assert(false);
            }
        } // omp section

#ifdef OPENMP
#pragma omp section
#endif
        {
            if (Verbose > VERBOSE_NFT_MESSAGES)
            {
                cerr << " 2/3";
                cerr.flush();
            }
            IMAGETYPE<SrcPixelType> diff21(size);
            combineTwoImages(src2_upperleft, src2_upperleft + size, sa2,
                             src1_upperleft, sa1,
                             diff21.upperLeft(), diff21.accessor(),
                             saturating_subtract<SrcPixelType>());
            switch (boundary)
            {
            case OpenBoundaries:
                distanceTransform(srcImageRange(diff21), destImage(dist21),
                                  background, norm);
                break;

            case HorizontalStrip: // FALLTHROUGH
            case VerticalStrip:   // FALLTHROUGH
            case DoubleStrip:
                periodicDistanceTransform(srcImageRange(diff21), destImage(dist21),
                                          background, norm, boundary);
                break;

            default:
                assert(false);
            }
        } // omp section
    } // omp parallel sections

    if (Verbose > VERBOSE_NFT_MESSAGES)
    {
        cerr << " 3/3";
        cerr.flush();
    }
    combineTwoImagesMP(dist12.upperLeft(), dist12.lowerRight(), dist12.accessor(),
                       dist21.upperLeft(), dist21.accessor(),
                       dest_upperleft, da,
                       ifThenElse(vigra::functor::Arg1() < vigra::functor::Arg2(),
                                  vigra::functor::Param(DestPixelTraits::max()),
                                  vigra::functor::Param(DestPixelTraits::zero())));

    if (Verbose > VERBOSE_NFT_MESSAGES)
    {
        cerr << endl;
    }
}


template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void
nearestFeatureTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src1,
                        pair<SrcImageIterator, SrcAccessor> src2,
                        pair<DestImageIterator, DestAccessor> dest,
                        nearest_neigbor_metric_t norm, boundary_t boundary)
{
    nearestFeatureTransform(src1.first, src1.second, src1.third,
                            src2.first, src2.second,
                            dest.first, dest.second,
                            norm, boundary);
}

} // namespace enblend

#endif /* __NEAREST_H__ */

// Local Variables:
// mode: c++
// End:
