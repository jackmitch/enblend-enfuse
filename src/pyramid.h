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

#include <assert.h>
#include <vector>
//#include <boost/static_assert.hpp>

#include "fixmath.h"

#include "vigra/numerictraits.hxx"

using std::cout;
using std::vector;
using vigra::NumericTraits;
using vigra::triple;

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

template <typename SrcImageIterator, typename SrcAccessor,
        typename MaskIterator, typename MaskAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename DestMaskIterator, typename DestMaskAccessor>
inline void reduce(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        MaskIterator mask_upperleft,
        MaskAccessor ma,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da,
        DestMaskIterator dest_mask_upperleft,
        DestMaskIterator dest_mask_lowerright,
        DestMaskAccessor dma) {

    typedef typename SrcAccessor::value_type PixelType;
    typedef typename NumericTraits<PixelType>::RealPromote RealPixelType;
    typedef typename DestMaskAccessor::value_type DestMaskPixelType;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    assert(src_w > 1 && src_h > 1);

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    SrcImageIterator sy = src_upperleft;
    MaskIterator my = mask_upperleft;
    DestMaskIterator dmy = dest_mask_upperleft;
    for (int srcy = 0; dy.y != dend.y; ++dy.y, ++dmy.y, sy.y+=2, my.y+=2, srcy+=2) {

        DestImageIterator dx = dy;
        SrcImageIterator sx = sy;
        MaskIterator mx = my;
        DestMaskIterator dmx = dmy;
        for (int srcx = 0; dx.x != dend.x; ++dx.x, ++dmx.x, sx.x+=2, mx.x+=2, srcx+=2) {

            RealPixelType p(NumericTraits<RealPixelType>::zero());
            unsigned int noContrib = 10000;

            for (int kx = -2; kx <= 2; kx++) {
                int bounded_kx = kx;

                if (wraparound) {
                    // Boundary condition: wrap around the image.
                    if (srcx + kx < 0) bounded_kx += src_w;
                    if (srcx + kx >= src_w) bounded_kx -= src_w;
                } else {
                    // Boundary condition: replicate first and last column.
                    if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                    if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                }

                for (int ky = -2; ky <= 2; ky++) {
                    int bounded_ky = ky;

                    // Boundary condition: replicate top and bottom rows.
                    if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                    if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                    if (mx(bounded_kx, bounded_ky)) {
                        p += (W[kx+2] * W[ky+2]) * sx(bounded_kx, bounded_ky);
                    } else {
                        // Transparent pixels don't count.
                        noContrib -= W100[kx+2] * W100[ky+2];
                    }

                }
            }

            // Adjust filter for any ignored transparent pixels.
            if (noContrib != 0) p /= ((double)noContrib / 10000.0);

            da.set(NumericTraits<PixelType>::fromRealPromote(p), dx);
            dma.set((noContrib == 0) ? NumericTraits<DestMaskPixelType>::zero()
                                     : NumericTraits<DestMaskPixelType>::nonZero(),
                    dmx);

        }
    }
};

template <typename SrcImageIterator, typename SrcAccessor,
        typename MaskIterator, typename MaskAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename DestMaskIterator, typename DestMaskAccessor>
inline void reduce(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<MaskIterator, MaskAccessor> mask,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest,
        triple<DestMaskIterator, DestMaskIterator, DestMaskAccessor> destMask) {
    reduce(wraparound,
            src.first, src.second, src.third,
            mask.first, mask.second,
            dest.first, dest.second, dest.third,
            destMask.first, destMask.second, destMask.third);
};

template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void reduce(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename SrcAccessor::value_type PixelType;
    typedef typename NumericTraits<PixelType>::RealPromote RealPixelType;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    assert(src_w > 1 && src_h > 1);

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    SrcImageIterator sy = src_upperleft;
    for (int srcy = 0; dy.y != dend.y; ++dy.y, sy.y+=2, srcy+=2) {

        DestImageIterator dx = dy;
        SrcImageIterator sx = sy;
        for (int srcx = 0; dx.x != dend.x; ++dx.x, sx.x+=2, srcx+=2) {

            RealPixelType p(NumericTraits<RealPixelType>::zero());

            for (int kx = -2; kx <= 2; kx++) {
                int bounded_kx = kx;

                if (wraparound) {
                    // Boundary condition: wrap around the image.
                    if (srcx + kx < 0) bounded_kx += src_w;
                    if (srcx + kx >= src_w) bounded_kx -= src_w;
                } else {
                    // Boundary condition: replicate first and last column.
                    if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                    if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                }

                for (int ky = -2; ky <= 2; ky++) {
                    int bounded_ky = ky;

                    // Boundary condition: replicate top and bottom rows.
                    if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                    if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                    p += (W[kx+2] * W[ky+2]) * sx(bounded_kx, bounded_ky);

                }
            }

            da.set(NumericTraits<PixelType>::fromRealPromote(p), dx);
        }
    }

};

template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void reduce(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest) {
    reduce(wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, dest.third);
};

template <typename SrcImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa) {

    vector<PyramidImageType*> *gp = new vector<PyramidImageType*>();

    // Size of pyramid level 0
    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    // Pyramid level 0
    PyramidImageType *gp0 = new PyramidImageType(w, h);

    if (Verbose > 0) {
        cout << "Generating Gaussian pyramid: g0";
        cout.flush();
    }

    // Copy src image into gp0, using fixed-point conversions.
    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src_upperleft, src_lowerright, sa,
            gp0->upperLeft(), gp0->accessor());

    gp->push_back(gp0);

    // Make remaining levels.
    PyramidImageType *lastGP = gp0;
    for (unsigned int l = 1; l < numLevels; l++) {

        if (Verbose > 0) {
            cout << " g" << l;
            cout.flush();
        }

        // Size of next level
        w = (w + 1) >> 1;
        h = (h + 1) >> 1;

        // Next pyramid level
        PyramidImageType *gpn = new PyramidImageType(w, h);

        reduce(wraparound, srcImageRange(*lastGP), destImageRange(*gpn));

        gp->push_back(gpn);
        lastGP = gpn;
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return gp;
};

template <typename SrcImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src) {
    return gaussianPyramid<SrcImageType, PyramidImageType>(numLevels,
            wraparound,
            src.first, src.second, src.third);
};

} // namespace enblend

#endif /* __PYRAMID_H__ */
