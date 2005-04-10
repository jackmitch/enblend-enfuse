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

#include "vigra/convolution.hxx"
#include "vigra/error.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"

#include "fixmath.h"

using std::cout;
using std::vector;

using vigra::linearRangeMapping;
using vigra::NumericTraits;
using vigra::transformImage;
using vigra::triple;
using vigra::USImage;
using vigra::USRGBImage;

namespace enblend {

// Pyramid filter coefficients.
static const double A = 0.4;
static const double W[] = {0.25 - A / 2.0, 0.25, A, 0.25, 0.25 - A / 2.0};
static const unsigned int A100 = 40;
static const unsigned int W100[] = {25 - A100 / 2, 25, A100, 25, 25 - A100 / 2};

/** Calculate the half-width of a n-level filter.
 *  Assumes that the input function is a left-handed function,
 *  and the last non-zero input is at location 0.
 *  Returns the location of the last non-zero output.
 */
template <typename PixelType>
unsigned int filterHalfWidth(const unsigned int levels) {
    // For levels >= 30, the full width will just barely fit in int32.
    // When this range is added to a bounding box it will certainly
    // overflow the Diff2D.
    vigra_precondition((levels >= 1 && levels <= 29),
            "filterHalfWidth: levels outside of range [1,29]");

    // This is the arithmetic half width.
    int halfWidth = (levels == 1) ? 0 : ((1 << (levels+1)) - 4);

    return halfWidth;
}

/** The Burt & Adelson Reduce operation.
 *  This version is for images with alpha channels.
 */
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
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in reduce");

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

// Version using argument object factories.
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

/** The Burt & Adelson Reduce operation.
 *  This version is for images that do not have alpha channels.
 */
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
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in reduce");

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

// Version using argument object factories.
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void reduce(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest) {
    reduce(wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, dest.third);
};

/** The Burt & Adelson Expand operation. */
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
void expand(bool add, bool wraparound,
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

    DestImageIterator dy = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    SrcImageIterator sy = src_upperleft;
    for (int srcy = 0, desty = 0; dy.y != dend.y; ++dy.y, desty++) {

        DestImageIterator dx = dy;
        SrcImageIterator sx = sy;
        for (int srcx = 0, destx = 0; dx.x != dend.x; ++dx.x, destx++) {

            RealPixelType p(NumericTraits<RealPixelType>::zero());
            unsigned int totalContrib = 0;

            if ((destx & 1) == 1) {
                // This is an odd-numbered column.
                for (int kx = 0; kx <= 1; kx++) {
                    // kx=0 -> W[-1]
                    // kx=1 -> W[ 1]
                    int wIndexA = (kx==0) ? 1 : 3;
                    int bounded_kx = kx;

                    if (wraparound) {
                        // Boundary condition - wrap around the image.
                        // First boundary case cannot happen when destx is odd.
                        //if (srcx + kx < 0) bounded_kx += src_w;
                        if (srcx + kx >= src_w) bounded_kx -= src_w;
                    } else {
                        // Boundary condition - replicate first and last column.
                        // First boundary case cannot happen when destx is odd.
                        //if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                        if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                    }

                    if ((desty & 1) == 1) {
                        // This is an odd-numbered row.
                        for (int ky = 0; ky <= 1; ky++) {
                            // ky=0 -> W[-1]
                            // ky=1 -> W[ 1]
                            int wIndexB = (ky==0) ? 1 : 3;
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            // First boundary condition cannot happen when desty is odd.
                            //if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                    else {
                        // This is an even-numbered row.
                        for (int ky = -1; ky <= 1; ky++) {
                            // ky=-1 -> W[-2]
                            // ky= 0 -> W[ 0]
                            // ky= 1 -> W[ 2]
                            int wIndexB = (ky==-1) ? 0 : ((ky==0) ? 2 : 4);
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                }
            }
            else {
                // This is an even-numbered column.
                for (int kx = -1; kx <= 1; kx++) {
                    // kx=-1 -> W[-2]
                    // kx= 0 -> W[ 0]
                    // kx= 1 -> W[ 2]
                    int wIndexA = (kx==-1) ? 0 : ((kx==0) ? 2 : 4); 
                    int bounded_kx = kx;

                    if (wraparound) {
                        // Boundary condition - wrap around the image.
                        if (srcx + kx < 0) bounded_kx += src_w;
                        if (srcx + kx >= src_w) bounded_kx -= src_w;
                    } else {
                        // Boundary condition - replicate first and last column.
                        if (srcx + kx < 0) bounded_kx -= (srcx + kx);
                        if (srcx + kx >= src_w) bounded_kx -= (srcx + kx - (src_w - 1));
                    }

                    if ((desty & 1) == 1) {
                        // This is an odd-numbered row.
                        for (int ky = 0; ky <= 1; ky++) {
                            // ky=0 -> W[-1]
                            // ky=1 -> W[ 1]
                            int wIndexB = (ky==0) ? 1 : 3;
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            // First boundary condition cannot happen when desty is odd.
                            //if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                    else {
                        // This is an even-numbered row.
                        for (int ky = -1; ky <= 1; ky++) {
                            // ky=-1 -> W[-2]
                            // ky= 0 -> W[ 0]
                            // ky= 1 -> W[ 2]
                            int wIndexB = (ky==-1) ? 0 : ((ky==0) ? 2 : 4);
                            int bounded_ky = ky;

                            // Boundary condition: replicate top and bottom rows.
                            if (srcy + ky < 0) bounded_ky -= (srcy + ky);
                            if (srcy + ky >= src_h) bounded_ky -= (srcy + ky - (src_h - 1));

                            p += (W[wIndexA] * W[wIndexB]) * sx(bounded_kx, bounded_ky);
                            totalContrib += (W100[wIndexA] * W100[wIndexB]);
                        }
                    }
                }
            }

            p /= (totalContrib / 10000.0);

            // Get dest pixel at dx.
            RealPixelType pOrig = NumericTraits<PixelType>::toRealPromote(da(dx));
            RealPixelType pMod;

            // Add or subtract p.
            if (add) {
                pMod = pOrig + p;
            }
            else {
                pMod = pOrig - p;
            }

            // Write back the result.
            da.set(NumericTraits<PixelType>::fromRealPromote(pMod), dx);

            if ((destx & 1) == 1) {
                sx.x++;
                srcx++;
            }
        }

        if ((desty & 1) == 1) {
            sy.y++;
            srcy++;
        }
    }

};

// Version using argument object factories.
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void expand(bool add, bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        triple<DestImageIterator, DestImageIterator, DestAccessor> dest) {
    expand(add, wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, dest.third);
};

/** Calculate the Gaussian pyramid for the given SrcImage/AlphaImage pair. */
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename AlphaImageType::const_traverser alpha_upperleft,
        typename AlphaImageType::ConstAccessor aa) {

    vector<PyramidImageType*> *gp = new vector<PyramidImageType*>();

    // Size of pyramid level 0
    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    // Pyramid level 0
    PyramidImageType *gp0 = new PyramidImageType(w, h);

    // Copy src image into gp0, using fixed-point conversions.
    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src_upperleft, src_lowerright, sa,
            gp0->upperLeft(), gp0->accessor());

    gp->push_back(gp0);

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Generating Gaussian pyramid: g0";
    }

    // Make remaining levels.
    PyramidImageType *lastGP = gp0;
    AlphaImageType *lastA = NULL;
    for (unsigned int l = 1; l < numLevels; l++) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " g" << l;
            cout.flush();
        }

        // Size of next level
        w = (w + 1) >> 1;
        h = (h + 1) >> 1;

        // Next pyramid level
        PyramidImageType *gpn = new PyramidImageType(w, h);
        AlphaImageType *nextA = new AlphaImageType(w, h);

        if (lastA == NULL) {
            reduce(wraparound,
                    srcImageRange(*lastGP), maskIter(alpha_upperleft, aa),
                    destImageRange(*gpn), destImageRange(*nextA));
        } else {
            reduce(wraparound,
                    srcImageRange(*lastGP), maskImage(*lastA),
                    destImageRange(*gpn), destImageRange(*nextA));
        }

        gp->push_back(gpn);
        lastGP = gpn;
        delete lastA;
        lastA = nextA;
    }

    delete lastA;

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << endl;
    }

    return gp;

};

// Version using argument object factories.
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        pair<typename AlphaImageType::const_traverser, typename AlphaImageType::ConstAccessor> alpha) {
    return gaussianPyramid<SrcImageType, AlphaImageType, PyramidImageType>(
            numLevels, wraparound,
            src.first, src.second, src.third,
            alpha.first, alpha.second);
};

/** Calculate the Gaussian pyramid for the given image (without an alpha channel). */
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

    // Copy src image into gp0, using fixed-point conversions.
    copyToPyramidImage<SrcImageType, PyramidImageType>(
            src_upperleft, src_lowerright, sa,
            gp0->upperLeft(), gp0->accessor());

    gp->push_back(gp0);

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Generating Gaussian pyramid: g0";
    }

    // Make remaining levels.
    PyramidImageType *lastGP = gp0;
    for (unsigned int l = 1; l < numLevels; l++) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
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

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << endl;
    }

    return gp;
};

// Version using argument object factories.
template <typename SrcImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *gaussianPyramid(unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src) {
    return gaussianPyramid<SrcImageType, PyramidImageType>(numLevels,
            wraparound,
            src.first, src.second, src.third);
};

/** Calculate the Laplacian pyramid of the given SrcImage/AlphaImage pair. */
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
vector<PyramidImageType*> *laplacianPyramid(unsigned int numLevels,
        bool wraparound,
        typename SrcImageType::const_traverser src_upperleft,
        typename SrcImageType::const_traverser src_lowerright,
        typename SrcImageType::ConstAccessor sa,
        typename AlphaImageType::const_traverser alpha_upperleft,
        typename AlphaImageType::ConstAccessor aa) {

    // First create a Gaussian pyramid.
    vector <PyramidImageType*> *gp =
            gaussianPyramid<SrcImageType, AlphaImageType, PyramidImageType>(
                    numLevels, wraparound,
                    src_upperleft, src_lowerright, sa,
                    alpha_upperleft, aa);

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Generating Laplacian pyramid:";
        cout.flush();
    }

    // For each level, subtract the expansion of the next level.
    // Stop if there is no next level.
    for (unsigned int l = 0; l < (numLevels-1); l++) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " l" << l;
            cout.flush();
        }

        expand(false, wraparound,
                srcImageRange(*((*gp)[l+1])),
                destImageRange(*((*gp)[l])));
    }

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << " l" << (numLevels-1) << endl;
    }

    return gp;
};

// Version using argument object factories.
template <typename SrcImageType, typename AlphaImageType, typename PyramidImageType>
inline vector<PyramidImageType*> *laplacianPyramid(unsigned int numLevels,
        bool wraparound,
        triple<typename SrcImageType::const_traverser, typename SrcImageType::const_traverser, typename SrcImageType::ConstAccessor> src,
        pair<typename AlphaImageType::const_traverser, typename AlphaImageType::ConstAccessor> alpha) {
    return laplacianPyramid<SrcImageType, AlphaImageType, PyramidImageType>(
            numLevels, wraparound,
            src.first, src.second, src.third,
            alpha.first, alpha.second);
};

/** Collapse the given Laplacian pyramid. */
template <typename PyramidImageType>
void collapsePyramid(bool wraparound, vector<PyramidImageType*> *p) {

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << "Collapsing Laplacian pyramid: "
             << "l" << p->size()-1;
        cout.flush();
    }

    // For each level, add the expansion of the next level.
    // Work backwards from the smallest level to the largest.
    for (int l = (p->size()-2); l >= 0; l--) {

        if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
            cout << " l" << l;
            cout.flush();
        }

        expand(true, wraparound,
                srcImageRange(*((*p)[l+1])),
                destImageRange(*((*p)[l])));
    }

    if (Verbose > VERBOSE_PYRAMID_MESSAGES) {
        cout << endl;
    }

};

// Export a scalar pyramid as a set of UINT16 tiff files.
template <typename PyramidImageType>
void exportPyramid(vector<PyramidImageType*> *v, char *prefix, VigraTrueType) {
    typedef typename PyramidImageType::value_type PyramidValueType;

    //for (unsigned int i = 0; i < (v->size() - 1); i++) {
    //    // Clear all levels except last.
    //    initImage(destImageRange(*((*v)[i])), NumericTraits<PyramidValueType>::zero());
    //}
    //collapsePyramid(false, v);

    for (unsigned int i = 0; i < v->size(); i++) {
        char filenameBuf[512];
        snprintf(filenameBuf, 512, "%s%04u.tif", prefix, i);

        // Rescale the pyramid values to fit in UINT16.
        USImage usPyramid((*v)[i]->width(), (*v)[i]->height());
        transformImage(srcImageRange(*((*v)[i])), destImage(usPyramid),
                linearRangeMapping(NumericTraits<PyramidValueType>::min(),
                                     NumericTraits<PyramidValueType>::max(),
                                     NumericTraits<unsigned short>::min(),
                                     NumericTraits<unsigned short>::max()));

        ImageExportInfo info(filenameBuf);
        exportImage(srcImageRange(usPyramid), info);
    }
};

// Export a vector pyramid as a set of UINT16 tiff files.
template <typename PyramidImageType>
void exportPyramid(vector<PyramidImageType*> *v, char *prefix, VigraFalseType) {
    typedef typename PyramidImageType::value_type PyramidVectorType;
    typedef typename PyramidVectorType::value_type PyramidValueType;

    //for (unsigned int i = 0; i < (v->size() - 1); i++) {
    //    // Clear all levels except last.
    //    initImage(destImageRange(*((*v)[i])), NumericTraits<PyramidValueType>::zero());
    //}
    //collapsePyramid(false, v);

    for (unsigned int i = 0; i < v->size(); i++) {
        char filenameBuf[512];
        snprintf(filenameBuf, 512, "%s%04u.tif", prefix, i);

        // Rescale the pyramid values to fit in UINT16.
        USRGBImage usPyramid((*v)[i]->width(), (*v)[i]->height());
        transformImage(srcImageRange(*((*v)[i])), destImage(usPyramid),
                linearRangeMapping(PyramidVectorType(NumericTraits<PyramidValueType>::min()),
                                   PyramidVectorType(NumericTraits<PyramidValueType>::max()),
                                   typename USRGBImage::value_type(NumericTraits<unsigned short>::min()),
                                   typename USRGBImage::value_type(NumericTraits<unsigned short>::max())));

        ImageExportInfo info(filenameBuf);
        exportImage(srcImageRange(usPyramid), info);
    }
};

// Export a pyramid as a set of UINT16 tiff files.
template <typename PyramidImageType>
void exportPyramid(vector<PyramidImageType*> *v, char *prefix) {
    typedef typename NumericTraits<typename PyramidImageType::value_type>::isScalar pyramid_is_scalar;
    exportPyramid(v, prefix, pyramid_is_scalar());
};


} // namespace enblend

#endif /* __PYRAMID_H__ */
