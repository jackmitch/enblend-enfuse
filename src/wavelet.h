/*
 * Copyright (C) 2005 Andrew Mihal
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
#ifndef __WAVELET_H__
#define __WAVELET_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vigra/basicimage.hxx"

using std::max;

using vigra::BasicImage;

namespace enblend {

// Forward CDF(2,2) wavelet transform with integer lifting
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void wavelet(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;
    typedef BasicImage<DestPixelType> TempImage;
    typedef typename TempImage::traverser TempImageIterator;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in wavelet");

    TempImage tempS((max(src_w, src_h) + 1) >> 1, 1);
    TempImage tempD((max(src_w, src_h) + 1) >> 1, 1);

    // Transform rows from source to dest
    SrcImageIterator sy = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dy = dest_upperleft;
    for (; sy.y != send.y; ++sy.y, ++dy.y) {

        SrcImageIterator sx = sy;

        // Splitting from source into tempS and tempD.
        TempImageIterator tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; sx.x < (send.x - 1); ++sx.x, ++tsi.x, ++tdi.x) {
            *tsi = sa(sx);
            ++sx.x;
            *tdi = sa(sx);
        }
        // Odd sample
        if (sx.x != send.x) {
            *tsi = sa(sx);
        }

        // Dual lifting, from tempS and tempD to second half of destImage
        DestImageIterator dx = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator dend = dest_lowerright;
        if ((src_w & 1) == 0) dend += Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++tdi.x) {
            *dx = *tdi - ((*tsi + tsi(1, 0)) >> 1);
        }
        if (wraparound) {
            *dx = *tdi - ((*tsi + *(tempS.upperLeft())) >> 1);
        } else if ((src_w & 1) == 0){
            // Mirror edge treatment
            *dx = *tdi - *tsi;
        }

        // Primal lifting - from tempS and second half of destImage to first half of destImage
        dx = dy;
        dend = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator dualLiftResult = dend;
        if ((src_w & 1) == 1) dend += Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        if (wraparound) {
            *dx = *tsi + ((*dualLiftResult + dualLiftResult(((src_w + 1) >> 1) - 1, 0)) >> 2);
        } else {
            // Mirror edge treatment
            *dx = *tsi + (*dualLiftResult >> 1);
        }
        ++dx.x;
        ++tsi.x;
        ++dualLiftResult.x;
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++dualLiftResult.x) {
            *dx = *tsi + ((dualLiftResult(-1, 0) + *dualLiftResult) >> 2);
        }
        if ((src_w & 1) == 1) {
            *dx = *tsi + (dualLiftResult(-1, 0) >> 1);
        }

    }

    // Transform columns from dest to dest
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; dx.x != dend.x; ++dx.x) {

        dy = dx;

        // Split from dest into tempS and tempD.
        TempImageIterator tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; dy.y < (dend.y - 1); ++dy.y, ++tsi.x, ++tdi.x) {
            *tsi = da(dy);
            ++dy.y;
            *tdi = da(dy);
        }
        // Odd sample
        if (dy.y != dend.y) {
            *tsi = da(dy);
        }

        // Dual lifting, from tempS and tempD to second half of destImage
        dy = dx + Diff2D(0, (src_h + 1) >> 1);
        DestImageIterator dend = dest_lowerright;
        if ((src_h & 1) == 0) dend += Diff2D(0, -1);
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++tdi.x) {
            *dy = *tdi - ((*tsi + tsi(1, 0)) >> 1);
        }
        if ((src_h & 1) == 0) {
            // Mirror edge treatment
            *dy = *tdi - *tsi;
        }

        // Primal lifting - from tempS and second half of destImage to first half of destImage
        dy = dx;
        dend = dx + Diff2D(0, (src_h + 1) >> 1);
        DestImageIterator dualLiftResult = dend;
        if ((src_h & 1) == 1) dend += Diff2D(0, -1);
        tsi = tempS.upperLeft();
        // Mirror edge treatment
        *dy = *tsi + (*dualLiftResult >> 1);
        ++dy.y;
        ++tsi.x;
        ++dualLiftResult.y;
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++dualLiftResult.y) {
            *dy = *tsi + ((dualLiftResult(0, -1) + *dualLiftResult) >> 2);
        }
        if ((src_h & 1) == 1) {
            *dy = *tsi + (dualLiftResult(0, -1) >> 1);
        }

    }

};

template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void wavelet(int levels, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    if (levels == 0) return;

    wavelet(wraparound,
            src_upperleft, src_lowerright, sa,
            dest_upperleft, dest_lowerright, da);

    for (int i = 1; i < levels; i++) {
        int dest_w = dest_lowerright.x - dest_upperleft.x;
        int dest_h = dest_lowerright.y - dest_upperleft.y;
        dest_lowerright = dest_upperleft + Diff2D((dest_w + 1) >> 1, (dest_h + 1) >> 1);

        wavelet(wraparound,
                dest_upperleft, dest_lowerright, da,
                dest_upperleft, dest_lowerright, da);
    }
};

// Inverse CDF(2,2) wavelet transform with integer lifting
template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void iwavelet(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    typedef typename DestAccessor::value_type DestPixelType;
    typedef BasicImage<DestPixelType> TempImage;
    typedef typename TempImage::traverser TempImageIterator;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in wavelet");

    TempImage tempS((max(src_w, src_h) + 1) >> 1, 1);
    TempImage tempD((max(src_w, src_h) + 1) >> 1, 1);

    // Transform columns from source to dest.
    SrcImageIterator sx = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; sx.x != send.x; ++sx.x, ++dx.x) {

        // Inverse primal lifting from source to tempS
        SrcImageIterator si = sx;
        SrcImageIterator siend = sx + Diff2D(0, (src_h + 1) >> 1);
        SrcImageIterator di = siend;
        if ((src_h & 1) == 1) siend += Diff2D(0, -1);
        TempImageIterator tsi = tempS.upperLeft();
        // Mirror edge treatment
        *tsi = *si - (*di >> 1);
        ++si.y;
        ++di.y;
        ++tsi.x;
        for (; si.y != siend.y; ++si.y, ++di.y, ++tsi.x) {
            *tsi = *si - ((di(0, -1) + *di) >> 2);
        }
        if ((src_h & 1) == 1) {
            *tsi = *si - (di(0, -1) >> 1);
        }

        // Inverse dual lifting from source and tempS to tempD
        di = sx + Diff2D(0, (src_h + 1) >> 1);
        SrcImageIterator diend = src_lowerright;
        if ((src_h & 1) == 0) diend += Diff2D(0, -1);
        tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; di.y != diend.y; ++di.y, ++tsi.x, ++tdi.x) {
            *tdi = *di + ((*tsi + tsi(1, 0)) >> 1);
        }
        if ((src_h & 1) == 0) {
            // Mirror edge treatment
            *tdi = *di + *tsi;
        }

        // Merging from tempS and tempD to destImage
        DestImageIterator dy = dx;
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dy.y < (dend.y - 1); ++dy.y, ++tsi.x, ++tdi.x) {
            *dy = *tsi;
            ++dy.y;
            *dy = *tdi;
        }
        // Odd sample
        if (dy.y != dend.y) {
            *dy = *tsi;
        }

    }

    // Transform rows from dest to dest.
    DestImageIterator dy = dest_upperleft;
    for (; dy.y != dend.y; ++dy.y) {

        // Inverse primal lifting from dest to tempS
        DestImageIterator si = dy;
        DestImageIterator siend = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator di = siend;
        if ((src_w & 1) == 1) siend += Diff2D(-1, 0);
        TempImageIterator tsi = tempS.upperLeft();
        if (wraparound) {
            *tsi = *si - ((di(((src_w + 1) >> 1) - 1, 0) + *di) >> 2);
        } else {
            // Mirror edge treatment
            *tsi = *si - (*di >> 1);
        }
        ++si.x;
        ++di.x;
        ++tsi.x;
        for (; si.x != siend.x; ++si.x, ++di.x, ++tsi.x) {
            *tsi = *si - ((di(-1, 0) + *di) >> 2);
        }
        if ((src_w & 1) == 1) {
            *tsi = *si - (di(-1, 0) >> 1);
        }

        // Inverse dual lifting from dest and tempS to tempD
        di = dy + Diff2D((src_w + 1) >> 1, 0);
        DestImageIterator diend = dend;
        if ((src_w & 1) == 0) diend += Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; di.x != diend.x; ++di.x, ++tsi.x, ++tdi.x) {
            *tdi = *di + ((*tsi + tsi(1, 0)) >> 1);
        }
        if (wraparound) {
            *tdi = *di + ((*tsi + *(tempS.upperLeft())) >> 1);
        } else if ((src_w & 1) == 0) {
            // Mirror edge treatment
            *tdi = *di + *tsi;
        }

        // Merging from tempS and tempD to destImage
        dx = dy;
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dx.x < (dend.x - 1); ++dx.x, ++tsi.x, ++tdi.x) {
            *dx = *tsi;
            ++dx.x;
            *dx = *tdi;
        }
        // Odd sample
        if (dx.x != dend.x) {
            *dx = *tsi;
        }
    }

};

template <typename SrcImageIterator, typename SrcAccessor,
        typename DestImageIterator, typename DestAccessor>
inline void iwavelet(int levels, bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestImageIterator dest_lowerright,
        DestAccessor da) {

    copyImage(src_upperleft, src_lowerright, sa, dest_upperleft, da);

    if (levels == 0) return;

    DestImageIterator *lowerRightArray = new DestImageIterator[levels];
    lowerRightArray[0] = dest_lowerright;
    for (int i = 1; i < levels; i++) {
        int dest_w = dest_lowerright.x - dest_upperleft.x;
        int dest_h = dest_lowerright.y - dest_upperleft.y;
        dest_lowerright = dest_upperleft + Diff2D((dest_w + 1) >> 1, (dest_h + 1) >> 1);
        lowerRightArray[i] = dest_lowerright;
    }

    for (int i = levels - 1; i >= 0; --i) {
        iwavelet(wraparound,
                dest_upperleft, lowerRightArray[i], da,
                dest_upperleft, lowerRightArray[i], da);
    }

    delete [] lowerRightArray;

}

template <typename SrcImageIterator, typename SrcAccessor,
        typename MaskIterator, typename MaskAccessor,
        typename DestImageIterator, typename DestAccessor,
        typename DestMaskIterator, typename DestMaskAccessor>
inline void wavelet(int levels, bool wraparound,
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

};

} // namespace enblend

#endif /* __WAVELET_H__ */
