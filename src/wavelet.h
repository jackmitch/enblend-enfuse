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
    vigra_precondition((src_w & 1) == 0 && (src_h & 1) == 0,
            "src image must have even dimensions");

    TempImage tempS(max(src_w, src_h)/2, 1);
    TempImage tempD(max(src_w, src_h)/2, 1);

    // Transform rows from source to dest
    SrcImageIterator sy = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dy = dest_upperleft;
    for (; sy.y != send.y; ++sy.y, ++dy.y) {

        SrcImageIterator sx = sy;

        // Splitting from source into tempS and tempD.
        TempImageIterator tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; sx.x != send.x; ++sx.x, ++tsi.x, ++tdi.x) {
            *tsi = sa(sx);
            ++sx.x;
            *tdi = sa(sx);
        }

        // Dual lifting, from tempS and tempD to second half of destImage
        DestImageIterator dx = dy + Diff2D(src_w/2, 0);
        DestImageIterator dend = dest_lowerright + Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++tdi.x) {
            *dx = *tdi - ((*tsi + tsi(1, 0)) >> 1);
        }
        // Mirror edge treatment
        //*dx = (*tdi - *tsi) >> 1;
        *dx = *tdi - *tsi;

        // Primal lifting - from tempS and second half of destImage to first half of destImage
        dx = dy;
        dend = dy + Diff2D(src_w/2, 0);
        DestImageIterator dualLiftResult = dend;
        tsi = tempS.upperLeft();
        // Mirror edge treatment
        //*dx = *tsi + ((*tsi + *dualLiftResult) >> 3);
        *dx = *tsi + (*dualLiftResult >> 1);
        ++dx.x;
        ++tsi.x;
        ++dualLiftResult.x;
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++dualLiftResult.x) {
            *dx = *tsi + ((dualLiftResult(-1, 0) + *dualLiftResult) >> 2);
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
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++tdi.x) {
            *tsi = da(dy);
            ++dy.y;
            *tdi = da(dy);
        }

        // Dual lifting, from tempS and tempD to second half of destImage
        dy = dx + Diff2D(0, src_h/2);
        DestImageIterator dend = dest_lowerright + Diff2D(0, -1);
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++tdi.x) {
            *dy = *tdi - ((*tsi + tsi(1, 0)) >> 1);
        }
        // Mirror edge treatment?
        //*dy = (*tdi - *tsi) >> 1;
        *dy = *tdi - *tsi;

        // Primal lifting - from tempS and second half of destImage to first half of destImage
        dy = dx;
        dend = dx + Diff2D(0, src_h/2);
        DestImageIterator dualLiftResult = dend;
        tsi = tempS.upperLeft();
        // Mirror edge treatment?
        //*dy = *tsi + ((*tsi + *dualLiftResult) >> 3);
        *dy = *tsi + (*dualLiftResult >> 1);
        ++dy.y;
        ++tsi.x;
        ++dualLiftResult.y;
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++dualLiftResult.y) {
            *dy = *tsi + ((dualLiftResult(0, -1) + *dualLiftResult) >> 2);
        }

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
    vigra_precondition((src_w & 1) == 0 && (src_h & 1) == 0,
            "src image must have even dimensions");

    TempImage tempS(max(src_w, src_h)/2, 1);
    TempImage tempD(max(src_w, src_h)/2, 1);

    // Tranform columns from source to dest.
    SrcImageIterator sx = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dx = dest_upperleft;
    DestImageIterator dend = dest_lowerright;
    for (; sx.x != send.x; ++sx.x, ++dx.x) {

        // Inverse primal lifting from source to tempS
        SrcImageIterator si = sx;
        SrcImageIterator siend = sx + Diff2D(0, src_h/2);
        SrcImageIterator di = siend;
        TempImageIterator tsi = tempS.upperLeft();
        // Mirror edge treatment
        //*tsi = *si - ((*si + *di) >> 2);
        *tsi = *si - (*di >> 1);
        ++si.y;
        ++di.y;
        ++tsi.x;
        for (; si.y != siend.y; ++si.y, ++di.y, ++tsi.x) {
            *tsi = *si - ((di(0, -1) + *di) >> 2);
        }

        // Inverse dual lifting from source and tempS to tempD
        di = siend;
        SrcImageIterator diend = send + Diff2D(0, -1);
        tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; di.y != diend.y; ++di.y, ++tsi.x, ++tdi.x) {
            *tdi = *di + ((*tsi + tsi(1, 0)) >> 1);
        }
        // Mirror edge treatment
        //*tdi = *di + ((*tsi + *di) >> 2);
        *tdi = *di + *tsi;

        // Merging from tempS and tempD to destImage
        DestImageIterator dy = dx;
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dy.y != dend.y; ++dy.y, ++tsi.x, ++tdi.x) {
            *dy = *tsi;
            ++dy.y;
            *dy = *tdi;
        }

    }

    // Transform rows from dest to dest.
    DestImageIterator dy = dest_upperleft;
    for (; dy.y != dend.y; ++dy.y) {

        // Inverse primal lifting from dest to tempS
        DestImageIterator si = dy;
        DestImageIterator siend = dy + Diff2D(src_w/2, 0);
        DestImageIterator di = siend;
        TempImageIterator tsi = tempS.upperLeft();
        // Mirror edge treatment
        //*tsi = *si - ((*si + *di) >> 2);
        *tsi = *si - (*di >> 1);
        ++si.x;
        ++di.x;
        ++tsi.x;
        for (; si.x != siend.x; ++si.x, ++di.x, ++tsi.x) {
            *tsi = *si - ((di(-1, 0) + *di) >> 2);
        }

        // Inverse dual lifting from dest and tempS to tempD
        di = siend;
        DestImageIterator diend = dend + Diff2D(-1, 0);
        tsi = tempS.upperLeft();
        TempImageIterator tdi = tempD.upperLeft();
        for (; di.x != diend.x; ++di.x, ++tsi.x, ++tdi.x) {
            *tdi = *di + ((*tsi + tsi(1, 0)) >> 1);
        }
        // Mirror edge treatment
        //*tdi = *di + ((*tsi + *di) >> 2);
        *tdi = *di + *tsi;

        // Merging from tempS and tempD to destImage
        dx = dy;
        tsi = tempS.upperLeft();
        tdi = tempD.upperLeft();
        for (; dx.x != dend.x; ++dx.x, ++tsi.x, ++tdi.x) {
            *dx = *tsi;
            ++dx.x;
            *dx = *tdi;
        }
    }

};

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
