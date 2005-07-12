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

    typedef typename DestAccessor::value_type DestPixelType;
    typedef BasicImage<DestPixelType> TempImage;
    typedef typename TempImage::traverser TempImageIterator;

    int src_w = src_lowerright.x - src_upperleft.x;
    int src_h = src_lowerright.y - src_upperleft.y;
    vigra_precondition(src_w > 1 && src_h > 1,
            "src image too small in reduce");
    vigra_precondition((src_w & 1) == 0 && (src_h & 1) == 0,
            "src image must have even dimensions");

    TempImage si(max(src_w, src_h)/2, 1);
    TempImage di(max(src_w, src_h)/2, 1);

    SrcImageIterator sy = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DestImageIterator dy = dest_upperleft;

    for (; sy.y != send.y; ++sy.y, ++dy.y) {

        SrcImageIterator sx = sy;

        // Splitting into si and di.
        TempImageIterator six = si.upperLeft();
        TempImageIterator dix = di.upperLeft();
        for (; sx.x != send.x; ++sx.x, ++six.x, ++dix.x) {
            *six = sa(sx);
            ++sx.x;
            *dix = sa(sx);
        }

        // Dual lifting.
        DestImageIterator dx = dy + Diff2D(src_w/2, 0);
        DestImageIterator dend = dest_lowerright + Diff2D(-1, 0);
        six = si.upperLeft();
        dix = di.upperLeft();
        for (; dx.x != dend.x; ++dx.x, ++six.x, ++dix.x) {
            // FIXME use shift instead of double divide
            *dx = *dix - ((*six + six(1, 0)) / 2);
        }
        // Mirror edge treatment?
        *dx = *dix - *six;

        // Primal lifting
        dx = dy;
        dend = dy + Diff2D(src_w/2, 0);
        DestImageIterator dualLiftResult = dend;
        six = si.upperLeft();
        // Mirror edge treatment?
        // FIXME use shift instead of double divide
        *dx = *six + (*dualLiftResult / 2);
        ++dx.x;
        ++six.x;
        ++dualLiftResult.x;
        for (; dx.x != dend.x; ++dx.x, ++six.x, ++dualLiftResult.x) {
            // FIXME use shift instead of double divide
            *dx = *six + ((dualLiftResult(-1, 0) + *dualLiftResult) / 4);
        }

    }

    SrcImageIterator sx = src_upperleft;
    SrcImageIterator dx = dest_upperleft;

    for (; sx.x != send.x; ++sx.x, ++dx.x) {

        sy = sx;

        // Split into si and di.
        TempImageIterator six = si.upperLeft();
        TempImageIterator dix = di.upperLeft();
        for (; sy.y != send.y; ++sy.y, ++six.x, ++dix.x) {
            *six = sa(sy);
            ++sy.y;
            *dix = sa(sy);
        }

        // Dual lifting.
        dy = dx + Diff2D(0, src_h/2);
        DestImageIterator dend = dest_lowerright + Diff2D(0, -1);
        six = si.upperLeft();
        dix = di.upperLeft();
        for (; dy.y != dend.y; ++dy.y, ++six.x, ++dix.x) {
            // FIXME use shift instead of double divide
            *dy = *dix - ((*six + six(1, 0)) / 2);
        }
        // Mirror edge treatment?
        *dy = *dix - *six;

        // Primal lifting
        dy = dx;
        dend = dx + Diff2D(0, src_h/2);
        DestImageIterator dualLiftResult = dend;
        six = si.upperLeft();
        // Mirror edge treatment?
        // FIXME use shift instead of double divide
        *dy = *six + (*dualLiftResult / 2);
        ++dy.y;
        ++six.x;
        ++dualLiftResult.y;
        for (; dy.y != dend.y; ++dy.y, ++six.x, ++dualLiftResult.y) {
            // FIXME use shift instead of double divide
            *dy = *six + ((dualLiftResult(0, -1) + *dualLiftResult) / 4);
        }

    }

};

} // namespace enblend

#endif /* __WAVELET_H__ */
