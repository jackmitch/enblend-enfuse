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
#ifndef __NEAREST_H__
#define __NEAREST_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <utility>

#include "vigra/numerictraits.hxx"
#include "vigra/stdimage.hxx"

using namespace std;

using vigra::NumericTraits;
using vigra::triple;
using vigra::UIImage;

#define EUCLIDEAN_METRIC

template <typename dist_t>
inline dist_t distance(dist_t deltaX, dist_t deltaY) {
    #ifdef EUCLIDEAN_METRIC
        return (deltaY == NumericTraits<dist_t>::max())
                ? deltaY
                : (deltaY + (deltaX * deltaX));
    #else
    #ifdef CHESSBOARD_METRIC
        return max(deltaX, deltaY);
    #else
    #ifdef MANHATTAN_METRIC
        return (deltaY == NumericTraits<dist_t>::max())
                ? deltaY
                : (deltaX + deltaY);
    #endif
    #endif
    #endif
};

// Distance to a pixel with the same x coordinate.
template <typename dist_t>
inline dist_t distance(dist_t deltaY) {
    #ifdef EUCLIDEAN_METRIC
        return deltaY * deltaY;
    #else
    #ifdef CHESSBOARD_METRIC
        return deltaY;
    #else
    #ifdef MANHATTAN_METRIC
        return deltaY;
    #endif
    #endif
    #endif
};

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void nearestFeatureTransform(SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestAccessor da) {

    typedef UIImage::Iterator DnfIterator;
    typedef typename SrcAccessor::value_type SrcValueType;

    int w = (src_lowerright.x - src_upperleft.x) + 1;
    int h = (src_lowerright.y - src_upperleft.y) + 1;

    UIImage dnfColumn(w, h);
    UIImage dnfLeft(w, h);

    // Initialize dnfColumn top-down. Store the distance to the nearest feature
    // in the same column and above us.
    if (Verbose > 0) {
        cout << "Creating blend mask: 1/4 ";
        cout.flush();
    }
    SrcImageIterator sx = src_upperleft;
    SrcImageIterator send = src_lowerright;
    DnfIterator dnfx = dnfColumn.upperLeft();
    DestImageIterator dx = dest_upperleft;
    for (; sx.x != send.x; ++sx.x, ++dnfx.x, ++dx.x) {
        SrcImageIterator sy = sx;
        DnfIterator dnfy = dnfx;
        DestImageIterator dy = dx;

        // Color of the last feature pixel.
        SrcValueType lastFeature;
        bool foundFirstFeature = false;
        unsigned int lastFeatureDeltaY = 0;

        for (; sy.y != send.y; ++sy.y, ++dnfy.y, ++dy.y, ++lastFeatureDeltaY) {
            if (sa(sy)) {
                // Source pixel is a feature pixel.
                lastFeature = sa(sy);
                foundFirstFeature = true;
                // Distance to feature pixel = 0
                *dnfy = 0;
                // Nearest feature color = source feature color.
                da.set(lastFeature, dy);
            } else if (foundFirstFeature) {
                // Source pixel is not a feature.
                *dnfy = distance(lastFeatureDeltaY);
                da.set(lastFeature, dy);
            } else {
                *dnfy = UINT_MAX;
            }
        }
    }

    // Initialize dnfColumn bottom-up. Caluclate the distance to the nearest
    // feature in the same column and below us.
    // If this is smaller than the value caluclated in the top-down pass,
    // overwrite that value.
    if (Verbose > 0) {
        cout << "2/4 ";
        cout.flush();
    }
    sx = src_upperleft;
    send = src_lowerright;
    dnfx = dnfColumn.upperLeft();
    dx = dest_upperleft;
    for (; sx.x != send.x; ++sx.x, ++dnfx.x, ++dx.x) {
        SrcImageIterator sy = sx;
        DnfIterator dnfy = dnfx;
        DestImageIterator dy = dx;

        // Color of the last feature pixel.
        SrcValueType lastFeature;
        bool foundFirstFeature = false;
        unsigned int lastFeatureDeltaY = 0;
    }
};

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void nearestFeatureTransform(
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<DestImageIterator, DestAccessor> dest) {

    nearestFeatureTransform(src.first, src.second, src.third,
            dest.first, dest.second);

};

#endif /* __NEAREST_H__ */
