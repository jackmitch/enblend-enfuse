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

using std::cout;
using std::cerr;
using std::endl;
using std::pair;

using vigra::NumericTraits;
using vigra::triple;
using vigra::UIImage;

namespace enblend {

#define EUCLIDEAN_METRIC

template <typename dist_t>
inline dist_t _nftDistance(dist_t deltaX, dist_t deltaY) {
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
inline dist_t _nftDistance(dist_t deltaY) {
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
void nearestFeatureTransform(bool wraparound,
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestAccessor da) {

    typedef UIImage::Iterator DnfIterator;
    typedef typename SrcAccessor::value_type SrcValueType;

    SrcImageIterator sx, sy, send, smidpoint;
    DnfIterator dnfcx, dnfcy;
    DnfIterator dnflx, dnfly;
    DestImageIterator dx, dy;

    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    UIImage dnfColumn(w, h);
    UIImage dnfLeft(w, h);

    // Initialize dnfColumn top-down. Store the distance to the nearest feature
    // in the same column and above us.
    if (Verbose > 0) {
        if (wraparound) cout << "Creating blend mask: 1/6";
        else cout << "Creating blend mask: 1/4";
        cout.flush();
    }
    sx = src_upperleft;
    send = src_lowerright;
    dnfcx = dnfColumn.upperLeft();
    dx = dest_upperleft;
    for (; sx.x != send.x; ++sx.x, ++dnfcx.x, ++dx.x) {
        sy = sx;
        dnfcy = dnfcx;
        dy = dx;

        // Color of the last feature pixel.
        SrcValueType lastFeature = sa(src_upperleft);
        bool foundFirstFeature = false;
        unsigned int lastFeatureDeltaY = 0;

        for (; sy.y != send.y; ++sy.y, ++dnfcy.y, ++dy.y, ++lastFeatureDeltaY) {
            if (sa(sy)) {
                // Source pixel is a feature pixel.
                lastFeature = sa(sy);
                foundFirstFeature = true;
                // Distance to feature pixel = 0
                *dnfcy = 0;
                lastFeatureDeltaY = 0;
                // Nearest feature color = source feature color.
                da.set(lastFeature, dy);
            }
            else if (foundFirstFeature) {
                // Source pixel is not a feature.
                *dnfcy = _nftDistance(lastFeatureDeltaY);
                da.set(lastFeature, dy);
            }
            else {
                *dnfcy = UINT_MAX;
            }
        }
    }

    // Initialize dnfColumn bottom-up. Caluclate the distance to the nearest
    // feature in the same column and below us.
    // If this is smaller than the value caluclated in the top-down pass,
    // overwrite that value.
    if (Verbose > 0) {
        if (wraparound) cout << " 2/6";
        else cout << " 2/4";
        cout.flush();
    }
    sx = src_lowerright;
    send = src_upperleft;
    dnfcx = dnfColumn.lowerRight();
    dx = dest_upperleft + (src_lowerright - src_upperleft);
    for (; sx.x != send.x;) {
        --sx.x;
        --dnfcx.x;
        --dx.x;

        sy = sx;
        dnfcy = dnfcx;
        dy = dx;

        // Color of the last feature pixel.
        SrcValueType lastFeature = sa(src_upperleft);
        bool foundFirstFeature = false;
        unsigned int lastFeatureDeltaY = 0;

        for (; sy.y != send.y; ++lastFeatureDeltaY) {
            --sy.y;
            --dnfcy.y;
            --dy.y;

            if (sa(sy)) {
                // Source pixel is a feature pixel.
                lastFeature = sa(sy);
                foundFirstFeature = true;
                // Distance to feature pixel = 0
                *dnfcy = 0;
                lastFeatureDeltaY = 0;
                // Nearest feature color = source feature color.
                da.set(lastFeature, dy);
            }
            else if (foundFirstFeature) {
                // Source pixel is not a feature
                unsigned int distLastFeature = _nftDistance(lastFeatureDeltaY);
                if (distLastFeature < *dnfcy) {
                    // Feature below us is closer than feature above us.
                    *dnfcy = distLastFeature;
                    da.set(lastFeature, dy);
                }
            }
        }
    }

    // Calculate dnfLeft for each pixel.
    if (Verbose > 0) {
        if (wraparound) cout << " 3/6";
        else cout << " 3/4";
        cout.flush();
    }
    sy = src_upperleft;
    send = src_lowerright;
    smidpoint = src_upperleft + Diff2D(w/2, h/2);
    dnfcy = dnfColumn.upperLeft();
    dnfly = dnfLeft.upperLeft();
    dy = dest_upperleft;
    for (; sy.y != send.y; ++sy.y, ++dnfcy.y, ++dnfly.y, ++dy.y) {

        // Indicate halfway mark when wraparound is true.
        if (Verbose > 0 && wraparound && (sy.y == smidpoint.y)) {
            cout << " 4/6";
            cout.flush();
        }

        // List of dnfcx's on the left that might be the closest features
        // to the current dnflx.
        list<DnfIterator> potentialFeatureList;

        for (int twiceAround = (wraparound?1:0); twiceAround >= 0; twiceAround--) {
            sx = sy;
            dnfcx = dnfcy;
            dnflx = dnfly;
            dx = dy;

            for (; sx.x != send.x; ++sx.x, ++dnfcx.x, ++dnflx.x, ++dx.x) {
                // First add ourself to the list.
                potentialFeatureList.push_back(dnfcx);

                // Iterate throught the list starting at the right. For each
                // potential feature, all of the potential features to the left
                // in the list must be strictly closer. If not delete them from
                // the list.
                list<DnfIterator>::iterator potentialFeature =
                        --(potentialFeatureList.end());
                // The last potential feature is dnfcx, just added above.
                // That is in the current column so the distance to that feature
                // is simply *dnfcx.
                unsigned int distPotentialFeature = *dnfcx;
                while (potentialFeature != potentialFeatureList.begin()) {
                    // Make an iterator that points to the predecessor.
                    list<DnfIterator>::iterator previousFeature = potentialFeature;
                    --previousFeature;

                    // Subtract the iterators .x components to find out how many
                    // columns to the left of dnfcx previousFeature is.
                    // DeltaX must be positive.
                    // modulo w to consider wraparound condition.
                    int deltaX = (dnfcx.x - (*previousFeature).x) % w;
                    if (deltaX < 0) deltaX += w;

                    // previousFeature is this far from dnfcx.
                    unsigned int distPreviousFeature =
                            _nftDistance((unsigned int)deltaX, **previousFeature);

                    if (distPreviousFeature >= distPotentialFeature) {
                        // previousFeature is not a candidate for dnflx
                        // or any dnflx further to the right.
                        potentialFeatureList.erase(previousFeature);
                    } else {
                        // previousFeature is a candidate.
                        potentialFeature = previousFeature;
                        distPotentialFeature = distPreviousFeature;
                    }
                }

                // The closest feature to dnflx in columns <= dnflx is the first
                // potential feature in the list.
                *dnflx = distPotentialFeature;

                // Set color of dx to be color of closest feature to the left.
                da.set(dx(((*potentialFeature).x - dnfcx.x), 0), dx);
            }
        }
    }

    // Final pass: calculate the distance to the nearest feature in the same
    // column or any column to the right. If this is smaller than dnflx,
    // Then recolor the pixel to the color of the nearest feature to the right.
    if (Verbose > 0) {
        if (wraparound) cout << " 5/6";
        else cout << " 4/4";
        cout.flush();
    }
    sy = src_lowerright;
    send = src_upperleft;
    smidpoint = src_upperleft + Diff2D(w/2, h/2);
    dnfcy = dnfColumn.lowerRight();
    dnfly = dnfLeft.lowerRight();
    dy = dest_upperleft + (src_lowerright - src_upperleft);
    for (; sy.y != send.y;) {
        --sy.y;
        --dnfcy.y;
        --dnfly.y;
        --dy.y;

        // Indicate halfway mark when wraparound is true.
        if (Verbose > 0 && wraparound && (sy.y == smidpoint.y)) {
            cout << " 6/6";
            cout.flush();
        }

        // List of dnfcx's on the right that might be the closest features
        // to the current dnflx.
        list<DnfIterator> potentialFeatureList;

        for (int twiceAround = (wraparound?1:0); twiceAround >= 0; twiceAround--) {
            sx = sy;
            dnfcx = dnfcy;
            dnflx = dnfly;
            dx = dy;

            for (; sx.x != send.x;) {
                --sx.x;
                --dnfcx.x;
                --dnflx.x;
                --dx.x;

                // First add ourself to the list.
                potentialFeatureList.push_back(dnfcx);

                // Iterate through list and prune as before.
                list<DnfIterator>::iterator potentialFeature =
                        --(potentialFeatureList.end());
                unsigned int distPotentialFeature = *dnfcx;
                while (potentialFeature != potentialFeatureList.begin()) {
                    // Iterator that points to predecessor.
                    list<DnfIterator>::iterator previousFeature = potentialFeature;
                    --previousFeature;

                    // Subtract the iterators .x components to find out how many
                    // columns to the right of dnfcx previousFeature is.
                    // DeltaX must be positive.
                    // modulo w to consider wraparound condition.
                    int deltaX = ((*previousFeature).x - dnfcx.x) % w;
                    if (deltaX < 0) deltaX += w;

                    // previousFeature is this far from dnfcx.
                    unsigned int distPreviousFeature =
                            _nftDistance((unsigned int)deltaX, **previousFeature);

                    if (distPreviousFeature >= distPotentialFeature) {
                        // previousFeature is not a candidate.
                        potentialFeatureList.erase(previousFeature);
                    } else {
                        // previousFeature is a candidate.
                        potentialFeature = previousFeature;
                        distPotentialFeature = distPreviousFeature;
                    }
                }

                // The closest feature on the right is potentialFeature.
                if (*dnflx > distPotentialFeature) {
                    // Recolor dx.
                    da.set(dx(((*potentialFeature).x - dx.x), 0), dx);
                }
            }
        }
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return;
};

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void nearestFeatureTransform(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<DestImageIterator, DestAccessor> dest) {

    nearestFeatureTransform(wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second);

};

} // namespace enblend

#endif /* __NEAREST_H__ */
