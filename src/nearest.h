/*
 * Copyright (C) 2004-2009 Andrew Mihal
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

#include "vigra/numerictraits.hxx"
#include "vigra/stdcachedfileimage.hxx"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;

using vigra::NumericTraits;
using vigra::triple;

namespace enblend {

// The metric to use for calculating distances.
#define EUCLIDEAN_METRIC
//#define CHESSBOARD_METRIC
//#define MANHATTAN_METRIC
//#define HYBRID_METRIC


template <typename dist_t>
inline dist_t _nftDistance(dist_t deltaX, dist_t deltaY)
{
#if defined(EUCLIDEAN_METRIC)
    const dist_t max = NumericTraits<dist_t>::max();
    const double z = rint(hypot(static_cast<double>(deltaX), static_cast<double>(deltaY)));
    return z > static_cast<double>(max) ? max : static_cast<dist_t>(z);
#elif defined(CHESSBOARD_METRIC)
    return std::max(deltaX, deltaY);
#elif defined(MANHATTAN_METRIC)
    const dist_t max = NumericTraits<dist_t>::max();
    const double sum = static_cast<double>(deltaX) + static_cast<double>(deltaY);
    return sum > static_cast<double>(max) ? max : static_cast<dist_t>(sum);
#elif defined(HYBRID_METRIC)
    const dist_t max = NumericTraits<dist_t>::max();
    const double sum =
        static_cast<double>(std::max(deltaX, deltaY)) +
        0.5 * static_cast<double>(std::min(deltaX, deltaY));
    return sum > static_cast<double>(max) ? max : static_cast<dist_t>(sum);
#else
#error no metric defined
#endif
}


// Distance to a pixel with the same x coordinate.
template <typename dist_t>
inline dist_t _nftDistance(dist_t deltaY)
{
#if defined(EUCLIDEAN_METRIC)
    return deltaY;
#elif defined(CHESSBOARD_METRIC)
    return deltaY;
#elif defined(MANHATTAN_METRIC)
    return deltaY;
#elif defined(HYBRID_METRIC)
    return deltaY;
#else
#error no metric defined
#endif
}


/** Data structure for potentialFeatureStack.
 *  Contribution from Fulvio Senore.
 *  Fast insert and delete, avoiding dynamic memory allocation.
 */
template <typename T>
class FeatureList
{
    struct Node {
        T value;
        Node *prev;
    };

    // Statically allocated memory for Nodes
    Node *array;

    // pointer to the first unused element of the array
    Node *firstUnused;

    // pointer to the end of the list
    Node *last;

    Node *iterator;

public:
    FeatureList(int size) :
        array(new Node[size]), firstUnused(array), last(NULL), iterator(NULL)
    {}

    ~FeatureList() {
        delete [] array;
    }

    void clear() {
        firstUnused = array;
        last = NULL;
        iterator = NULL;
    }

    void push_back(const T& value) {
        firstUnused->value = value;
        firstUnused->prev = last;
        last = firstUnused;
        firstUnused++;
    }

    void move_to_beginning() {
        iterator = array;
    }

    void move_to_end() {
        iterator = last;
    }

    bool has_previous() const {
        return iterator->prev != NULL;
    }

    void erase_previous() {
        iterator->prev = iterator->prev->prev;
    }

    T get_current() const {
        return iterator->value;
    }

    T get_previous() const {
        return iterator->prev->value;
    }

    void move_backwards() {
        iterator = iterator->prev;
    }

    bool empty() const {
        return last == NULL;
    }
};


/** Compute the nearest feature transform.
 *  A non-zero pixel in the src image is considered a feature.  Each
 *  pixel in the dest image is given the value of the nearest feature
 *  to that pixel.
 */
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void nearestFeatureTransform(bool wraparound,
                             SrcImageIterator src_upperleft,
                             SrcImageIterator src_lowerright,
                             SrcAccessor sa,
                             DestImageIterator dest_upperleft,
                             DestAccessor da,
                             typename SrcAccessor::value_type featurelessPixel)
{
    typedef typename EnblendNumericTraits<UInt32>::ImageType DNFImage;
    typedef typename DNFImage::traverser DnfIterator;
    typedef typename SrcAccessor::value_type SrcValueType;

    SrcImageIterator sx, sy, send, smidpoint;
    DnfIterator dnfcx, dnfcy;
    DnfIterator dnflx, dnfly;
    DestImageIterator dx, dy;

    const int w = src_lowerright.x - src_upperleft.x;
    const int h = src_lowerright.y - src_upperleft.y;

    // Distance to the nearest feature in the current column.
    DNFImage* dnfColumn = new DNFImage(w, h);
    // Distance to the nearest feature in the current column, or any
    // column to the left of this column.
    DNFImage* dnfLeft = new DNFImage(w, h);

    // Data structures for initializing dnfColumn.
    // These let us initialize all of the columns in one pass
    // over the rows of the image.  Cache-friendly.
    SrcValueType* lastFeature = new SrcValueType[w];
    bool* foundFirstFeature = new bool[w];
    UInt32* lastFeatureDeltaY = new UInt32[w];

    // Initialize dnfColumn top-down.  Store the distance to the
    // nearest feature in the same column and above us.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) {
            cerr << command << ": info: creating blend mask: 1/6";
        } else {
            cerr << command << ": info: creating blend mask: 1/4";
        }
        cerr.flush();
    }

    // Initialization.
    for (int i = 0; i < w; i++) {
        // before commenting out: 18.23
        // after commenting out:  18.35
        //lastFeature[i] = sa(src_upperleft);
        foundFirstFeature[i] = false;
        //lastFeatureDeltaY[i] = 0;
    }
    sy = src_upperleft;
    send = src_lowerright;
    dnfcy = dnfColumn->upperLeft();
    dy = dest_upperleft;
    for (; sy.y < send.y; ++sy.y, ++dnfcy.y, ++dy.y) {
        sx = sy;
        dnfcx = dnfcy;
        dx = dy;

        for (int xIndex = 0; sx.x < send.x; ++sx.x, ++dnfcx.x, ++dx.x, ++xIndex) {
            if (sa(sx) != featurelessPixel) {
                // Source pixel is a feature pixel.
                lastFeature[xIndex] = sa(sx);
                foundFirstFeature[xIndex] = true;
                // Distance to feature pixel = 0
                *dnfcx = 0;
                lastFeatureDeltaY[xIndex] = 0;
                // Nearest feature color = source feature color.
                da.set(lastFeature[xIndex], dx);
            } else if (foundFirstFeature[xIndex]) {
                // Source pixel is not a feature.
                *dnfcx = _nftDistance(lastFeatureDeltaY[xIndex]);
                da.set(lastFeature[xIndex], dx);
            } else {
                *dnfcx = NumericTraits<UInt32>::max();
            }

            ++lastFeatureDeltaY[xIndex];
        }
    }

    // Initialize dnfColumn bottom-up. Calculate the distance to the
    // nearest feature in the same column and below us.  If this is
    // smaller than the value calculated in the top-down pass,
    // overwrite that value.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) {
            cerr << " 2/6";
        } else {
            cerr << " 2/4";
        }
        cerr.flush();
    }

    // Initialization.
    for (int i = 0; i < w; i++) {
        // before commenting out: 18.35
        // after commenting out:  18.69
        //lastFeature[i] = sa(src_upperleft);
        foundFirstFeature[i] = false;
        // before commenting out this and previous lastFeatureDeltaY init: 18.69
        // after:
        //lastFeatureDeltaY[i] = 0;
    }
    sy = src_lowerright;
    send = src_upperleft;
    dnfcy = dnfColumn->lowerRight();
    dy = dest_upperleft + Diff2D(w, h);
    for (; sy.y > send.y;) {
        --sy.y;
        --dnfcy.y;
        --dy.y;

        sx = sy;
        dnfcx = dnfcy;
        dx = dy;

        for (int xIndex = w - 1; sx.x > send.x; --xIndex) {
            --sx.x;
            --dnfcx.x;
            --dx.x;

            if (sa(sx) != featurelessPixel) {
                // Source pixel is a feature pixel.
                lastFeature[xIndex] = sa(sx);
                foundFirstFeature[xIndex] = true;
                // Distance to feature pixel = 0
                *dnfcx = 0;
                lastFeatureDeltaY[xIndex] = 0;
                // Nearest feature color = source feature color.
                // time before commenting out this line: 19.52
                // time after commenting out this line: 18.23
                //da.set(lastFeature[xIndex], dx);
            } else if (foundFirstFeature[xIndex]) {
                // Source pixel is not a feature
                const UInt32 distLastFeature = _nftDistance(lastFeatureDeltaY[xIndex]);
                if (distLastFeature < *dnfcx) {
                    // Feature below us is closer than feature above us.
                    *dnfcx = distLastFeature;
                    da.set(lastFeature[xIndex], dx);
                }
            }

            ++lastFeatureDeltaY[xIndex];
        }
    }

    // List of dnfcx's on the left that might be the closest features
    // to the current dnflx.
    FeatureList<typename DnfIterator::MoveX> potentialFeatureList(wraparound ? 2 * w : w);

    // Calculate dnfLeft for each pixel.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) {
            cerr << " 3/6";
        } else {
            cerr << " 3/4";
        }
        cerr.flush();
    }
    sy = src_upperleft;
    send = src_lowerright;
    smidpoint = src_upperleft + Diff2D(w / 2, h / 2);
    dnfcy = dnfColumn->upperLeft();
    dnfly = dnfLeft->upperLeft();
    dy = dest_upperleft;
    for (; sy.y < send.y; ++sy.y, ++dnfcy.y, ++dnfly.y, ++dy.y) {

        // Indicate halfway mark when wraparound is true.
        if (Verbose > VERBOSE_NFT_MESSAGES && wraparound && sy.y == smidpoint.y) {
            cerr << " 4/6";
            cerr.flush();
        }

        potentialFeatureList.clear();

        // If wraparound is true, we must go across the row twice.
        // This takes care of the case when the nearest feature is reached by
        // wrapping around the image.
        for (int twiceAround = wraparound ? 1 : 0; twiceAround >= 0; twiceAround--) {
            sx = sy;
            dnfcx = dnfcy;
            dnflx = dnfly;
            dx = dy;

            for (; sx.x < send.x; ++sx.x, ++dnfcx.x, ++dnflx.x, ++dx.x) {
                // Distance to nearest feature in current column.
                UInt32 distPotentialFeature = *dnfcx;

                if (distPotentialFeature == NumericTraits<UInt32>::max()) {
                    // No feature in current column.

                    if (potentialFeatureList.empty()) {
                        // No features to the left either.
                        *dnflx = distPotentialFeature;
                        continue;
                    }

                    potentialFeatureList.move_to_end();
                    typename DnfIterator::MoveX firstFeature =
                        potentialFeatureList.get_current();
                    int deltaX = (dnfcx.x - firstFeature) % w;
                    if (deltaX < 0) {
                        deltaX += w;
                    }
                    distPotentialFeature = _nftDistance(static_cast<UInt32>(deltaX),
                                                        dnfcx(firstFeature - dnfcx.x, 0));
                } else {
                    // First add ourself to the list.
                    potentialFeatureList.push_back(dnfcx.x);

                    // Iterate through the list starting at the right. For each
                    // potential feature, all of the potential features to the left
                    // in the list must be strictly closer. If not delete them from
                    // the list.
                    potentialFeatureList.move_to_end();
                }

                while (potentialFeatureList.has_previous()) {
                    // X coordinate of the predecessor.
                    typename DnfIterator::MoveX previousFeature =
                        potentialFeatureList.get_previous();

                    // Subtract the X coordinates to find out how many
                    // columns to the left of dnfcx previousFeature is.
                    // DeltaX must be positive.
                    // modulo w to consider wraparound condition.
                    int deltaX = (dnfcx.x - previousFeature) % w;
                    if (deltaX < 0) {
                        deltaX += w;
                    }

                    // previousFeature is this far from dnfcx.
                    UInt32 distPreviousFeature =
                        _nftDistance(static_cast<UInt32>(deltaX),
                                     dnfcx(previousFeature - dnfcx.x, 0));

                    if (distPreviousFeature >= distPotentialFeature) {
                        // previousFeature is not a candidate for dnflx
                        // or any dnflx further to the right.
                        potentialFeatureList.erase_previous();
                    } else {
                        // previousFeature is a candidate.
                        potentialFeatureList.move_backwards();
                        distPotentialFeature = distPreviousFeature;
                    }
                }

                // The closest feature to dnflx in columns <= dnflx is the first
                // potential feature in the list.
                *dnflx = distPotentialFeature;

                // Set color of dx to be color of closest feature to the left.
                da.set(dx(potentialFeatureList.get_current() - dnfcx.x, 0), dx);
            }
        }
    }

    // Final pass: calculate the distance to the nearest feature in
    // the same column or any column to the right. If this is smaller
    // than dnflx, Then recolor the pixel to the color of the nearest
    // feature to the right.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) {
            cerr << " 5/6";
        } else {
            cerr << " 4/4";
        }
        cerr.flush();
    }
    sy = src_lowerright;
    send = src_upperleft;
    smidpoint = src_upperleft + Diff2D(w / 2, h / 2);
    dnfcy = dnfColumn->lowerRight();
    dnfly = dnfLeft->lowerRight();
    dy = dest_upperleft + Diff2D(w, h);
    for (; sy.y > send.y;) {
        --sy.y;
        --dnfcy.y;
        --dnfly.y;
        --dy.y;

        // Indicate halfway mark when wraparound is true.
        if (Verbose > VERBOSE_NFT_MESSAGES && wraparound && sy.y == smidpoint.y) {
            cerr << " 6/6";
            cerr.flush();
        }

        potentialFeatureList.clear();

        // If wraparound is true, we must go across the row twice.
        // This takes care of the case when the nearest feature is
        // reached by wrapping around the image.
        for (int twiceAround = wraparound ? 1 : 0; twiceAround >= 0; twiceAround--) {
            sx = sy;
            dnfcx = dnfcy;
            dnflx = dnfly;
            dx = dy;

            for (; sx.x > send.x;) {
                --sx.x;
                --dnfcx.x;
                --dnflx.x;
                --dx.x;

                UInt32 distPotentialFeature = *dnfcx;

                if (distPotentialFeature == NumericTraits<UInt32>::max()) {
                    // No feature in current column.

                    if (potentialFeatureList.empty()) {
                        // No features to the right. Nearest feature
                        // must be to the left.
                        continue;
                    }

                    potentialFeatureList.move_to_end();
                    typename DnfIterator::MoveX firstFeature =
                        potentialFeatureList.get_current();
                    int deltaX = (firstFeature - dnfcx.x) % w;
                    if (deltaX < 0) {
                        deltaX += w;
                    }
                    distPotentialFeature = _nftDistance(static_cast<UInt32>(deltaX),
                                                        dnfcx(firstFeature - dnfcx.x, 0));
                } else {
                    // First add ourself to the list.
                    potentialFeatureList.push_back(dnfcx.x);

                    // Iterate through list and prune as before.
                    potentialFeatureList.move_to_end();
                }

                while (potentialFeatureList.has_previous()) {
                    // X coordinate of the predecessor.
                    typename DnfIterator::MoveX previousFeature =
                        potentialFeatureList.get_previous();

                    // Subtract the X coordinates to find out how many
                    // columns to the right of dnfcx previousFeature is.
                    // DeltaX must be positive.
                    // modulo w to consider wraparound condition.
                    int deltaX = (previousFeature - dnfcx.x) % w;
                    if (deltaX < 0) {
                        deltaX += w;
                    }

                    // previousFeature is this far from dnfcx.
                    UInt32 distPreviousFeature =
                        _nftDistance(static_cast<UInt32>(deltaX),
                                     dnfcx(previousFeature - dnfcx.x, 0));

                    if (distPreviousFeature >= distPotentialFeature) {
                        // previousFeature is not a candidate.
                        potentialFeatureList.erase_previous();
                    } else {
                        // previousFeature is a candidate.
                        potentialFeatureList.move_backwards();
                        distPotentialFeature = distPreviousFeature;
                    }
                }

                // The closest feature on the right is potentialFeature.
                if (*dnflx > distPotentialFeature) {
                    // Following line only necessary for advanced mask generation.
                    //*dnflx = distPotentialFeature;
                    // Recolor dx.
                    da.set(dx(potentialFeatureList.get_current() - dnfcx.x, 0), dx);
                }
            }
        }
    }

    delete dnfColumn;
    delete dnfLeft;
    delete [] lastFeature;
    delete [] foundFirstFeature;
    delete [] lastFeatureDeltaY;

    if (Verbose > VERBOSE_NFT_MESSAGES) {
        cerr << endl;
    }
}


// Version using argument object factories.
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void nearestFeatureTransform(bool wraparound,
                                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                    pair<DestImageIterator, DestAccessor> dest,
                                    typename SrcAccessor::value_type featurelessPixel)
{
    nearestFeatureTransform(wraparound,
                            src.first, src.second, src.third,
                            dest.first, dest.second, featurelessPixel);
}

} // namespace enblend

#endif /* __NEAREST_H__ */

// Local Variables:
// mode: c++
// End:
