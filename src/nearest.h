/*
 * Copyright (C) 2004-2007 Andrew Mihal
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

/** Data structure for potentialFeatureStack.
 *  Contribution from Fulvio Senore.
 *  Fast insert and delete, avoiding dynamic memory allocation.
 */
template <typename T>
class FeatureList {

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

    FeatureList(int size) {
        array = new Node[size];
        firstUnused = array;
        last = NULL;
        iterator = NULL;
    }

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
        return (iterator->prev != NULL);
    }

    void erase_previous() {
        iterator->prev = iterator->prev->prev;
    }

    T get_current() {
        return iterator->value;
    }

    T get_previous() {
        return iterator->prev->value;
    }

    void move_backwards() {
        iterator = iterator->prev;
    }

    bool empty() {
        return (last == NULL);
    }

};

template<class FeatureType> class Chain {
public:
    struct Node {
        bool active;
        FeatureType f;
        int x;
        int y;
        Node* left;
        Node* right;
    };

    Chain(int w) {
        nodeArray = new Node[w];
        activeNodes = 0;
        leftmostNode = NULL;
        rightmostNode = NULL;
        for (int i = 0; i < w; ++i) {
            nodeArray[i].active = false;
            nodeArray[i].x = i;
            nodeArray[i].left = NULL;
            nodeArray[i].right = NULL;
        }
    }

    ~Chain() {
        delete [] nodeArray;
    }

    void insertAfter(Node* n, int x, int y, FeatureType f) {
        if (nodeArray[x].active) {
            nodeArray[x].y = y;
            nodeArray[x].f = f;
        }
        else if (n == NULL) {
            // Insert at beginning
            nodeArray[x].active = true;
            nodeArray[x].f = f;
            nodeArray[x].y = y;
            nodeArray[x].left = NULL;
            nodeArray[x].right = leftmostNode;
            if (leftmostNode != NULL) {
                leftmostNode->left = &nodeArray[x];
            }
            leftmostNode = &nodeArray[x];
            if (rightmostNode == NULL) {
                rightmostNode = &nodeArray[x];
            }
            ++activeNodes;
        } else {
            // Insert after given node.
            nodeArray[x].active = true;
            nodeArray[x].f = f;
            nodeArray[x].y = y;
            nodeArray[x].left = n;
            nodeArray[x].right = n->right;
            if (n->right == NULL) {
                rightmostNode = &nodeArray[x];
            } else {
                nodeArray[x].right->left = &nodeArray[x];
            }
            n->right = &nodeArray[x];
            ++activeNodes;
        }
    }

    void erase(Node* n) {
        n->active = false;
        if (n->left == NULL) {
            leftmostNode = n->right;
        } else {
            n->left->right = n->right;
        }
        if (n->right == NULL) {
            rightmostNode = n->left;
        } else {
            n->right->left = n->left;
        }
        n->left = NULL;
        n->right = NULL;
        --activeNodes;
    }

    // Compute the y-coordinate of the intersection of the
    // perpendicular bisector of pq with the vertical line at x.
    pair<double, bool> pbY(Node* p, double x) {
        Node* q = p->right;
        assert(q != NULL);
        const double px = p->x;
        const double py = p->y;
        const double qx = q->x;
        const double qy = q->y;
        if (q->y == p->y) return std::make_pair<double, bool>(0.0, false);
        const double midX = (qx + px) / 2.0;
        const double midY = (qy + py) / 2.0;
        const double m = (qy - py) / (qx - px);
        const double mPB = -1.0 / m;
        // y = mPB * (x - midX) + midY
        return std::make_pair<double, bool>((mPB * (x - midX) + midY), true);
    }

    // Compute the x-coordinate of the intersection of the
    // perpendicular bisector of pq with the horizontal line at y.
    double pbX(Node* p, double y) {
        Node* q = p->right;
        assert(q != NULL);
        const double px = p->x;
        const double py = p->y;
        const double qx = q->x;
        const double qy = q->y;
        const double midX = (qx + px) / 2.0;
        if (q->y == p->y) return midX;
        const double midY = (qy + py) / 2.0;
        const double m = (qy - py) / (qx - px);
        const double mPB = -1.0 / m;
        // y = mPB * (x - midX) + midY
        // y - midY = mPB * (x - midX)
        // (y - midY) / mPB = x - midX
        // x = (y - midY) / mPB + midX
        return ((y - midY) / mPB + midX);
    }

    // Compute the y-coordinate of the intersection of the
    // perpendicular bisectors of pq and qr.
    pair<double, bool> intersectionY(Node* p) {
        Node* q = p->right;
        assert(q != NULL);
        Node* r = q->right;
        assert(r != NULL);
        const int px = p->x;
        const int py = p->y;
        const int qx = q->x;
        const int qy = q->y;
        const int rx = r->x;
        const int ry = r->y;

        const double midPQX = (qx + px) / 2.0;
        if (qy == py) return pbY(q, midPQX);

        const double midQRX = (rx + qx) / 2.0;
        if (qy == ry) return pbY(p, midQRX);

        const double midPQY = (qy + py) / 2.0;
        const double slopePQ = static_cast<double>(qy - py) / static_cast<double>(qx - px);
        const double pbSlopePQ = -1.0 / slopePQ;

        const double midQRY = (ry + qy) / 2.0;
        const double slopeQR = static_cast<double>(ry - qy) / static_cast<double>(rx - qx);
        const double pbSlopeQR = -1.0 / slopeQR;

        if (slopePQ == slopeQR) return std::make_pair<double, bool>(0.0, false);

        // pbPQ: y = pbSlopePQ * (x - midPQX) + midPQY;
        // pbQR: y = pbSlopeQR * (x - midQRX) + midQRY;
        //       y = (pbSlopeQR * x) - (pbSlopeQR * midQRX) + midQRY
        //       pbSlopeQR * x = y + (pbSlopeQR * midQRX) - midQRY
        //       x = (y + (pbSlopeQR * midQRX) - midQRY) / pbSlopeQR
        // substituting into pbPQ:
        //       y = pbSlopePQ * (y + (pbSlopeQR * midQRX) - midQRY) / pbSlopeQR - pbSlopePQ*midPQX + midPQY;
        //         = (pbSlopePQ/pbSlopeQR) * (y + pbSlopeQR*midQRX - midQRY) - pbSlopePQ*midPQX + midPQY;
        //       y - (pbSlopePQ/pbSlopeQR)*y = (pbSlopePQ/pbSlopeQR)*(pbSlopeQR*midQRX - midQRY) - pbSlopePQ*midPQX + midPQY
        //       y * (1 - pbSlopePQ/pbSlopeQR) = ...
        //       y = ... / (1 - pbSlopePQ/pbSlopeQR)
        const double iy = ((pbSlopePQ/pbSlopeQR)*(pbSlopeQR*midQRX - midQRY) - pbSlopePQ*midPQX + midPQY)
                          / (1.0 - pbSlopePQ / pbSlopeQR);
        return std::make_pair<double, bool>(iy, true);
    }

    Node* operator[](int idx) {
        return &(nodeArray[idx]);
    }

    Node* nodeArray;
    int activeNodes;
    Node* leftmostNode;
    Node* rightmostNode;
};

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void nearestFeatureTransformExact(
        SrcImageIterator src_upperleft,
        SrcImageIterator src_lowerright,
        SrcAccessor sa,
        DestImageIterator dest_upperleft,
        DestAccessor da,
        typename SrcAccessor::value_type featurelessPixel) {

    typedef typename EnblendNumericTraits<UInt32>::ImageType DNFImage;
    typedef typename DNFImage::traverser DnfIterator;
    typedef typename SrcAccessor::value_type SrcValueType;

    SrcImageIterator sx, sy, send;
    DnfIterator dnfAboveX, dnfAboveY;
    DestImageIterator dx, dy;

    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    DNFImage dnfAbove(w, h, NumericTraits<UInt32>::max());
    Chain<SrcValueType> chain(w);

    // Top-down pass
    sy = src_upperleft;
    send = src_lowerright;
    dy = dest_upperleft;
    dnfAboveY = dnfAbove.upperLeft();

    for (int yIndex = 0; sy.y < send.y; ++yIndex, ++sy.y, ++dy.y, ++dnfAboveY.y) {
        // Insert new features into chain.
        sx = sy;
        typename Chain<SrcValueType>::Node* lastActiveNode = NULL;
        for (int xIndex = 0; sx.x < send.x; ++xIndex, ++sx.x) {
            SrcValueType pixel = sa(sx);
            if (pixel != featurelessPixel) {
                //cout << "new feature xIndex=" << xIndex << " yIndex=" << yIndex << " f=" << (int)pixel << endl;
                chain.insertAfter(lastActiveNode, xIndex, yIndex, pixel);
            }
            if (chain[xIndex]->active) lastActiveNode = chain[xIndex];
        }

        //cout << "after insert: ";
        //for (lastActiveNode = chain.leftmostNode; lastActiveNode != NULL; lastActiveNode = lastActiveNode->right) {
        //    cout << "[" << lastActiveNode->x << "," << lastActiveNode->y << ": " << (int)lastActiveNode->f << "] ";
        //}
        //cout << endl;

        // Remove leftmost extreme features
        while (chain.activeNodes > 1) {
            // If intersection of perpendicular bisector of pq and x=0 is < yIndex, delete p.
            pair<double, bool> r = chain.pbY(chain.leftmostNode, 0);
            if (r.second && (r.first < yIndex) && (r.first > chain.leftmostNode->y)) chain.erase(chain.leftmostNode);
            else break;
        }

        //cout << "after leftmost: ";
        //for (lastActiveNode = chain.leftmostNode; lastActiveNode != NULL; lastActiveNode = lastActiveNode->right) {
        //    cout << "[" << lastActiveNode->x << "," << lastActiveNode->y << ": " << (int)lastActiveNode->f << "] ";
        //}

        // Remove rightmost extreme features
        while (chain.activeNodes > 1) {
            // If intersection of perpendicular bisector of pq and x=w-1 is <= yIndex, delete q.
            pair<double, bool> r = chain.pbY(chain.rightmostNode->left, w-1);
            if (r.second && (r.first <= yIndex) && (r.first > chain.rightmostNode->y)) chain.erase(chain.rightmostNode);
            else break;
        }

        //cout << "after rightmost: ";
        //for (lastActiveNode = chain.leftmostNode; lastActiveNode != NULL; lastActiveNode = lastActiveNode->right) {
        //    cout << "[" << lastActiveNode->x << "," << lastActiveNode->y << ": " << (int)lastActiveNode->f << "] ";
        //}

        // Remove internal features
        typename Chain<SrcValueType>::Node* p = chain.leftmostNode;
        while (chain.activeNodes > 2) {
            typename Chain<SrcValueType>::Node* q = p->right;
            typename Chain<SrcValueType>::Node* r = q->right;
            if (r == NULL) break;
            pair<double, bool> iy = chain.intersectionY(p);
            if (iy.second && (iy.first < yIndex) && (iy.first > q->y)) {
                // remove q and backtrack
                chain.erase(q);
                if (p->left != NULL) p = p->left;
            } else {
                // advance
                p = q;
            }
        }

        //cout << "after internal: ";
        //for (lastActiveNode = chain.leftmostNode; lastActiveNode != NULL; lastActiveNode = lastActiveNode->right) {
        //    cout << "[" << lastActiveNode->x << "," << lastActiveNode->y << ": " << (int)lastActiveNode->f << "] ";
        //}
        //cout << endl;

        // assign
        dx = dy;
        dnfAboveX = dnfAboveY;
        lastActiveNode = chain.leftmostNode;
        int xIndex = 0;
        while (lastActiveNode != NULL) {
            int xLimit = (lastActiveNode->right) ? static_cast<int>(floor(chain.pbX(lastActiveNode, yIndex))) : w;
            //cout << "[" << lastActiveNode->x << "," << lastActiveNode->y << ": " << (int)lastActiveNode->f << "] xLimit=" << xLimit << endl;
            for (; xIndex < w && xIndex <= xLimit; ++xIndex) {
                dx(xIndex, 0) = lastActiveNode->f;
                dnfAboveX(xIndex, 0) = (xIndex - lastActiveNode->x) * (xIndex - lastActiveNode->x) + (yIndex - lastActiveNode->y) * (yIndex - lastActiveNode->y);
            }
            lastActiveNode = lastActiveNode->right;
        }

    } /* end of top-down pass */
/*
    sy = src_lowerright;
    send = src_upperleft;
    dy = dest_upperleft + Diff2D(w, h);
    dnfAboveY = dnfAbove.lowerRight();

    for (int yIndex = h; sy.y > send.y;) {
        --yIndex;
        --sy.y;
        --dy.y;
        --dnfAboveY.y;

        sx = sy;
        typename Chain<SrcValueType>::Node* lastActiveNode = NULL;
        for (int xIndex = w; sx.x > send.x;) {
            --xIndex;
            --sx.x;
*/

}
    
/** Compute the nearest feature transform.
 *  A non-zero pixel in the src image is considered a feature.
 *  Each pixel in the dest image is given the value of the nearest feature
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
        typename SrcAccessor::value_type featurelessPixel) {

    typedef typename EnblendNumericTraits<UInt32>::ImageType DNFImage;
    typedef typename DNFImage::traverser DnfIterator;
    typedef typename SrcAccessor::value_type SrcValueType;

    SrcImageIterator sx, sy, send, smidpoint;
    DnfIterator dnfcx, dnfcy;
    DnfIterator dnflx, dnfly;
    DestImageIterator dx, dy;

    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;

    // Distance to the nearest feature in the current column.
    DNFImage *dnfColumn = new DNFImage(w, h);
    // Distance to the nearest feature in the current column, or any
    // column to the left of this column.
    DNFImage *dnfLeft = new DNFImage(w, h);

    // Data structures for initializing dnfColumn.
    // These let us initialize all of the columns in one pass
    // over the rows of the image. Cache-friendly.
    SrcValueType* lastFeature = new SrcValueType[w];
    bool* foundFirstFeature = new bool[w];
    UInt32* lastFeatureDeltaY = new UInt32[w];

    // Initialize dnfColumn top-down. Store the distance to the nearest feature
    // in the same column and above us.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) cout << "Creating blend mask: 1/6";
        else cout << "Creating blend mask: 1/4";
        cout.flush();
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
            }
            else if (foundFirstFeature[xIndex]) {
                // Source pixel is not a feature.
                *dnfcx = _nftDistance(lastFeatureDeltaY[xIndex]);
                da.set(lastFeature[xIndex], dx);
            }
            else {
                *dnfcx = NumericTraits<UInt32>::max();
            }
            ++lastFeatureDeltaY[xIndex];
        }
    }

    // Initialize dnfColumn bottom-up. Caluclate the distance to the nearest
    // feature in the same column and below us.
    // If this is smaller than the value caluclated in the top-down pass,
    // overwrite that value.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) cout << " 2/6";
        else cout << " 2/4";
        cout.flush();
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

        for (int xIndex = w-1; sx.x > send.x; --xIndex) {
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
            }
            else if (foundFirstFeature[xIndex]) {
                // Source pixel is not a feature
                UInt32 distLastFeature = _nftDistance(lastFeatureDeltaY[xIndex]);
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
    FeatureList<typename DnfIterator::MoveX> potentialFeatureList((wraparound) ? w*2 : w);

    // Calculate dnfLeft for each pixel.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) cout << " 3/6";
        else cout << " 3/4";
        cout.flush();
    }
    sy = src_upperleft;
    send = src_lowerright;
    smidpoint = src_upperleft + Diff2D(w/2, h/2);
    dnfcy = dnfColumn->upperLeft();
    dnfly = dnfLeft->upperLeft();
    dy = dest_upperleft;
    for (; sy.y < send.y; ++sy.y, ++dnfcy.y, ++dnfly.y, ++dy.y) {

        // Indicate halfway mark when wraparound is true.
        if (Verbose > VERBOSE_NFT_MESSAGES
                && wraparound && (sy.y == smidpoint.y)) {
            cout << " 4/6";
            cout.flush();
        }

        potentialFeatureList.clear();

        // If wraparound is true, we must go across the row twice.
        // This takes care of the case when the nearest feature is reached by
        // wrapping around the image.
        for (int twiceAround = (wraparound?1:0); twiceAround >= 0; twiceAround--) {
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
                    typename DnfIterator::MoveX firstFeature = potentialFeatureList.get_current();
                    int deltaX = (dnfcx.x - firstFeature) % w;
                    if (deltaX < 0) deltaX += w;
                    distPotentialFeature = _nftDistance((UInt32)deltaX, dnfcx(firstFeature - dnfcx.x, 0));
                }
                else {
                    // First add ourself to the list.
                    potentialFeatureList.push_back(dnfcx.x);

                    // Iterate throught the list starting at the right. For each
                    // potential feature, all of the potential features to the left
                    // in the list must be strictly closer. If not delete them from
                    // the list.
                    potentialFeatureList.move_to_end();
                }

                while (potentialFeatureList.has_previous()) {
                    // X coordinate of the predecessor.
                    typename DnfIterator::MoveX previousFeature = potentialFeatureList.get_previous();

                    // Subtract the X coordinates to find out how many
                    // columns to the left of dnfcx previousFeature is.
                    // DeltaX must be positive.
                    // modulo w to consider wraparound condition.
                    int deltaX = (dnfcx.x - previousFeature) % w;
                    if (deltaX < 0) deltaX += w;

                    // previousFeature is this far from dnfcx.
                    UInt32 distPreviousFeature =
                            _nftDistance((UInt32)deltaX, dnfcx(previousFeature - dnfcx.x, 0));

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
                da.set(dx((potentialFeatureList.get_current() - dnfcx.x), 0), dx);
            }
        }
    }

    // Final pass: calculate the distance to the nearest feature in the same
    // column or any column to the right. If this is smaller than dnflx,
    // Then recolor the pixel to the color of the nearest feature to the right.
    if (Verbose > VERBOSE_NFT_MESSAGES) {
        if (wraparound) cout << " 5/6";
        else cout << " 4/4";
        cout.flush();
    }
    sy = src_lowerright;
    send = src_upperleft;
    smidpoint = src_upperleft + Diff2D(w/2, h/2);
    dnfcy = dnfColumn->lowerRight();
    dnfly = dnfLeft->lowerRight();
    dy = dest_upperleft + Diff2D(w, h);
    for (; sy.y > send.y;) {
        --sy.y;
        --dnfcy.y;
        --dnfly.y;
        --dy.y;

        // Indicate halfway mark when wraparound is true.
        if (Verbose > VERBOSE_NFT_MESSAGES
                && wraparound && (sy.y == smidpoint.y)) {
            cout << " 6/6";
            cout.flush();
        }

        potentialFeatureList.clear();

        // If wraparound is true, we must go across the row twice.
        // This takes care of the case when the nearest feature is reached by
        // wrapping around the image.
        for (int twiceAround = (wraparound?1:0); twiceAround >= 0; twiceAround--) {
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
                        // No features to the right. Nearest feature must be to the left.
                        continue;
                    }

                    potentialFeatureList.move_to_end();
                    typename DnfIterator::MoveX firstFeature = potentialFeatureList.get_current();
                    int deltaX = (firstFeature - dnfcx.x) % w;
                    if (deltaX < 0) deltaX += w;
                    distPotentialFeature = _nftDistance((UInt32)deltaX, dnfcx(firstFeature - dnfcx.x, 0));
                }
                else {
                    // First add ourself to the list.
                    potentialFeatureList.push_back(dnfcx.x);

                    // Iterate through list and prune as before.
                    potentialFeatureList.move_to_end();
                }

                while (potentialFeatureList.has_previous()) {
                    // X coordinate of the predecessor.
                    typename DnfIterator::MoveX previousFeature = potentialFeatureList.get_previous();

                    // Subtract the X coordinates to find out how many
                    // columns to the right of dnfcx previousFeature is.
                    // DeltaX must be positive.
                    // modulo w to consider wraparound condition.
                    int deltaX = (previousFeature - dnfcx.x) % w;
                    if (deltaX < 0) deltaX += w;

                    // previousFeature is this far from dnfcx.
                    UInt32 distPreviousFeature =
                            _nftDistance((UInt32)deltaX, dnfcx(previousFeature - dnfcx.x, 0));

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
                    da.set(dx((potentialFeatureList.get_current() - dnfcx.x), 0), dx);
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
        cout << endl;
    }

    return;
};

// Version using argument object factories.
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void nearestFeatureTransform(bool wraparound,
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<DestImageIterator, DestAccessor> dest,
        typename SrcAccessor::value_type featurelessPixel) {

    nearestFeatureTransform(wraparound,
            src.first, src.second, src.third,
            dest.first, dest.second, featurelessPixel);

};

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void nearestFeatureTransformExact(
        triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<DestImageIterator, DestAccessor> dest,
        typename SrcAccessor::value_type featurelessPixel) {

    nearestFeatureTransformExact(
            src.first, src.second, src.third,
            dest.first, dest.second, featurelessPixel);

};

} // namespace enblend

#endif /* __NEAREST_H__ */
