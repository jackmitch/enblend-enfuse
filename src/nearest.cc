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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <math.h>

#include "enblend.h"

using namespace std;

extern int Verbose;
extern uint32 OutputWidth;
extern uint32 OutputHeight;

// Region of interest for this operation.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

#define EUCLIDEAN_METRIC

#ifdef EUCLIDEAN_METRIC
    //typedef float dist_t;
    //#define DIST_MAX INFINITY
    //#define DIST_MIN 0.0f
    typedef uint32 dist_t;
    #define DIST_MAX ((dist_t)-1)
    #define DIST_MIN 0
    inline dist_t distance(uint32 deltaX, dist_t *deltaY) {
        return *deltaY + deltaX * deltaX;
    }
    inline dist_t distance(uint32 deltaY) {
        return deltaY * deltaY;
    }
#else
#ifdef CHESSBOARD_METRIC
    typedef uint32 dist_t;
    #define DIST_MAX ((dist_t)-1)
    #define DIST_MIN 0
    inline dist_t distance(uint32 deltaX, dist_t *deltaY) {
        return max(deltaX, *deltaY);
    }
    inline dist_t distance(uint32 deltaY) {
        return deltaY;
    }
#else
#ifdef MANHATTAN_METRIC
    typedef uint32 dist_t;
    #define DIST_MAX ((dist_t)-1)
    #define DIST_MIN 0
    inline dist_t distance(uint32 deltaX, dist_t *deltaY) {
        return deltaX + *deltaY;
    }
    inline dist_t distance(uint32 deltaY) {
        return deltaY;
    }
#endif
#endif
#endif

/** Perform a nearest feature transform on the input image within the region
 *  of interest. For each thinnable pixel, determine if the pixel is closer
 *  to a green pixel or a blue pixel. Make the thinnable pixel the same color
 *  as the closest green or blue pixel.
 */
void nearestFeatureTransform(MaskPixel *mask) {

    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;
    uint32 roiPixels = roiWidth * roiHeight;

    // For each pixel, store the distance to the nearest feature in the same
    // column.
    dist_t *dnfColumn = (dist_t*)malloc(roiPixels * sizeof(dist_t));
    if (dnfColumn == NULL) {
        cerr << "nearestFeatureTransform: malloc failed for dnfColumn" << endl;
        exit(1);
    }

    // For each pixel, store the distance to the nearest feature in the same
    // column or any column to the left.
    dist_t *dnfLeft = (dist_t*)malloc(roiPixels * sizeof(dist_t));
    if (dnfLeft == NULL) {
        cerr << "nearestFeatureTransform: malloc failed for dnfLeft" << endl;
        exit(1);
    }

    // Initialize dnfColumn top-down. Store the distance to the nearest feature
    // in the same column and above us.
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: top-down pass" << endl;
    }
    MaskPixel *firstMaskP = &mask[ROIFirstY * OutputWidth + ROIFirstX];
    for (uint32 x = 0; x < roiWidth; x++) {
        dist_t *dnfColumnP = &dnfColumn[x];
        MaskPixel *maskP = firstMaskP + x;

        // Color of the last feature pixel.
        uint8 lastFeatureG = 0;
        uint8 lastFeatureB = 0;
        // Distance to the last feature pixel in pixels.
        uint32 lastFeatureDeltaY = 0;
        bool foundFirstFeature = false;

        for (uint32 y = 0; y < roiHeight; y++) {
            if (maskP->r == 0) {
                // maskP is a feature pixel.
                *dnfColumnP = DIST_MIN;
                lastFeatureDeltaY = 0;
                lastFeatureG = maskP->g;
                lastFeatureB = maskP->b;
                foundFirstFeature = true;
            } else if (foundFirstFeature) {
                // maskP is not a feature.
                *dnfColumnP = distance(lastFeatureDeltaY);
                maskP->g = lastFeatureG;
                maskP->b = lastFeatureB;
            } else {
                *dnfColumnP = DIST_MAX;
            }

            lastFeatureDeltaY++;

            // Move pointers down one row.
            dnfColumnP += roiWidth;
            maskP += OutputWidth;
        }
    }

    // Initialize dnfColumn bottom-up. Caluclate the distance to the nearest
    // feature in the same column and below us.
    // If this is smaller than the value caluclated in the top-down pass,
    // overwrite that value.
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: bottom-up pass" << endl;
    }
    MaskPixel *lastMaskP = &mask[ROILastY * OutputWidth + ROILastX];
    dist_t *lastDNFColumnP = &dnfColumn[roiWidth * roiHeight - 1];
    for (uint32 x = 0; x < roiWidth; x++) {
        dist_t *dnfColumnP = lastDNFColumnP - x;
        MaskPixel *maskP = lastMaskP - x;

        // Color of the last feature pixel.
        uint8 lastFeatureG = 0;
        uint8 lastFeatureB = 0;
        // Distance to the last feature pixel in pixels.
        uint32 lastFeatureDeltaY = 0;
        bool foundFirstFeature = false;

        for (uint32 y = 0; y < roiHeight; y++) {
            if (maskP->r == 0) {
                // maskP is a feature pixel.
                //*dnfColumnP = DIST_MIN; don't need to do this again.
                lastFeatureDeltaY = 0;
                lastFeatureG = maskP->g;
                lastFeatureB = maskP->b;
                foundFirstFeature = true;
            } else if (foundFirstFeature) {
                // maskP is not a feature.
                dist_t distLastFeature = distance(lastFeatureDeltaY);
                // If last feature is closer than nearest feature above,
                // change distance and color to match last feature.
                if (distLastFeature < *dnfColumnP) {
                    *dnfColumnP = distLastFeature;
                    maskP->g = lastFeatureG;
                    maskP->b = lastFeatureB;
                }
            }

            lastFeatureDeltaY++;

            // Move pointers up one row.
            dnfColumnP -= roiWidth;
            maskP -= OutputWidth;
        }
    }

    //size_t maxListSize = 0;

    // Calculate dnfLeft for each pixel.
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: left-right pass" << endl; //... ";
        //cout.flush();
    }
    for (uint32 y = 0; y < roiHeight; y++) {
        dist_t *dnfLeftP = &dnfLeft[y * roiWidth];
        dist_t *dnfColumnP = &dnfColumn[y * roiWidth];
        MaskPixel *maskP = firstMaskP + (y * OutputWidth);

        // List of dnfColumnP's on the left that might be the closest features
        // to the current dnfColumnP.
        list<dist_t*> potentialFeatureList;

        for (uint32 x = 0; x < roiWidth; x++) {
            // First add ourself to the list.
            potentialFeatureList.push_back(dnfColumnP);

            // Iterate through the list starting at the right. For each
            // potential feature, all of the potential features to the left
            // in the list must be strictly closer. If not delete them from
            // the list.
            list<dist_t*>::iterator potentialFeature =
                    --(potentialFeatureList.end());
            // The last potential feature is dnfColumnP, just added above.
            // That is in the current column so the distance to that feature
            // is simply *dnfColumnP.
            dist_t distPotentialFeature = *dnfColumnP;
            while (potentialFeature != potentialFeatureList.begin()) {
                // Make an iterator that points to the predecessor.
                list<dist_t*>::iterator previousFeature = potentialFeature;
                previousFeature--;

                // previousFeature is this many columns to the left of (x,y).
                uint32 deltaX = dnfColumnP - *previousFeature;

                // previousFeature is this far from (x,y).
                dist_t distPreviousFeature = distance(deltaX, *previousFeature);

                if (distPreviousFeature >= distPotentialFeature) {
                    // previousFeature is not a candidate for dnfLeftP
                    // or anything further to the right.
                    potentialFeatureList.erase(previousFeature);
                } else {
                    // previousFeature is a candidate.
                    potentialFeature = previousFeature;
                    distPotentialFeature = distPreviousFeature;
                }
            }

            // The closest feature to (x,y) in columns <= x is the first
            // potential feature in the list.
            *dnfLeftP = distPotentialFeature;

            // Set color of maskP to be color of closest feature to the left.
            MaskPixel *maskPLeft = maskP - (dnfColumnP - *potentialFeature);
            maskP->g = maskPLeft->g;
            maskP->b = maskPLeft->b;

            // Move pointers right one column.
            dnfLeftP++;
            dnfColumnP++;
            maskP++;

            //maxListSize = max(maxListSize, potentialFeatureList.size());
        }
    }
    //if (Verbose > 0) {
    //    cout << "max feature list size=" << maxListSize << endl;
    //}
    //maxListSize = 0;

    // Final pass: calculate the distance to the nearest feature in the same
    // column or any column to the right. If this is smaller than dnfLeftP,
    // Then recolor the pixel to the color of the nearest feature to the right.
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: right-left pass" << endl; //... ";
        //cout.flush();
    }
    dist_t *lastDNFLeftP = &dnfLeft[roiWidth * roiHeight - 1];
    for (uint32 y = 0; y < roiHeight; y++) {
        dist_t *dnfColumnP = lastDNFColumnP - (y * roiWidth);
        dist_t *dnfLeftP = lastDNFLeftP - (y * roiWidth);
        MaskPixel *maskP = lastMaskP - (y * OutputWidth);

        // List of dnfColumnP's on the right that might be the closest features
        // to the current dnfColumnP.
        list<dist_t*> potentialFeatureList;

        for (uint32 x = 0; x < roiWidth; x++) {
            // First add ourself to the list.
            potentialFeatureList.push_back(dnfColumnP);

            // Iterate through list and prune as before.
            list<dist_t*>::iterator potentialFeature =
                    --(potentialFeatureList.end());
            dist_t distPotentialFeature = *dnfColumnP;
            while (potentialFeature != potentialFeatureList.begin()) {
                // Iterator that points to predecessor.
                list<dist_t*>::iterator previousFeature = potentialFeature;
                previousFeature--;

                // previousFeature is this many columns to the right of (x,y);
                uint32 deltaX = *previousFeature - dnfColumnP;

                // previousFeature is this far from (x,y);
                dist_t distPreviousFeature = distance(deltaX, *previousFeature);

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
            if (*dnfLeftP > distPotentialFeature) {
                // Recolor maskP.
                MaskPixel *maskPRight = maskP + (*potentialFeature - dnfColumnP);
                maskP->g = maskPRight->g;
                maskP->b = maskPRight->b;
            }

            // Move pointers left one column.
            dnfLeftP--;
            dnfColumnP--;
            maskP--;

            //maxListSize = max(maxListSize, potentialFeatureList.size());
        }
    }
    //if (Verbose > 0) {
    //    cout << "max feature list size=" << maxListSize << endl;
    //}

    free(dnfColumn);
    free(dnfLeft);

    return;
}

