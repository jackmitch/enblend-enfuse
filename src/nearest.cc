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
    typedef float dist_t;
    #define DIST_MAX INFINITY
    #define DIST_MIN 0.0f
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

    dist_t *g = (dist_t*)malloc(roiPixels * sizeof(dist_t));
    if (g == NULL) {
        cerr << "nearestFeatureTransform: malloc failed for g" << endl;
        exit(1);
    }
    dist_t *d = (dist_t*)malloc(roiPixels * sizeof(dist_t));
    if (d == NULL) {
        cerr << "nearestFeatureTransform: malloc failed for d" << endl;
        exit(1);
    }

    // Top-down g initialization
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: preinit top-down." << endl;
    }
    for (uint32 x = 0; x < roiWidth; x++) {
        dist_t *gPixel = &g[x];
        uint32 maskX = x + ROIFirstX;
        uint8 lastG = 0;
        uint8 lastB = 0;
        uint32 lastY = 0;
        bool foundFirstFeature = false;
        for (uint32 y = 0; y < roiHeight; y++) {
            uint32 maskY = y + ROIFirstY;

            MaskPixel *maskPixel = &mask[maskY * OutputWidth + maskX];

            if (maskPixel->r == 0) {
                // maskPixel is a feature pixel.
                *gPixel = DIST_MIN;
                lastY = y;
                lastG = maskPixel->g;
                lastB = maskPixel->b;
                foundFirstFeature = true;
            } else if (foundFirstFeature) {
                // maskPixel is not a feature.
                *gPixel = distance(y - lastY);
                maskPixel->g = lastG;
                maskPixel->b = lastB;
            } else {
                *gPixel = DIST_MAX;
            }

            gPixel += roiWidth;
        }
    }

    // Bottom-up g initialization
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: preinit bottom-up." << endl;
    }
    for (uint32 x = 1; x <= roiWidth; x++) {
        dist_t *gPixel = &g[roiPixels - x];
        uint32 maskX = ROILastX - x + 1;
        uint8 lastG = 0;
        uint8 lastB = 0;
        uint32 lastY = 0;
        bool foundFirstFeature = false;
        for (uint32 y = 0; y < roiHeight; y++) {
            uint32 maskY = ROILastY - y;

            MaskPixel *maskPixel = &mask[maskY * OutputWidth + maskX];

            if (maskPixel->r == 0) {
                // maskPixel is a feature pixel.
                lastY = y;
                lastG = maskPixel->g;
                lastB = maskPixel->b;
                foundFirstFeature = true;
            } else if (foundFirstFeature) {
                // maskPixel is not a feature.
                dist_t dist = distance(y - lastY);
                if (dist < *gPixel) {
                    *gPixel = dist;
                    maskPixel->g = lastG;
                    maskPixel->b = lastB;
                }
            }

            gPixel -= roiWidth;
        }
    }

    // Left-to-right transform
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: LR pass." << endl;
    }
    for (uint32 y = 0; y < roiHeight; y++) {
        uint32 maskY = ROIFirstY + y;
        dist_t *gPixel = &g[y * roiWidth];
        dist_t *dPixel = &d[y * roiWidth];
        list<dist_t*> featureList;

        for (uint32 x = 0; x < roiWidth; x++) {
            uint32 maskX = ROIFirstX + x;
            MaskPixel *maskPixel = &mask[maskY * OutputWidth + maskX];

            // First add ourself to the list.
            featureList.push_back(gPixel);

            // Iterate backwards through list, prune
            list<dist_t*>::iterator featureIterator = featureList.end();
            featureIterator--;
            dist_t featureIteratorDistance = **featureIterator;
            while (featureIterator != featureList.begin()) {
                list<dist_t*>::iterator previous = featureIterator;
                previous--;
                // previous is this many columns to the left of the origin.
                uint32 colDiff = gPixel - *previous;
                // previous feature is this far from origin.
                dist_t dist = distance(colDiff, *previous);
                if (dist >= featureIteratorDistance) {
                    // previous is no candidate.
                    featureList.erase(previous);
                } else {
                    // previous is a candidate.
                    featureIterator = previous;
                    featureIteratorDistance = dist;
                }
            }

            // Now choose the feature
            dist_t *gChoose = featureList.front();
            int32 colDiff = gChoose - gPixel;
            *dPixel = distance(colDiff, gChoose);
            *gPixel = distance(colDiff, gChoose);
            maskPixel->g = maskPixel[colDiff].g;
            maskPixel->b = maskPixel[colDiff].b;

            gPixel++;
            dPixel++;
        }
    }

    // right-to-left transform
    if (Verbose > 0) {
        cout << "nearestFeatureTransform: RL pass." << endl;
    }
    for (uint32 y = 0; y < roiHeight; y++) {
        uint32 maskY = ROIFirstY + y;
        dist_t *gPixel = &g[(y+1) * roiWidth - 1];
        dist_t *dPixel = &d[(y+1) * roiWidth - 1];
        list<dist_t*> featureList;

        for (uint32 x = 0; x < roiWidth; x++) {
            uint32 maskX = ROILastX - x;
            MaskPixel *maskPixel = &mask[maskY * OutputWidth + maskX];

            // First add ourself to the list.
            featureList.push_back(gPixel);

            // Iterate backwards through list, prune
            list<dist_t*>::iterator featureIterator = featureList.end();
            featureIterator--;
            dist_t featureIteratorDistance = **featureIterator;
            while (featureIterator != featureList.begin()) {
                list<dist_t*>::iterator previous = featureIterator;
                previous--;
                uint32 colDiff = *previous - gPixel;
                dist_t dist = distance(colDiff, *previous);
                if (dist >= featureIteratorDistance) {
                    // previous is no candidate.
                    featureList.erase(previous);
                } else {
                    // previous is a candidate.
                    featureIterator = previous;
                    featureIteratorDistance = dist;
                }
            }

            // Now choose the feature
            dist_t *gChoose = featureList.front();
            int32 colDiff = gChoose - gPixel;
            if (*dPixel > distance(colDiff, gChoose)) {
                maskPixel->g = maskPixel[colDiff].g;
                maskPixel->b = maskPixel[colDiff].b;
            }

            gPixel--;
            dPixel--;
        }
    }

    free(g);
    free(d);
    return;
}

