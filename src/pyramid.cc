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
#include <stdlib.h>
#include <tiffio.h>
#include <stdio.h>

#include "enblend.h"

using namespace std;

extern int Verbose;
extern bool Wraparound;
extern uint32 OutputWidth;
extern uint32 OutputHeight;
extern uint16 PlanarConfig;
extern uint16 Photometric;

// Union bounding box.
extern uint32 UBBFirstX;
extern uint32 UBBLastX;
extern uint32 UBBFirstY;
extern uint32 UBBLastY;

// Region of interest for this operation.
extern uint32 ROIFirstX;
extern uint32 ROILastX;
extern uint32 ROIFirstY;
extern uint32 ROILastY;

// Gaussian filter coefficients
static const double A = 0.4;
static const double W[] = {0.25 - A / 2.0, 0.25, A, 0.25, 0.25 - A / 2.0};
static const uint32 A100 = 40;
static const uint32 W100[] = {25 - A100 / 2, 25, A100, 25, 25 - A100 / 2};

/** Calculate the half width of a level n filter, taking into account
 *  pixel precision and rounding method.
 */
uint32 filterHalfWidth(uint32 level, uint32 maxPixelValue) {

    // This is the arithmetic half width (true for level  > 0).
    uint32 length = 1 + (1 << level);

    // Use internal LPPixel precision data type.
    int16 *f = (int16*)calloc(length, sizeof(int16));
    if (f == NULL) {
        cerr << endl
             << "enblend: out of memory (in filterHalfWidth for f)" << endl;
        exit(1);
    }

    // input f(x) is the step function u(-x)
    f[0] = maxPixelValue;

    for (uint32 l = 1; l <= level; l++) {
        // sample 0 from level l-1
        double pZero = f[0];
        // sample 1 from level l-1
        double pOne = f[1 << (l-1)];

        // sample 0 on level l
        double nZero = (pZero * W[2]) + (pOne * W[3])
                + (maxPixelValue * W[0]) + (maxPixelValue * W[1]);
        f[0] = (int16)rint(nZero);

        // sample 1 on level l
        double nOne = (pZero * W[0]) + (pOne * W[1]);
        f[1 << l] = (int16)rint(nOne);

        // remaining samples on level l are zero.

        // If sample 1 was rounded down to zero, then sample 1 on
        // level l-1 is the rightmost nonzero value.
        if (f[1 << l] == 0) {
            free(f);
            // return the index of the rightmost nonzero value.
            return (1 << (l-1));
        }
    }

    // Else there is no round-to-zero issue.
    free(f);
    return (length - 1);
}

/** The Burt & Adelson Expand operation.
 *  Expand in, and either add or subtract from out.
 */
void expand(LPPixel *in, uint32 inW, uint32 inH,
        LPPixel *out, uint32 outW, uint32 outH,
        bool add, bool wraparound) {

    for (uint32 outY = 0; outY < outH; outY++) {
        for (uint32 outX = 0; outX < outW; outX++) {

            double r = 0.0;
            double g = 0.0;
            double b = 0.0;
            uint32 totalContrib = 0;

            for (int32 m = 0; m < 5; m++) {
                int32 inX = (int32)outX - (m-2);

                if (wraparound) {
                    // Boundary condition - wrap around the image.
                    if (inX < 0) inX = (int32)outW + inX;
                    if (inX >= (int32)outW) inX = inX - (int32)outW;
                }

                // Skip non-integral values of inX index.
                if (inX & 1 == 1) continue;

                inX = inX >> 1;

                if (!wraparound) {
                    // Boundary condition - replicate first and last column.
                    if (inX >= (int32)inW) inX = (int32)inW - 1;
                    if (inX < 0) inX = 0;
                }

                //if (outW == (ROILastX - ROIFirstX + 1) && (outX == 0 || outX == outW - 1 || outX == outW - 2)) {
                //    cout << "outX=" << outX << " inX=" << inX
                //         << " m=" << m << " outW=" << outW << endl;
                //}

                for (int32 n = 0; n < 5; n++) {
                    // Skip non-integral values of inY index.
                    if (((int32)outY - (n-2)) & 1 == 1) continue;

                    int32 inY = ((int32)outY - (n-2)) >> 1;

                    // Boundary condition: replicate top and bottom rows.
                    if (inY >= (int32)inH) inY = inH - 1;
                    if (inY < 0) inY = 0;

                    LPPixel *inPixel = &(in[inY * inW + inX]);

                    r += W[m] * W[n] * inPixel->r;
                    g += W[m] * W[n] * inPixel->g;
                    b += W[m] * W[n] * inPixel->b;
                    totalContrib += W100[m] * W100[n];
                }
            }

            if (add) {
                out->r += (int16)rint(r / (totalContrib / 10000.0));
                out->g += (int16)rint(g / (totalContrib / 10000.0));
                out->b += (int16)rint(b / (totalContrib / 10000.0));
            } else {
                out->r -= (int16)rint(r / (totalContrib / 10000.0));
                out->g -= (int16)rint(g / (totalContrib / 10000.0));
                out->b -= (int16)rint(b / (totalContrib / 10000.0));
            }

            out++;

        }
    }

    return;
}

/** The Burt & Adelson Reduce operation.
 *  Allocates a new LPPixel array one quarter the size of in.
 */
LPPixel *reduce(LPPixel *in, uint32 w, uint32 h, bool wraparound) {
    uint32 outW = w >> 1;
    uint32 outH = h >> 1;

    LPPixel *out = (LPPixel*)calloc(outW * outH, sizeof(LPPixel));
    if (out == NULL) {
        cerr << endl
             << "enblend: out of memory (in reduce for out)" << endl;
        exit(1);
    }

    LPPixel *outIndex = out;
    for (uint32 outY = 0; outY < outH; outY++) {
        for (uint32 outX = 0; outX < outW; outX++) {

            double r = 0.0;
            double g = 0.0;
            double b = 0.0;
            uint32 noContrib = 10000;

            for (int32 m = 0; m < 5; m++) {
                int32 inX = 2 * (int32)outX + m - 2;

                if (wraparound) {
                    // Boundary condition: wrap around the image.
                    if (inX < 0) inX = (int32)w + inX;
                    if (inX >= (int32)w) inX = inX - (int32)w;
                } else {
                    // Boundary condition: replicate first and last column.
                    if (inX >= (int32)w) inX = w - 1;
                    if (inX < 0) inX = 0;
                }

                //if (w == (ROILastX - ROIFirstX + 1) && (outX == 0 || outX == outW - 1)) {
                //    cout << "r outX=" << outX << " inX=" << inX
                //         << " m=" << m << " outW=" << outW << endl;
                //}

                for (int32 n = 0; n < 5; n++) {
                    int32 inY = 2 * (int32)outY + n - 2;

                    // Boundary condition: replicate first and last column.
                    if (inY >= (int32)h) inY = h - 1;
                    if (inY < 0) inY = 0;

                    LPPixel *inPixel = &(in[inY * w + inX]);

                    if (inPixel->a != 255) {
                        // Transparent pixels don't count.
                        noContrib -= W100[m] * W100[n];
                    } else {
                        r += W[m] * W[n] * inPixel->r;
                        g += W[m] * W[n] * inPixel->g;
                        b += W[m] * W[n] * inPixel->b;
                    }
                }
            }

            // Adjust filter for any ignored transparent pixels.
            r = (noContrib == 0) ? 0.0 : r / (noContrib / 10000.0);
            g = (noContrib == 0) ? 0.0 : g / (noContrib / 10000.0);
            b = (noContrib == 0) ? 0.0 : b / (noContrib / 10000.0);

            outIndex->r = (int16)rint(r);
            outIndex->g = (int16)rint(g);
            outIndex->b = (int16)rint(b);
            outIndex->a = (noContrib == 0) ? 0 : 255;

            outIndex++;

        }
    }

    return out;
}

/** Create a Gaussian pyramid for the mask in the specified file.
 *  maskFile is ubbWidth * ubbHeight * MaskPixel.
 *  Saves the result in a temporary file.
 */
FILE *gaussianPyramidFile(FILE *maskFile, uint32 levels) {

    vector<LPPixel*> v;

    // Only consider the region-of-interest within image.
    uint32 ubbWidth = UBBLastX - UBBFirstX + 1;
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // Allocate space for level 0.
    LPPixel *g = (LPPixel*)malloc(roiWidth * roiHeight * sizeof(LPPixel));
    if (g == NULL) {
        cerr << endl
             << "enblend: out of memory (in gaussianPyramidFile for g)" << endl;
        exit(1);
    }

    // Allocate space for a scanline.
    MaskPixel *line = (MaskPixel*)malloc(roiWidth * sizeof(MaskPixel));
    if (line == NULL) {
        cerr << endl
             << "enblend: out of memory (in gaussianPyramidFile for line)" << endl;
        exit(1);
    }

    if (Verbose > 0) {
        cout << "Generating Gaussian pyramid: g0";
        cout.flush();
    }

    // Make sure we're at the beginning of the file.
    rewind(maskFile);

    // Copy image region-of-interest verbatim into g.
    LPPixel *gIndex = g;
    for (uint32 y = ROIFirstY - UBBFirstY; y <= ROILastY - UBBFirstY; y++) {

        fseek(maskFile, (y * ubbWidth + ROIFirstX - UBBFirstX) * sizeof(MaskPixel), SEEK_SET);
        readFromTmpfile(line, sizeof(MaskPixel), roiWidth, maskFile);

        MaskPixel *pixel = line;
        for (uint32 x = 0; x < roiWidth; x++) {
            gIndex->r = pixel->r;
            gIndex->g = pixel->g;
            gIndex->b = pixel->b;
            gIndex->a = pixel->a;
            gIndex++;
            pixel++;
        }
    }

    free(line);

    v.push_back(g);

    // If the user specified a 360 degree panorama, and the transition line
    // is near the left or right edge, then enable the wraparound boundary
    // condition.
    // If the transition line is near the edge, then the roiWidth is set to
    // the full output width in roiBounds.
    bool wraparound = Wraparound && (roiWidth == OutputWidth);

    // Make remaining levels.
    uint32 l = 1;
    while (l < levels) {
        if (Verbose > 0) {
            cout << " g" << v.size();
            cout.flush();
        }

        g = reduce(g, roiWidth >> (l-1), roiHeight >> (l-1), wraparound);
        v.push_back(g);

        l++;
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return dumpPyramidToTmpfile(v);
}

/** Create a Gaussian pyramid for the image in the specified file.
 *  uint32File is OutputWidth * OutputHeight * uint32.
 */
vector<LPPixel*> *gaussianPyramid(FILE *uint32File, uint32 levels) {

    vector<LPPixel*> *v = new vector<LPPixel*>();

    // Only consider the region-of-interest within image.
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // Allocate space for level 0.
    LPPixel *g = (LPPixel*)malloc(roiWidth * roiHeight * sizeof(LPPixel));
    if (g == NULL) {
        cerr << endl
             << "enblend: out of memory (in gaussianPyramid for g)" << endl;
        exit(1);
    }

    // Allocate space for a scanline.
    uint32 *line = (uint32*)malloc(roiWidth * sizeof(uint32));
    if (line == NULL) {
        cerr << endl
             << "enblend: out of memory (in gaussianPyramid for line)" << endl;
        exit(1);
    }

    if (Verbose > 0) {
        cout << "Generating Gaussian pyramid: g0";
        cout.flush();
    }

    // Copy image region-of-interest verbatim into g.
    LPPixel *gIndex = g;
    for (uint32 y = ROIFirstY; y <= ROILastY; y++) {

        fseek(uint32File, (y * OutputWidth + ROIFirstX) * sizeof(uint32), SEEK_SET);
        readFromTmpfile(line, sizeof(uint32), roiWidth, uint32File);

        uint32 *pixel = line;
        for (uint32 x = 0; x < roiWidth; x++) {
            gIndex->r = TIFFGetR(*pixel);
            gIndex->g = TIFFGetG(*pixel);
            gIndex->b = TIFFGetB(*pixel);
            gIndex->a = TIFFGetA(*pixel);
            gIndex++;
            pixel++;
        }
    }

    free(line);

    v->push_back(g);

    // If the user specified a 360 degree panorama, and the transition line
    // is near the left or right edge, then enable the wraparound boundary
    // condition.
    // If the transition line is near the edge, then the roiWidth is set to
    // the full output width in roiBounds.
    bool wraparound = Wraparound && (roiWidth == OutputWidth);

    // Make remaining levels.
    uint32 l = 1;
    while (l < levels) {
        if (Verbose > 0) {
            cout << " g" << v->size();
            cout.flush();
        }

        g = reduce(g, roiWidth >> (l-1), roiHeight >> (l-1), wraparound);
        v->push_back(g);

        l++;
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return v;
}

/** Create a Laplacian pyramid from imageFile with the specified number of
 *  levels.
 *  Returns a vector of pyramid levels.
 */
vector<LPPixel*> *laplacianPyramid(FILE *imageFile, uint32 levels) {
    // Only consider the region-of-interest within image.
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    // First create a Gaussian pyramid.
    vector<LPPixel*> *gp = gaussianPyramid(imageFile, levels);

    if (Verbose > 0) {
        cout << "Generating Laplacian pyramid:";
        cout.flush();
    }

    // If the user specified a 360 degree panorama, and the transition line
    // is near the left or right edge, then enable the wraparound boundary
    // condition.
    // If the transition line is near the edge, then the roiWidth is set to
    // the full output width in roiBounds.
    bool wraparound = Wraparound && (roiWidth == OutputWidth);

    // For each level, subtract the expansion of the next level.
    // Stop if there is no next level.
    for (uint32 l = 0; l < (levels-1); l++) {
        if (Verbose > 0) {
            cout << " l" << l;
            cout.flush();
        }
        expand((*gp)[l + 1], roiWidth >> (l+1), roiHeight >> (l+1),
                (*gp)[l], roiWidth >> l, roiHeight >> l,
                false,
                wraparound);
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return gp;
}

/** Create a Laplacian pyramid from imageFile with the specified number of
 *  levels.
 *  Returns the pyramid in a temp file.
 */
FILE *laplacianPyramidFile(FILE *imageFile, uint32 levels) {
    vector<LPPixel*> *v = laplacianPyramid(imageFile, levels);
    FILE *pyramidFile = dumpPyramidToTmpfile(*v);
    delete v;
    return pyramidFile;
}

/** Collapse the Laplacian pyramid given in p.
 *  Copy the result into the region-of-interest of dest.
 *  Use mask to set full transparency on pixels within the region-of-interest
 *  but outside the union of the images that make up the region-of-interest.
 */
void collapsePyramid(vector<LPPixel*> &p) {

    // Only operate within the region-of-interest.
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;

    if (Verbose > 0) {
        cout << "Collapsing Laplacian pyramid:";
        cout.flush();
    }

    // If the user specified a 360 degree panorama, and the transition line
    // is near the left or right edge, then enable the wraparound boundary
    // condition.
    // If the transition line is near the edge, then the roiWidth is set to
    // the full output width in roiBounds.
    bool wraparound = Wraparound && (roiWidth == OutputWidth);

    // For each level, add the expansion of the next level.
    // Work backwards from the smallest level to the largest.
    for (int l = (p.size()-2); l >= 0; l--) {
        if (Verbose > 0) {
            cout << " l" << l;
            cout.flush();
        }
        expand(p[l + 1], roiWidth >> (l+1), roiHeight >> (l+1),
                p[l], roiWidth >> l, roiHeight >> l,
                true,
                wraparound);
    }

    if (Verbose > 0) {
        cout << endl;
    }

    return;
}

void savePyramid(vector<LPPixel*> &p, char *prefix) {
    uint32 *image = (uint32*)calloc(OutputWidth * OutputHeight, sizeof(uint32));

    vector<LPPixel*> pCopy;
    uint32 roiWidth = ROILastX - ROIFirstX + 1;
    uint32 roiHeight = ROILastY - ROIFirstY + 1;
    for (unsigned int i = 0; i < p.size(); i++) {
        LPPixel *level = (LPPixel*)calloc((roiWidth >> i) * (roiHeight >> i), sizeof(LPPixel));
        pCopy.push_back(level);
    }

    // If the user specified a 360 degree panorama, and the transition line
    // is near the left or right edge, then enable the wraparound boundary
    // condition.
    // If the transition line is near the edge, then the roiWidth is set to
    // the full output width in roiBounds.
    bool wraparound = Wraparound && (roiWidth == OutputWidth);

    for (unsigned int i = 0; i < p.size(); i++) {
        char buf[512];
        snprintf(buf, 512, "%s%u.tif", prefix, i);
        cout << buf << endl;

        TIFF *outputTIFF = TIFFOpen(buf, "w");
        TIFFSetField(outputTIFF, TIFFTAG_ORIENTATION, 1);
        TIFFSetField(outputTIFF, TIFFTAG_SAMPLESPERPIXEL, 4);
        TIFFSetField(outputTIFF, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(outputTIFF, TIFFTAG_IMAGEWIDTH, OutputWidth);
        TIFFSetField(outputTIFF, TIFFTAG_IMAGELENGTH, OutputHeight);
        TIFFSetField(outputTIFF, TIFFTAG_PHOTOMETRIC, Photometric);
        TIFFSetField(outputTIFF, TIFFTAG_PLANARCONFIG, PlanarConfig);
 
        // Clear layers up to i
        for (unsigned int j = 0; j < i; j++) {
            memset(pCopy[j], 0, (roiWidth >> j) * (roiHeight >> j) * sizeof(LPPixel));
        }
        // Copy layer i from p
        memcpy(pCopy[i], p[i], (roiWidth >> i) * (roiHeight >> i) * sizeof(LPPixel));
        for (int j = i-1; j >= 0; j--) {
            expand(pCopy[j + 1], roiWidth >> (j+1), roiHeight >> (j+1),
                    pCopy[j], roiWidth >> j, roiHeight >> j,
                    true,
                    wraparound);
        }

        LPPixel *pixel = pCopy[0];
        for (uint32 y = ROIFirstY; y <= ROILastY; y++) {
            for (uint32 x = ROIFirstX; x <= ROILastX; x++) {
                pixel->r = min(255, max(0, (int)abs(pixel->r)));
                pixel->g = min(255, max(0, (int)abs(pixel->g)));
                pixel->b = min(255, max(0, (int)abs(pixel->b)));
                image[y * OutputWidth + x] =
                        (pixel->r & 0xFF)
                        | ((pixel->g & 0xFF) << 8)
                        | ((pixel->b & 0xFF) << 16)
                        | (0xFF << 24);
                pixel++;
            }
        }
        
        for (uint32 scan = 0; scan < OutputHeight; scan++) {
            TIFFWriteScanline(outputTIFF,
                    &(image[scan * OutputWidth]),
                    scan,
                    8);
        }

        TIFFClose(outputTIFF);

    }

    free(image);
    for (unsigned int i = 0; i < pCopy.size(); i++) {
        free(pCopy[i]);
    }

}
