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
#ifndef __ENBLEND_H__
#define __ENBLEND_H__

#include <tiffio.h>
#include <list>
#include <vector>

typedef struct {
    uint8 r;
    uint8 g;
    uint8 b;
    uint8 a;
} MaskPixel;

typedef struct {
    int16 r;
    int16 g;
    int16 b;
    int16 a;
} LPPixel;

// assemble.cc
FILE *assemble(std::list<char*> &filenames, bool pickOne);

// blend.cc
void blend(std::vector<LPPixel*> &whiteLP, FILE *blackLPFile, FILE *maskGPFile);

// bounds.cc
void ubbBounds(FILE *uint32File1, FILE *uint32File2);
uint32 roiBounds(FILE *maskFile);
void copyExcludedPixels(FILE *dst, FILE *src);
void copyROIToOutputWithMask(LPPixel *roi, FILE *uint32File, FILE *maskFile);

// io.cc
void closeTmpfile(FILE *f);
void readFromTmpfile(void *ptr, size_t size, size_t nmemb, FILE *stream);
void writeToTmpfile(void *ptr, size_t size, size_t nmemb, FILE *stream);
FILE *dumpToTmpfile(void *ptr, size_t size, size_t nmemb);
FILE *dumpPyramidToTmpfile(std::vector<LPPixel*> &v);
void saveMaskAsTIFF(MaskPixel *mask);

// mask.cc
FILE *createMask(FILE *whiteImageFile, FILE *blackImageFile);

// nearest.cc
void nearestFeatureTransform(MaskPixel *mask);

// pyramid.cc
uint32 filterHalfWidth(uint32 level, uint32 maxPixelValue);
FILE *gaussianPyramidFile(FILE *maskFile, uint32 levels);
std::vector<LPPixel*> *laplacianPyramid(FILE *imageFile, uint32 levels);
FILE *laplacianPyramidFile(FILE *imageFile, uint32 levels);
void collapsePyramid(std::vector<LPPixel*> &p);

#endif /* __ENBLEND_H__ */
