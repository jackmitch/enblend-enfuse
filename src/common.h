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
#ifndef __COMMON_H__
#define __COMMON_H__

#include <tiffio.h>
#include <list>
#include <vector>

#include "vigra_ext/ROI.h"

// Defines to control how many -v flags are required for each type
// of message to be produced on stdout.
#define VERBOSE_ASSEMBLE_MESSAGES           0
#define VERBOSE_BLEND_MESSAGES              0
#define VERBOSE_NUMLEVELS_MESSAGES          0
#define VERBOSE_ROIBB_SIZE_MESSAGES         0
#define VERBOSE_MEMORY_ESTIMATION_MESSAGES  0
#define VERBOSE_CHECKPOINTING_MESSAGES      0
#define VERBOSE_INPUT_IMAGE_INFO_MESSAGES   0
#define VERBOSE_INPUT_UNION_SIZE_MESSAGES   0
#define VERBOSE_COLOR_CONVERSION_MESSAGES   0
#define VERBOSE_NFT_MESSAGES                0
#define VERBOSE_PYRAMID_MESSAGES            0

namespace enblend {

//typedef struct {
//    uint8 r;
//    uint8 g;
//    uint8 b;
//    uint8 a;
//} MaskPixel;
//
//typedef struct {
//    int16 r;
//    int16 g;
//    int16 b;
//    int16 a;
//} LPPixel;

//template<typename T1> T1 GetMaxAlpha();

typedef vigra_ext::ROI<vigra::Diff2D> EnblendROI;

enum Overlap {NoOverlap, PartialOverlap, CompleteOverlap};

//// assemble.cc
//FILE *assemble(std::list<char*> &filenames, bool pickOne);
//
//// blend.cc
//void blend(std::vector<LPPixel*> &whiteLP, FILE *blackLPFile, FILE *maskGPFile);
//
//// bounds.cc
//void ubbBounds(FILE *uint32File1, FILE *uint32File2);
//uint32 roiBounds(FILE *maskFile);
//void copyExcludedPixels(FILE *dst, FILE *src, FILE *mask);
//void copyROIToOutputWithMask(LPPixel *roi, FILE *uint32File, FILE *maskFile);
//
//// io.cc
//void closeTmpfile(FILE *f);
//void readFromTmpfile(void *ptr, size_t size, size_t nmemb, FILE *stream);
//void writeToTmpfile(void *ptr, size_t size, size_t nmemb, FILE *stream);
//FILE *dumpToTmpfile(void *ptr, size_t size, size_t nmemb);
//FILE *dumpPyramidToTmpfile(std::vector<LPPixel*> &v);
//void saveMaskAsTIFF(MaskPixel *mask, char *filename);
//void savePyramidAsTIFF(std::vector<LPPixel*> &p, char *filename);
//
//// mask.cc
//FILE *createMask(FILE *whiteImageFile, FILE *blackImageFile);
//
//// nearest.cc
//void nearestFeatureTransform(MaskPixel *mask);
//
//// pyramid.cc
//uint32 filterHalfWidth(uint32 level, uint32 maxPixelValue);
//FILE *gaussianPyramidFile(FILE *maskFile, uint32 levels);
//std::vector<LPPixel*> *laplacianPyramid(FILE *imageFile, uint32 levels);
//FILE *laplacianPyramidFile(FILE *imageFile, uint32 levels);
//void collapsePyramid(std::vector<LPPixel*> &p);
//
//// Macros for accessing 8-bit color fields in 32-bit words on different machines
//#ifdef WORDS_BIGENDIAN
//#define GetR(rgba)      (((rgba) >> 24) & 0xff)
//#define GetG(rgba)      (((rgba) >> 16) & 0xff)
//#define GetB(rgba)      (((rgba) >> 8) & 0xff)
//#define GetA(rgba)      ((rgba) & 0xff)
//#define PACK(a, b, g, r)    (((a) & 0xff) | (((b) & 0xff) << 8) | (((g) & 0xff) << 16) | (((r) & 0xff) << 24))
//#else
//#define GetA(abgr)      (((abgr) >> 24) & 0xff)
//#define GetB(abgr)      (((abgr) >> 16) & 0xff)
//#define GetG(abgr)      (((abgr) >> 8) & 0xff)
//#define GetR(abgr)      ((abgr) & 0xff)
//#define PACK(a, b, g, r)    (((r) & 0xff) | (((g) & 0xff) << 8) | (((b) & 0xff) << 16) | (((a) & 0xff) << 24))
//#endif

} // namespace enblend

#endif /* __COMMON_H__ */
