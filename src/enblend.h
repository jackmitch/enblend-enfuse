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
 * along with Foobar; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef __ENBLEND_H__
#define __ENBLEND_H__

#include <tiffio.h>
#include <vector>

uint32 *createMask(uint32 *outputBuf, TIFF *inputTIFF);

void thinMask(uint32 *mask);

std::vector<uint32*> *gaussianPyramid(uint32 *mask, int32 levels);
std::vector<uint32*> *laplacianPyramid(uint32 *image, int32 levels);
std::vector<uint32*> *laplacianPyramid(TIFF *image, int32 levels);
void collapsePyramid(std::vector<uint32*> *p, uint32 *dest);

void blend(std::vector<uint32*> *mask,
        std::vector<uint32*> *inputLP,
        std::vector<uint32*> *outputLP);

#define WHITE 0xFFFFFFFF
#define BLACK 0xFF000000
#define BLUE  0xFFFF0000
#define GREEN 0xFF00FF00
#define RED   0xFF0000FF
#define TRANS 0x00000000

#endif /* __ENBLEND_H__ */
