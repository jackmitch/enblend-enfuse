/*
 * Copyright (C) 2004-2005 Andrew Mihal
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

// Defines to control how many -v flags are required for each type
// of message to be produced on stdout.
#define VERBOSE_ASSEMBLE_MESSAGES           0
#define VERBOSE_ABB_MESSAGES                1
#define VERBOSE_UBB_MESSAGES                1
#define VERBOSE_IBB_MESSAGES                1
#define VERBOSE_BLEND_MESSAGES              0
#define VERBOSE_NUMLEVELS_MESSAGES          0
#define VERBOSE_ROIBB_SIZE_MESSAGES         1
#define VERBOSE_MEMORY_ESTIMATION_MESSAGES  1
#define VERBOSE_CHECKPOINTING_MESSAGES      0
#define VERBOSE_INPUT_IMAGE_INFO_MESSAGES   1
#define VERBOSE_INPUT_UNION_SIZE_MESSAGES   1
#define VERBOSE_COLOR_CONVERSION_MESSAGES   0
#define VERBOSE_NFT_MESSAGES                0
#define VERBOSE_MASK_MESSAGES               0
#define VERBOSE_PYRAMID_MESSAGES            0
#define VERBOSE_CFI_MESSAGES                2

#define GDA_KMAX    32

namespace enblend {

/** The different image overlap classifications. */
enum Overlap {NoOverlap, PartialOverlap, CompleteOverlap};

} // namespace enblend

#endif /* __COMMON_H__ */
