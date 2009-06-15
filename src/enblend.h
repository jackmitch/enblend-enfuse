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
#ifndef __ENBLEND_H__
#define __ENBLEND_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <list>
#include <stdio.h>

#include <boost/static_assert.hpp>

#include "common.h"
#include "numerictraits.h"
#include "fixmath.h"
#include "assemble.h"
#include "blend.h"
#include "bounds.h"
#include "mask.h"
#include "pyramid.h"

#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/transformimage.hxx"

using std::cout;
using std::endl;
using std::list;
using std::pair;

using vigra::BasicImage;
using vigra::CachedFileImage;
using vigra::CachedFileImageDirector;
using vigra::FindMinMax;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::initImage;
using vigra::initImageIf;
using vigra::inspectImage;
using vigra::NumericTraits;
using vigra::VigraFalseType;
using vigra::VigraTrueType;

namespace enblend {

/** Enblend's main blending loop. Templatized to handle different image types.
 */
template <typename ImagePixelType>
void enblendMain(list<ImageImportInfo*> &imageInfoList,
                 ImageExportInfo &outputImageInfo,
                 Rect2D &inputUnion) {

    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePixelComponentType ImagePixelComponentType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImageType ImageType;
    typedef typename EnblendNumericTraits<ImagePixelType>::AlphaPixelType AlphaPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::AlphaType AlphaType;
    typedef typename EnblendNumericTraits<ImagePixelType>::MaskPixelType MaskPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::MaskType MaskType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePyramidPixelType ImagePyramidPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::ImagePyramidType ImagePyramidType;
    typedef typename EnblendNumericTraits<ImagePixelType>::MaskPyramidPixelType MaskPyramidPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::MaskPyramidType MaskPyramidType;

    enum {ImagePyramidIntegerBits = EnblendNumericTraits<ImagePixelType>::ImagePyramidIntegerBits};
    enum {ImagePyramidFractionBits = EnblendNumericTraits<ImagePixelType>::ImagePyramidFractionBits};
    enum {MaskPyramidIntegerBits = EnblendNumericTraits<ImagePixelType>::MaskPyramidIntegerBits};
    enum {MaskPyramidFractionBits = EnblendNumericTraits<ImagePixelType>::MaskPyramidFractionBits};
    typedef typename EnblendNumericTraits<ImagePixelType>::SKIPSMImagePixelType SKIPSMImagePixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::SKIPSMAlphaPixelType SKIPSMAlphaPixelType;
    typedef typename EnblendNumericTraits<ImagePixelType>::SKIPSMMaskPixelType SKIPSMMaskPixelType;

    //sigset_t oldsigmask;
    //sigprocmask(SIG_BLOCK, &SigintMask, &oldsigmask);

    // Create the initial black image.
    Rect2D blackBB;
    pair<ImageType*, AlphaType*> blackPair =
        assemble<ImageType, AlphaType>(imageInfoList, inputUnion, blackBB);

    //sigprocmask(SIG_SETMASK, &oldsigmask, NULL);

    if (Checkpoint) {
        checkpoint(blackPair, outputImageInfo);
    }

    // mem usage before = 0
    // mem xsection = OneAtATime: inputUnion*imageValueType + inputUnion*AlphaValueType
    //                !OneAtATime: 2*inputUnion*imageValueType + 2*inputUnion*AlphaValueType
    // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType

    #ifdef ENBLEND_CACHE_IMAGES
    if (Verbose > VERBOSE_CFI_MESSAGES) {
        CachedFileImageDirector& v = CachedFileImageDirector::v();
        cerr << command
             << ": info: image cache statistics after loading black image\n";
        v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
        v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
        v.printStats(cerr, command + ": info: ");
        v.resetCacheMisses();
    }
    #endif

    // Main blending loop.
    while (!imageInfoList.empty()) {
        // Create the white image.
        Rect2D whiteBB;
        pair<ImageType*, AlphaType*> whitePair =
            assemble<ImageType, AlphaType>(imageInfoList, inputUnion, whiteBB);

        // mem usage before = inputUnion*ImageValueType + inputUnion*AlphaValueType
        // mem xsection = OneAtATime: inputUnion*imageValueType + inputUnion*AlphaValueType
        //                !OneAtATime: 2*inputUnion*imageValueType + 2*inputUnion*AlphaValueType
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 <<": info: image cache statistics after loading white image\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            v.printStats(cerr, command + ": info:     whiteImage", whitePair.first);
            v.printStats(cerr, command + ": info:     whiteAlpha", whitePair.second);
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif

        // Union bounding box of whiteImage and blackImage.
        Rect2D uBB = blackBB | whiteBB;

        if (Verbose > VERBOSE_UBB_MESSAGES) {
            cerr << command
                 << ": info: image union bounding box: "
                 << uBB
                 << endl;
        }

        // Intersection bounding box of whiteImage and blackImage.
        Rect2D iBB = blackBB & whiteBB;
        bool iBBValid = !iBB.isEmpty();

        if (Verbose > VERBOSE_IBB_MESSAGES) {
            cerr << command << ": info: image intersection bounding box: ";
            if (iBBValid) {
                cerr << iBB;
            } else {
                cerr << "(no intersection)";
            }
            cerr << endl;
        }

        // Determine what kind of overlap we have.
        Overlap overlap = inspectOverlap(uBB.apply(srcImageRange(*(blackPair.second))),
                                         uBB.apply(srcImage(*(whitePair.second))));

        // If white image is redundant, skip it and go to next images.
        if (overlap == CompleteOverlap) {
            // White image is redundant.
            delete whitePair.first;
            delete whitePair.second;
            cerr << command << ": warning: some images are redundant and will not be blended"
                 << endl;
            continue;
        }
        else if (overlap == NoOverlap && ExactLevels == 0) {
            // Images do not actually overlap.
            cerr << command << ": images do not overlap - they will be combined without blending\n"
                 << command << ": use the \"-l\" flag to force blending with a certain number of levels"
                 << endl;

            // Copy white image into black image verbatim.
            copyImageIf(srcImageRange(*(whitePair.first)),
                        maskImage(*(whitePair.second)),
                        destImage(*(blackPair.first)));
            copyImageIf(srcImageRange(*(whitePair.second)),
                        maskImage(*(whitePair.second)),
                        destImage(*(blackPair.second)));

            delete whitePair.first;
            delete whitePair.second;

            // Checkpoint results.
            if (Checkpoint) {
                if (Verbose > VERBOSE_CHECKPOINTING_MESSAGES) {
                    cerr << command << ": info: ";
                    if (imageInfoList.empty()) {
                        cerr << "writing final output" << endl;
                    } else {
                        cerr << "checkpointing" << endl;
                    }
                }
                checkpoint(blackPair, outputImageInfo);
            }

            blackBB = uBB;
            continue;
        }

        // Estimate memory requirements.
        if (Verbose > VERBOSE_MEMORY_ESTIMATION_MESSAGES) {
            long long bytes = 0;

            // Input images
            bytes += 2 * uBB.area() * sizeof(ImagePixelType);

            // Input alpha channels
            bytes += 2 * uBB.area() * sizeof(AlphaPixelType);

            // Mem used during mask generation:
            long long nftBytes = 0;
            if (!LoadMaskFileName.empty()) {
                nftBytes = 0;
            } else if (CoarseMask) {
                nftBytes =
                    2 * 1/8 * uBB.area() * sizeof(MaskPixelType)
                    + 2 * 1/8 * uBB.area() * sizeof(UInt32);
            } else {
                nftBytes =
                    2 * uBB.area() * sizeof(MaskPixelType)
                    + 2 * uBB.area() * sizeof(UInt32);
            }

            long long optBytes = 0;
            if (!LoadMaskFileName.empty()) {
                optBytes = 0;
            } else if (!OptimizeMask) {
                optBytes = 0;
            } else if (CoarseMask) {
                optBytes = 1/2 * iBB.area() * sizeof(UInt8);
            } else {
                optBytes = iBB.area() * sizeof(UInt8);
            }
            if (!VisualizeMaskFileName.empty()) optBytes *= 2;

            long long bytesDuringMask = bytes + std::max(nftBytes, optBytes);
            long long bytesAfterMask = bytes + uBB.area() * sizeof(MaskPixelType);

            bytes = std::max(bytesDuringMask, bytesAfterMask);

            cerr << command << ": info: estimated space required for mask generation: "
                 << static_cast<int>(ceil(bytes / 1000000.0))
                 << "MB" << endl;
        }

        // Create the blend mask.
        bool wraparoundForMask = Wraparound && (uBB.width() == inputUnion.width());

        MaskType* mask =
            createMask<ImageType, AlphaType, MaskType>(whitePair.first, blackPair.first,
                                                       whitePair.second, blackPair.second,
                                                       uBB, iBB, wraparoundForMask);

        // Calculate bounding box of seam line.
        Rect2D mBB;
        maskBounds(mask, uBB, mBB);

        if (!SaveMaskFileName.empty()) {
            ImageExportInfo maskInfo(SaveMaskFileName.c_str());
            maskInfo.setPosition(uBB.upperLeft());
            exportImage(srcImageRange(*mask), maskInfo);
        }

        // mem usage here = MaskType*ubb +
        //                  2*inputUnion*ImageValueType +
        //                  2*inputUnion*AlphaValueType

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after mask generation\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            v.printStats(cerr, command + ": info:     whiteImage", whitePair.first);
            v.printStats(cerr, command + ": info:     whiteAlpha", whitePair.second);
            v.printStats(cerr, command + ": info:     mask", mask);
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif

        // Calculate ROI bounds and number of levels from mBB.
        // ROI bounds must be at least mBB but not to extend uBB.
        Rect2D roiBB;
        unsigned int numLevels =
            roiBounds<ImagePixelComponentType>(inputUnion,
                                               iBB, mBB, uBB, roiBB,
                                               wraparoundForMask);
        bool wraparoundForBlend = Wraparound && (roiBB.width() == inputUnion.width());

        // Estimate memory requirements for this blend iteration
        if (Verbose > VERBOSE_MEMORY_ESTIMATION_MESSAGES) {
            // Maximum utilization is when all three pyramids have been built
            // mem xsection = 4 * roiBB.width() * SKIPSMImagePixelType
            //                + 4 * roiBB.width() * SKIPSMAlphaPixelType
            // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
            //      + (4/3)*roiBB*MaskPyramidType
            //      + 2*(4/3)*roiBB*ImagePyramidType
            long long bytes =
                inputUnion.area() * (sizeof(ImagePixelType) + 2 * sizeof(AlphaPixelType))
                + (4/3) * roiBB.area() * (sizeof(MaskPyramidPixelType)
                                          + 2 * sizeof(ImagePyramidPixelType))
                + (4 * roiBB.width()) * (sizeof(SKIPSMImagePixelType)
                                         + sizeof(SKIPSMAlphaPixelType));

            cerr << command << ": info: estimated space required for this blend step: "
                 << static_cast<int>(ceil(bytes / 1000000.0))
                 << "MB" << endl;
        }

        // Create a version of roiBB relative to uBB upperleft corner.
        // This is to access roi within images of size uBB.
        // For example, the mask.
        Rect2D roiBB_uBB = roiBB;
        roiBB_uBB.moveBy(-uBB.upperLeft());

        // Build Gaussian pyramid from mask.
        vector<MaskPyramidType*> *maskGP =
            gaussianPyramid<MaskType, MaskPyramidType,
            MaskPyramidIntegerBits, MaskPyramidFractionBits,
            SKIPSMMaskPixelType>(numLevels, wraparoundForBlend,
                                 roiBB_uBB.apply(srcImageRange(*mask)));
        //exportPyramid(maskGP, "mask");
        // mem usage before = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        // mem usage xsection = 3 * roiBB.width * MaskPyramidType
        // mem usage after = MaskType*ubb + 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType + (4/3)*roiBB*MaskPyramidType

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after calculating mask pyramid\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            v.printStats(cerr, command + ": info:     whiteImage", whitePair.first);
            v.printStats(cerr, command + ": info:     whiteAlpha", whitePair.second);
            v.printStats(cerr, command + ": info:     mask", mask);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats(cerr, command + ": info:     maskGP", i, (*maskGP)[i]);
            }
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif

        // Now it is safe to make changes to mask image.
        // Black out the ROI in the mask.
        // Make an roiBounds relative to uBB origin.
        initImage(roiBB_uBB.apply(destImageRange(*mask)),
                  NumericTraits<MaskPyramidPixelType>::zero());

        // Copy pixels inside whiteBB and inside white part of mask into black image.
        // These are pixels where the white image contributes outside of the ROI.
        // We cannot modify black image inside the ROI yet because we haven't built the
        // black pyramid.
        copyImageIf(uBB.apply(srcImageRange(*(whitePair.first))),
                    maskImage(*mask),
                    uBB.apply(destImage(*(blackPair.first))));

        // We no longer need the mask.
        delete mask;
        // mem usage after = 2*inputUnion*ImageValueType +
        //                   2*inputUnion*AlphaValueType +
        //                   (4/3)*roiBB*MaskPyramidType

        // Build Laplacian pyramid from white image.
        vector<ImagePyramidType*>* whiteLP =
            laplacianPyramid<ImageType, AlphaType, ImagePyramidType,
            ImagePyramidIntegerBits, ImagePyramidFractionBits,
            SKIPSMImagePixelType, SKIPSMAlphaPixelType>("whiteGP",
                                                        numLevels, wraparoundForBlend,
                                                        roiBB.apply(srcImageRange(*(whitePair.first))),
                                                        roiBB.apply(maskImage(*(whitePair.second))));

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after calculating white pyramid\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            v.printStats(cerr, command + ": info:     whiteImage", whitePair.first);
            v.printStats(cerr, command + ": info:     whiteAlpha", whitePair.second);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats(cerr, command + ": info:     maskGP", i, (*maskGP)[i]);
            }
            for (unsigned int i = 0; i < whiteLP->size(); i++) {
                v.printStats(cerr, command + ": info:     whiteLP", i, (*whiteLP)[i]);
            }
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif
        // mem usage after = 2*inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        //                   + (4/3)*roiBB*MaskPyramidType + (4/3)*roiBB*ImagePyramidType
        // mem xsection = 4 * roiBB.width() * SKIPSMImagePixelType
        //                + 4 * roiBB.width() * SKIPSMAlphaPixelType

        // We no longer need the white rgb data.
        delete whitePair.first;
        // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        //                   + (4/3)*roiBB*MaskPyramidType + (4/3)*roiBB*ImagePyramidType

        // Build Laplacian pyramid from black image.
        vector<ImagePyramidType*>* blackLP =
            laplacianPyramid<ImageType, AlphaType, ImagePyramidType,
            ImagePyramidIntegerBits, ImagePyramidFractionBits,
            SKIPSMImagePixelType, SKIPSMAlphaPixelType>("blackGP",
                                                        numLevels, wraparoundForBlend,
                                                        roiBB.apply(srcImageRange(*(blackPair.first))),
                                                        roiBB.apply(maskImage(*(blackPair.second))));

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after calculating black pyramid\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            v.printStats(cerr, command + ": info:     whiteAlpha", whitePair.second);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats(cerr, command + ": info:     maskGP", i, (*maskGP)[i]);
            }
            for (unsigned int i = 0; i < whiteLP->size(); i++) {
                v.printStats(cerr, command + ": info:     whiteLP", i, (*whiteLP)[i]);
            }
            for (unsigned int i = 0; i < blackLP->size(); i++) {
                v.printStats(cerr, command + ": info:     blackLP", i, (*blackLP)[i]);
            }
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif
        //exportPyramid(blackLP, "enblend_black_lp");

        // Peak memory xsection is here!
        // mem xsection = 4 * roiBB.width() * SKIPSMImagePixelType
        //                + 4 * roiBB.width() * SKIPSMAlphaPixelType
        // mem usage after = inputUnion*ImageValueType + 2*inputUnion*AlphaValueType
        //      + (4/3)*roiBB*MaskPyramidType
        //      + 2*(4/3)*roiBB*ImagePyramidType

        // Make the black image alpha equal to the union of the
        // white and black alpha channels.
        initImageIf(whiteBB.apply(destImageRange(*(blackPair.second))),
                    whiteBB.apply(maskImage(*(whitePair.second))),
                    NumericTraits<AlphaPixelType>::max());

        // We no longer need the white alpha data.
        delete whitePair.second;

        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType
        //      + (4/3)*roiBB*MaskPyramidType + 2*(4/3)*roiBB*ImagePyramidType

        // Blend pyramids
        ConvertScalarToPyramidFunctor<MaskPixelType, MaskPyramidPixelType, MaskPyramidIntegerBits, MaskPyramidFractionBits> whiteMask;
        blend(maskGP, whiteLP, blackLP, whiteMask(NumericTraits<MaskPixelType>::max()));
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after blending pyramids\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            for (unsigned int i = 0; i < maskGP->size(); i++) {
                v.printStats(cerr, command + ": info:     maskGP", i, (*maskGP)[i]);
            }
            for (unsigned int i = 0; i < whiteLP->size(); i++) {
                v.printStats(cerr, command + ": info:     whiteLP", i, (*whiteLP)[i]);
            }
            for (unsigned int i = 0; i < blackLP->size(); i++) {
                v.printStats(cerr, command + ": info:     blackLP", i, (*blackLP)[i]);
            }
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif

        // delete mask pyramid
        //exportPyramid(maskGP, "enblend_mask_gp");
        for (unsigned int i = 0; i < maskGP->size(); i++) {
            delete (*maskGP)[i];
        }
        delete maskGP;

        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + 2*(4/3)*roiBB*ImagePyramidType

        // delete white pyramid
        //exportPyramid(whiteLP, "enblend_white_lp");
        for (unsigned int i = 0; i < whiteLP->size(); i++) {
            delete (*whiteLP)[i];
        }
        delete whiteLP;

        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType + (4/3)*roiBB*ImagePyramidType

        //exportPyramid(blackLP, "enblend_blend_lp");
        // collapse black pyramid
        collapsePyramid<SKIPSMImagePixelType>(wraparoundForBlend, blackLP);
        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after collapsing black pyramid\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            for (unsigned int i = 0; i < blackLP->size(); i++) {
                v.printStats(cerr, command + ": info:     blackLP", i, (*blackLP)[i]);
            }
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif

        // copy collapsed black pyramid into black image ROI, using black alpha mask.
        copyFromPyramidImageIf<ImagePyramidType, MaskType, ImageType,
                               ImagePyramidIntegerBits, ImagePyramidFractionBits>(
                srcImageRange(*((*blackLP)[0])),
                roiBB.apply(maskImage(*(blackPair.second))),
                roiBB.apply(destImage(*(blackPair.first))));

        // delete black pyramid
        for (unsigned int i = 0; i < blackLP->size(); i++) {
            delete (*blackLP)[i];
        }
        delete blackLP;

        // mem usage after = inputUnion*ImageValueType + inputUnion*AlphaValueType

        // Checkpoint results.
        if (Checkpoint) {
            if (Verbose > VERBOSE_CHECKPOINTING_MESSAGES) {
                cerr << command << ": info: ";
                if (imageInfoList.empty()) {
                    cerr << "writing final output" << endl;
                } else {
                    cerr << "checkpointing" << endl;
                }
            }
            checkpoint(blackPair, outputImageInfo);
        }

        #ifdef ENBLEND_CACHE_IMAGES
        if (Verbose > VERBOSE_CFI_MESSAGES) {
            CachedFileImageDirector& v = CachedFileImageDirector::v();
            cerr << command
                 << ": info: image cache statistics after checkpointing\n";
            v.printStats(cerr, command + ": info:     blackImage", blackPair.first);
            v.printStats(cerr, command + ": info:     blackAlpha", blackPair.second);
            v.printStats(cerr, command + ": info: ");
            v.resetCacheMisses();
        }
        #endif

        // Now set blackBB to uBB.
        blackBB = uBB;

    }

    if (!Checkpoint) {
        if (Verbose > VERBOSE_CHECKPOINTING_MESSAGES) {
            cerr << command << ": info: writing final output" << endl;
        }
        checkpoint(blackPair, outputImageInfo);
    }

    delete blackPair.first;
    delete blackPair.second;
};

} // namespace enblend

#endif /* __ENBLEND_H__ */

// Local Variables:
// mode: c++
// End:
