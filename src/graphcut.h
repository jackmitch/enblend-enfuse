/* 
 * File:   graphCut.h
 * Author: rosomack
 *
 * Created on May 30, 2011, 3:23 PM
 */

#ifndef GRAPHCUT_H
#define	GRAPHCUT_H

#include <iostream>
#ifdef _WIN32
#include <cmath>
#else
#include <math.h>
#endif
#include <stdlib.h>
#include <utility>

#include "vigra/functorexpression.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdcachedfileimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/stdconvolution.hxx"
#include "vigra/bordertreatment.hxx"
#include "common.h"
#include "maskcommon.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;

using vigra::NumericTraits;
using vigra::triple;
using vigra::stride;
using vigra::StridedImageIterator;
using vigra::kernel2d;
using vigra::BorderTreatmentMode;

//using vigra::convolveImage;

namespace enblend{
    
    
    template<class MaskImageIterator, class MaskAccessor, typename GraphmaskPixelType>
    void findExtremePoints(MaskImageIterator upperleft, MaskImageIterator lowerright, MaskAccessor ma, 
            MaskImageIterator* output){
        typedef typename IMAGETYPE<GraphmaskPixelType>::traverser IteratorType;
        IteratorType iter = upperleft;
        IteratorType iter2 = upperleft;
        bool iterDone = false, iter2Done = false;
        Diff2D dim(lowerright.x - upperleft.x,
                          lowerright.y - upperleft.y);
        dim -= Diff2D(1,1);
        //look for upper left corner of overlap region
        while(!iterDone && !iter2Done){
            
            if(iter.x == lowerright.x)
                iterDone = true;
            if(iter2.y == lowerright.y)
                iter2Done = true;
            if(!iter2Done && ma(iter2) == VISUALIZE_RGB_COLOR_GRAY255){
                output[0] = iter2;
                break;
            }
            if(!iterDone && ma(iter) == VISUALIZE_RGB_COLOR_GRAY255){
                output[0] = iter;
                break;
            }
            iter.x++;
            iter2.y++;
        }
         
        //look for upper right corner of overlap region    
        iter.x = upperleft.x;
        iter.x+=dim.x;
        iter.y = upperleft.y;
        iter2 = iter;  
        iterDone = false; iter2Done = false;
        
        while(!iterDone && !iter2Done){
            
            if(iter.x < upperleft.x)
                iterDone = true;
            if(iter2.y == lowerright.y)
                iter2Done = true;
            if(!iter2Done && ma(iter2) == VISUALIZE_RGB_COLOR_GRAY255){
                output[1] = iter2;
                break;
            }
            if(!iterDone && ma(iter) == VISUALIZE_RGB_COLOR_GRAY255){
                output[1] = iter;
                break;
            }
            iter.x--;
            iter2.y++;
        }
            
            
        //look for lower right corner of overlap region   
        iter.x = upperleft.x; 
        iter.x += dim.x;
        iter.y = upperleft.y;
        iter.y += dim.y;
        iter2 = iter;  
        iterDone = false; iter2Done = false;
        
        while(!iterDone && !iter2Done){
            
            if(iter.x < upperleft.x)
                iterDone = true;
            if(iter2.y < upperleft.y)
                iter2Done = true;
            if(!iter2Done && ma(iter2) == VISUALIZE_RGB_COLOR_GRAY255){
                output[2] = iter2;
                break;
            }
            if(!iterDone && ma(iter) == VISUALIZE_RGB_COLOR_GRAY255){
                output[2] = iter;
                break;
            }
            iter.x--;
            iter2.y--;
        }
            
        //look for lower left corner of overlap region    
        iter.x = upperleft.x;
        iter.y = upperleft.y;
        iter.y += dim.y;
        iter2 = iter;  
        iterDone = false; iter2Done = false;
        
        while(!iterDone && !iter2Done){
            
            if(iter.x == lowerright.x)
                iterDone = true;
            if(iter2.y < upperleft.y)
                iter2Done = true;
            if(!iter2Done && ma(iter2) == VISUALIZE_RGB_COLOR_GRAY255){
                output[3] = iter2;
                break;
            }
            if(!iterDone && ma(iter) == VISUALIZE_RGB_COLOR_GRAY255){
                output[3] = iter;
                break;
            }
            iter.x++;
            iter2.y--;
        }
            
            
            
    }

    template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor,
            class MaskImageIterator, class MaskAccessor>
    inline void
    graphCut(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src1,
                            pair<SrcImageIterator, SrcAccessor> src2,
                            pair<DestImageIterator, DestAccessor> dest,
                            triple<MaskImageIterator, MaskImageIterator, MaskAccessor> mask1,
                            pair<MaskImageIterator, MaskAccessor> mask2,
                            nearest_neigbor_metric_t norm, boundary_t boundary)
    {
        graphCut(src1.first, src1.second, src1.third,
                                src2.first, src2.second,
                                dest.first, dest.second,
                                mask1.first, mask1.second, mask1.third,
                                mask2.first, mask2.second,
                                norm, boundary);
    }
    
    template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor,
            class MaskImageIterator, class MaskAccessor>
    void
    graphCut(SrcImageIterator src1_upperleft, SrcImageIterator src1_lowerright, SrcAccessor sa1,
                            SrcImageIterator src2_upperleft, SrcAccessor sa2,
                            DestImageIterator dest_upperleft, DestAccessor da,
                            MaskImageIterator mask1_upperleft, MaskImageIterator mask1_lowerright, MaskAccessor ma1,
                            MaskImageIterator mask2_upperleft, MaskAccessor ma2,
                            nearest_neigbor_metric_t norm, boundary_t boundary)
    {
        typedef typename SrcAccessor::value_type SrcPixelType;
        typedef vigra::NumericTraits<SrcPixelType> SrcPixelTraits;
        typedef typename SrcPixelTraits::Promote SrcPromoteType;

        typedef typename DestAccessor::value_type DestPixelType;
        typedef vigra::NumericTraits<DestPixelType> DestPixelTraits;
        
        typedef typename SrcAccessor::value_type MaskPixelType;
        typedef vigra::NumericTraits<SrcPixelType> MaskPixelTraits;
        
        typedef UInt8 GraphPixelType;
        typedef RGBValue<UInt8> GraphmaskPixelType;

        const SrcPixelType background = SrcPixelTraits::zero();
        const Diff2D size(src1_lowerright.x - src1_upperleft.x,
                          src1_lowerright.y - src1_upperleft.y);

        IMAGETYPE<GraphPixelType> tmp(size), graphtmp(size + size - Diff2D(1,1));
        IMAGETYPE<GraphPixelType> graph(size + size - Diff2D(1,1));
        IMAGETYPE<GraphmaskPixelType> graphmask(size);
        
        
        //difference image calculation
        combineTwoImagesMP(src1_upperleft, src1_lowerright, sa1, 
                src2_upperleft, sa2,
                tmp.upperLeft(), tmp.accessor(),
                PixelDifferenceFunctor<SrcPixelType, GraphPixelType>()
                );        

        //copying to a grid
        copyImage(srcImageRange(tmp), stride(2,2,destImage(graphtmp)));
        //mask image calculation (for storing visited pixels etc.)
        combineTwoImagesMP(mask1_upperleft, mask1_lowerright, ma1, 
                mask2_upperleft, ma2,
                graphmask.upperLeft(), graphmask.accessor(),
                ifThenElse(Arg1() & Arg2(), Param(VISUALIZE_RGB_COLOR_GRAY255), Param(VISUALIZE_RGB_COLOR_GRAY0))
                );
        
        //DEBUG
        exportImage(srcImageRange(tmp), ImageExportInfo("gc1.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(graphmask), ImageExportInfo("mask.tif"));
        exportImage(srcImageRange(graphtmp), ImageExportInfo("gc2.tif").setPixelType("UINT8"));
        //exportImage(mask1_upperleft, mask1_lowerright, ma1, ImageExportInfo("src1.tif").setPixelType("UINT8"));
        //exportImage(src2_upperleft, sa2, ImageExportInfo("src2.tif").setPixelType("UINT8"));
        //
        
        //calculating differences between pixels that are adjacent in the original image
        vigra::Kernel2D<GraphPixelType> edgeWeightKernel;
        edgeWeightKernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
                         0, 1, 0,
                         1, 1, 1,
                         0, 1, 0;
        edgeWeightKernel.setBorderTreatment(vigra::BORDER_TREATMENT_CLIP);//BORDER_TREATMENT_CLIP);
        convolveImage(srcImageRange(graphtmp), destImage(graph), kernel2d(edgeWeightKernel));
        
        //DEBUG
        exportImage(srcImageRange(graph), ImageExportInfo("gc3.tif").setPixelType("UINT8"));
        //
        
        //look for possible start and end points
        IMAGETYPE<GraphmaskPixelType>::traverser pts[4];
        findExtremePoints<IMAGETYPE<GraphmaskPixelType>::traverser,
        IMAGETYPE<GraphmaskPixelType>::ConstAccessor, GraphmaskPixelType>
        (graphmask.upperLeft(),graphmask.lowerRight(), graphmask.accessor(), pts);
        
        
        
    }
    

}/* namespace enblend */
#endif	/* GRAPHCUT_H */

