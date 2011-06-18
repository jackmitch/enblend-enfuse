/*
 * Copyright (C) 2011 Mikolaj Leszczynski
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
#include <queue>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

#include "vigra/functorexpression.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdcachedfileimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/stdconvolution.hxx"
#include "vigra/bordertreatment.hxx"
#include "vigra/labelimage.hxx"
#include <vigra/transformimage.hxx>
#include "common.h"
#include "maskcommon.h"
#include "masktypedefs.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::priority_queue;

using vigra::NumericTraits;
using vigra::triple;
using vigra::stride;
using vigra::StridedImageIterator;
using vigra::kernel2d;
using vigra::BorderTreatmentMode;
using vigra::labelImage;

//using vigra::convolveImage;
#define BIT_MASK_DIR 0x03
#define BIT_MASK_OPDIR 0x02
#define BIT_MASK_OPEN 0x04
#define BIT_MASK_DIVIDE 0x0F
//#define GRAPHCUT_DBG

namespace enblend{
    
    
    template<class MaskImageIterator, class MaskAccessor, typename GraphmaskPixelType>
    void findExtremePoints(MaskImageIterator upperleft, MaskImageIterator lowerright, MaskAccessor ma, 
            MaskImageIterator* output){
        typedef typename IMAGETYPE<GraphmaskPixelType>::traverser IteratorType;
        typedef NumericTraits<GraphmaskPixelType> SrcPixelTraits;
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
            if(!iter2Done && ma(iter2) != SrcPixelTraits::max()){
                output[0] = iter2;
                break;
            }
            if(!iterDone && ma(iter) != SrcPixelTraits::max()){
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
            if(!iter2Done && ma(iter2) != SrcPixelTraits::max()){
                output[1] = iter2;
                break;
            }
            if(!iterDone && ma(iter) != SrcPixelTraits::max()){
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
            if(!iter2Done && ma(iter2) != SrcPixelTraits::max()){
                output[2] = iter2;
                break;
            }
            if(!iterDone && ma(iter) != SrcPixelTraits::max()){
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
            if(!iter2Done && ma(iter2) != SrcPixelTraits::max()){
                output[3] = iter2;
                break;
            }
            if(!iterDone && ma(iter) != SrcPixelTraits::max()){
                output[3] = iter;
                break;
            }
            iter.x++;
            iter2.y--;
        }
         
    }
    
    template <typename ImageType>
    struct CostComparer
    {
        public:
            CostComparer(ImageType* image):img(image){}
            bool operator()(Point2D& a, Point2D& b){
                return (*img)[a] > (*img)[b];
            }
        protected:
            ImageType* img;
    };
    
    
    struct pointHash{
        std::size_t operator()(vigra::Point2D const& p) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, p.x);
            boost::hash_combine(seed, p.y);
            return seed;
        }
    };
    
    struct PathEqualityFunctor
    {
        public:
            PathEqualityFunctor(boost::unordered_set<Point2D, pointHash>* a_, 
                    boost::unordered_set<Point2D, pointHash>* b_):left(a_), right(b_){}
            bool operator()(Diff2D a2, Diff2D b2){
                Point2D a(a2), b(b2);
                //add border to detect seams close to border
                a-=Point2D(1,1);b-=Point2D(1,1);
                if((left->find(a)) != left->end() && (right->find(b)) != right->end() ||
                   (right->find(a)) != right->end() && (left->find(b) != left->end()))
                        return false;
                else return true;
            }
        protected:
            boost::unordered_set<Point2D, pointHash>* left;
            boost::unordered_set<Point2D, pointHash>* right;
    };
    template<class MaskPixelType, class DestPixelType>
    struct FinalMaskFunctor
    {
        public:
            
            FinalMaskFunctor(int l, int w):leftMaskLabel(l),leftMaskColor(w){
            if(leftMaskColor == 255)
                rightMaskColor = 0;
            else rightMaskColor = 255;}
            
            DestPixelType operator()(MaskPixelType const& arg1, MaskPixelType const& arg2,  
                            MaskPixelType const& arg3) const{
                
                if(arg3 == 0){
                    if(arg1)
                        return leftMaskColor;
                    else return rightMaskColor;
                } 
                else {
                    if(arg3 == leftMaskLabel){
                        return leftMaskColor;
                    }
                    else return rightMaskColor;
                }
                
            }
        protected:
            int leftMaskLabel;
            int leftMaskColor;
            int rightMaskColor;
    };

    template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor,
            class MaskImageIterator, class MaskAccessor>
    inline void
    graphCut(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src1,
                            pair<SrcImageIterator, SrcAccessor> src2,
                            pair<DestImageIterator, DestAccessor> dest,
                            triple<MaskImageIterator, MaskImageIterator, MaskAccessor> mask1,
                            pair<MaskImageIterator, MaskAccessor> mask2,
                            nearest_neigbor_metric_t norm, boundary_t boundary, const Rect2D& iBB)
    {
        graphCut(src1.first, src1.second, src1.third,
                                src2.first, src2.second,
                                dest.first, dest.second,
                                mask1.first, mask1.second, mask1.third,
                                mask2.first, mask2.second,
                                norm, boundary, iBB);
    }
    
    template <class ImageType>
    uint heuristic_dummy(Point2D pt, Point2D destpt, ImageType* img){
        return 0;
    }
    
    template <class ImageType, class MaskPixelType>
    void getNeighbourList(Point2D src, ImageType* img, Point2D* list, Diff2D bounds){
        MaskPixelType maskVal = NumericTraits<MaskPixelType>::max();
        //return neighbour points from top to left in clockwise order
        if(src.y == 1 || (*img)[src.y-2][src.x] == maskVal)
            list[0] = Point2D(-1,-1);
        else
            list[0] = src(0,-2);
        
        if(src.x == bounds.x - 1|| (*img)[src.y][src.x+2] == maskVal)
            list[1] = Point2D(-1,-1);
        else
            list[1] = src(2,0);
        
        if(src.y == bounds.y - 1|| (*img)[src.y+2][src.x] == maskVal)
            list[2] = Point2D(-1,-1);
        else
            list[2] = src(0,2);
        
        if(src.x == 1 || (*img)[src.y][src.x-2] == maskVal)
            list[3] = Point2D(-1,-1);
        else
            list[3] = src(-2,0);
        
    }
    
    template <class ImageType>
    vector<Point2D>* path(Point2D pt, Point2D srcpt, ImageType* img){
        vector<Point2D>* vec = new vector<Point2D>;
        Point2D current = pt;
        vec->push_back(pt);
        do{
            switch((*img)[current(1,1)] & BIT_MASK_DIR){
                case 0:
                    vec->push_back(current(0,-2));
                    break;
                case 1:
                    vec->push_back(current(2,0));
                    break;
                case 2:
                    vec->push_back(current(0,2));
                    break;
                case 3:
                    vec->push_back(current(-2,0));
                    break;
                default:
                    #ifdef GRAPHCUT_DBG
                        cout<<"path tracing error"<<endl;
                    #endif
                    break;
            }
            current = vec->back();
        }while(current != srcpt);
            
        return vec;
    }
    
    template <class ImageType>
    uint getEdgeWeight(int dir, Point2D pt, ImageType* img){
        switch(dir){
            case 0:
                return (*img)[pt(1,0)];
                break;
            case 1:
                return (*img)[pt(2,1)];
                break;
            case 2:
                return (*img)[pt(1,2)];
                break;
            case 3:
                return (*img)[pt(0,1)];
                break;
        }
    }
    
    template <class ImageType, class MaskPixelType>
    vector<Point2D>* A_star(Point2D srcpt, Point2D destpt, ImageType* img, Diff2D bounds){
        MaskPixelType zeroVal = NumericTraits<MaskPixelType>::zero();
        typedef priority_queue<Point2D, vector<Point2D>, CostComparer<ImageType> > Queue;
        Queue *openset = new Queue(CostComparer<ImageType>(img));
        //uint *heur(Point2D, Point2D, ImageType*) = &heuristic_dummy;
        long score;
        long count = 0;
        bool scoreIsBetter;
        Point2D list[4], current, neighbour;
        openset->push(srcpt);
        
        while(!openset->empty()){
            current = openset->top();
            openset->pop();
            count++;
            if(current == destpt){
                cout<<"Graphcut completed after visiting "<<count<<" nodes"<<endl;
                delete openset;
                return path<ImageType>(destpt, srcpt, img);
            }
            
            getNeighbourList<ImageType, MaskPixelType>(current, img, list, bounds);
            for(int i = 0; i < 4; i++){
                score = 0;
                scoreIsBetter = false;
                neighbour = list[i];
                if(neighbour != Point2D(-1,-1) && (*img)[neighbour(1,1)] == zeroVal){
                    score = (*img)[neighbour] + getEdgeWeight(i, current, img);
                    if(((*img)[neighbour(1,1)] & BIT_MASK_OPEN) == 0){
                        openset->push(neighbour);
                        (*img)[neighbour(1,1)] += BIT_MASK_OPEN;
                        scoreIsBetter = true;
                    }
                    else if(score < (*img)[neighbour])
                        scoreIsBetter = true;
                    
                    if(scoreIsBetter){
                        (*img)[neighbour(1,1)] &= BIT_MASK_OPEN;
                        (*img)[neighbour(1,1)] += i;
                        (*img)[neighbour(1,1)] ^= BIT_MASK_OPDIR;
                        (*img)[neighbour] = score;
                    }
                    
                }
            }
        }
        cout<<"Graphcut failed after visiting "<<count<<" nodes"<<endl;
        delete openset;
        return new vector<Point2D>();
    }
   
    
    Point2D convertFromDual(const Point2D& dualPixel){
        int stride = 2;
        return (Point2D((dualPixel->x - 1)/stride, (dualPixel->y - 1)/stride));
    }
    
    
    template<class ImageType>
    void dividePath(ImageType* img, vector<Point2D>* cut,
                        boost::unordered_set<Point2D, pointHash>* left,
                        boost::unordered_set<Point2D, pointHash>* right
                        ){
        Point2D previous, current;
        std::vector<Point2D>::iterator previousDual;
        for (std::vector<Point2D>::iterator currentDual = cut->begin();
                                            currentDual != cut->end();
                                            ++currentDual) {
            
            current = convertFromDual(*currentDual);
            
            if(currentDual == cut->begin()){
                left->insert(current(0,1));
                right->insert(current(1,1));
                
                switch((*img)[(*currentDual)(1,1)] & BIT_MASK_DIR){
                case 0:
                    left->insert(current);
                    right->insert(current(1,0));
                    break;
                case 1:
                    left->insert(current);
                    left->insert(current(1,0));
                    break;
                case 3:
                    right->insert(current);
                    right->insert(current(1,0));
                    break;
                case 2:
                default:
                    #ifdef GRAPHCUT_DBG
                        cout<<"path dividing error"<<endl;
                    #endif
                    break;
                }
            } else if(*currentDual == cut->back()){
            
                if(current.y == 0){
                    left->insert(current);
                    right->insert(current(1, 0));
                    left->insert(current(0,-1));
                    right->insert(current(1, -1));
                }
                else if(current.x == 0){
                    left->insert(current(0,1));
                    right->insert(current);
                    left->insert(current(-1,1));
                    right->insert(current(-1,0));
                }
                else{
                    left->insert(current(1,0));
                    right->insert(current(1,1));
                    left->insert(current(2,0));
                    right->insert(current(2,1));
                }
                    
            } else{
                
                int a = ((*img)[(*previousDual)(1,1)] & BIT_MASK_DIR) << 2,
                    b = ((*img)[(*currentDual)(1,1)] & BIT_MASK_DIR);
                switch(a+b){
                case 0: // top -> top
                    left->insert(current);
                    right->insert(current(1,0));
                    left->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 1:// top -> right
                    left->insert(current);
                    left->insert(current(1,0));
                    left->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 3:// top -> left
                    right->insert(current);
                    right->insert(current(1,0));
                    left->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 4://right->top
                    left->insert(current);
                    right->insert(current(1,0));
                    right->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 5://right->right
                    left->insert(current);
                    left->insert(current(1,0));
                    right->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 6://right->bottom
                    left->insert(current);
                    left->insert(current(1,0));
                    right->insert(current(0,1));
                    left->insert(current(1,1));
                    break;
                case 9://bottom->right
                    right->insert(current);
                    left->insert(current(1,0));
                    right->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 10://bottom->bottom
                    right->insert(current);
                    left->insert(current(1,0));
                    right->insert(current(0,1));
                    left->insert(current(1,1));
                    break;
                case 11://bottom->left
                    right->insert(current);
                    left->insert(current(1,0));
                    left->insert(current(0,1));
                    left->insert(current(1,1));
                    break;
                case 12://left->top
                    left->insert(current);
                    right->insert(current(1,0));
                    left->insert(current(0,1));
                    left->insert(current(1,1));
                    break;
                case 14://left->bottom
                    right->insert(current);
                    right->insert(current(1,0));
                    right->insert(current(0,1));
                    left->insert(current(1,1));
                    break;
                case 15://left->left
                    right->insert(current);
                    right->insert(current(1,0));
                    left->insert(current(0,1));
                    left->insert(current(1,1));
                    break;
                case 2:// top -> bottom
                case 7://right->left
                case 8://bottom->top
                case 13://left->right
                default:
                    #ifdef GRAPHCUT_DBG
                        cout<<"path dividing error"<<endl;
                    #endif
                    break;
                }
            }
            
            previous = current;
            previousDual = currentDual;
        }
    }
    
    /* Graph-cut
     * implementation of an algorithm based on an article by:
     * V.Kwatra, A.Schoedl, I.Essa, G.Turk, A.Bobick
     * "Graphcut Textures: Image and Video Synthesis Using Graph Cuts"
     * 
     * The algorithm finds the best seam for two images using a graph-based approach.
    */
    
    template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor,
            class MaskImageIterator, class MaskAccessor>
    void
    graphCut(SrcImageIterator src1_upperleft, SrcImageIterator src1_lowerright, SrcAccessor sa1,
                            SrcImageIterator src2_upperleft, SrcAccessor sa2,
                            DestImageIterator dest_upperleft, DestAccessor da,
                            MaskImageIterator mask1_upperleft, MaskImageIterator mask1_lowerright, MaskAccessor ma1,
                            MaskImageIterator mask2_upperleft, MaskAccessor ma2,
                            nearest_neigbor_metric_t norm, boundary_t boundary, const Rect2D& iBB)
    {
        typedef typename SrcAccessor::value_type SrcPixelType;
        typedef vigra::NumericTraits<SrcPixelType> SrcPixelTraits;
        typedef typename SrcPixelTraits::Promote SrcPromoteType;

        typedef typename DestAccessor::value_type DestPixelType;
        typedef vigra::NumericTraits<DestPixelType> DestPixelTraits;
        
        typedef typename MaskAccessor::value_type MaskPixelType;
        typedef vigra::NumericTraits<MaskPixelType> MaskPixelTraits;
        
        typedef UInt8 BasePixelType;
        typedef vigra::NumericTraits<BasePixelType> BasePixelTraits;
        typedef typename BasePixelTraits::Promote BasePromotePixelType;
        typedef vigra::NumericTraits<BasePromotePixelType> BasePromotePixelTraits;
        typedef typename BasePromotePixelTraits::Promote GraphPixelType;

        const SrcPixelType background = SrcPixelTraits::zero();
        const Diff2D size(src1_lowerright.x - src1_upperleft.x,
                          src1_lowerright.y - src1_upperleft.y),
                          masksize(mask1_lowerright.x - mask1_upperleft.x,
                          mask1_lowerright.y - mask1_upperleft.y);

        IMAGETYPE<BasePixelType> tmp(size), mask(size+Diff2D(2,2));
        IMAGETYPE<BasePromotePixelType> graphtmp(size + size + Diff2D(1,1));
        IMAGETYPE<GraphPixelType> graph(size + size + Diff2D(1,1));
        IMAGETYPE<MaskPixelType> testmask(masksize);
        IMAGETYPE<BasePixelType>::traverser pts[4];
        vector<Point2D>* dualPath;
        boost::unordered_set<Point2D, pointHash> a, b;
        
        
        const Diff2D graphsize(graph.lowerRight().x - graph.upperLeft().x,
                          graph.lowerRight().y - graph.upperLeft().y);
        const Rect2D gBB(Point2D(1,1), Size2D(graphsize)),
                     bBB(Point2D(1,1), Size2D(size));
        
        
        
        //difference image calculation
        combineTwoImagesMP(src1_upperleft, src1_lowerright, sa1, 
                src2_upperleft, sa2,
                tmp.upperLeft(), tmp.accessor(),
                PixelDifferenceFunctor<SrcPixelType, BasePixelType>()
                );     
        
        //masking overlap region borders
        combineThreeImagesMP(mask1_upperleft + iBB.upperLeft(), mask1_upperleft + iBB.lowerRight(), ma1, 
                mask2_upperleft + iBB.upperLeft(), ma2,
                tmp.upperLeft(), tmp.accessor(),
                tmp.upperLeft(), tmp.accessor(),
                ifThenElse(!(Arg1() & Arg2()), Param(BasePixelTraits::max()), Arg3())
                );
        
        //look for possible start and end points     
        findExtremePoints<IMAGETYPE<BasePixelType>::traverser,
        IMAGETYPE<BasePixelType>::Accessor, BasePixelType>
        (tmp.upperLeft(), tmp.lowerRight(), tmp.accessor(), pts);
        
        //copying to a grid
        copyImage(srcImageRange(tmp), stride(2,2,gBB.apply(destImage(graphtmp))));
                     
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(tmp), ImageExportInfo("diff.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(graphtmp), ImageExportInfo("diff2.tif").setPixelType("UINT8"));
        exportImage(mask1_upperleft, mask1_lowerright, ma1, ImageExportInfo("mask1.tif").setPixelType("UINT8"));
        exportImage(mask2_upperleft, mask1_lowerright, ma2, ImageExportInfo("mask2.tif").setPixelType("UINT8"));
#endif
        //calculating differences between pixels that are adjacent in the original image
        vigra::Kernel2D<BasePromotePixelType> edgeWeightKernel;
        edgeWeightKernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
                         0, 1, 0,
                         1, 1, 1,
                         0, 1, 0;
        edgeWeightKernel.setBorderTreatment(vigra::BORDER_TREATMENT_CLIP);
        convolveImage(srcImageRange(graphtmp), destImage(graph), kernel2d(edgeWeightKernel));
        
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(graph), ImageExportInfo("graph.tif").setPixelType("UINT8"));
#endif
        //find optimal cut in dual graph
        dualPath = A_star<IMAGETYPE<GraphPixelType>, BasePixelType>(
            Point2D((pts[0] - tmp.upperLeft())*2) + Diff2D(1,1), 
            Point2D((pts[2] - tmp.upperLeft())*2) + Diff2D(1,1),
            &graph, graphsize - Diff2D(1,1));
        
        dividePath<IMAGETYPE<GraphPixelType> >(&graph, dualPath, &a, &b);
        
        //labels areas to belong to left/right images
        //adds a 1-pixel border to catch any area that is cut off by the seam
        labelImage(srcIterRange(Diff2D(), Diff2D() + size + Diff2D(2,2)), 
                        destImage(mask), 
                        false, PathEqualityFunctor(&a, &b));
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(mask), ImageExportInfo("labels.tif").setPixelType("UINT8"));
#endif
        copyImage(bBB.apply(srcImageRange(mask)), destIter(testmask.upperLeft() += iBB.upperLeft()));
        
        if(ma1(mask1_upperleft) == 255)
        {
                combineThreeImagesMP(mask1_upperleft, mask1_lowerright, ma1,
                             mask2_upperleft, ma2,
                             testmask.upperLeft(), testmask.accessor(),
                             dest_upperleft, da,
                             FinalMaskFunctor<MaskPixelType, DestPixelType>(2, 255));
        } else{
                combineThreeImagesMP(mask2_upperleft, mask1_lowerright, ma2,
                             mask1_upperleft, ma1,
                             testmask.upperLeft(), testmask.accessor(),
                             dest_upperleft, da,
                             FinalMaskFunctor<MaskPixelType, DestPixelType>(2, 0));
        }
    }
    

}/* namespace enblend */
#endif	/* GRAPHCUT_H */

