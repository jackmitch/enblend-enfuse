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
#include "vigra/transformimage.hxx"
#include "vigra/contourcirculator.hxx"
#include <vigra/mathutil.hxx>
#include "common.h"
#include "maskcommon.h"
#include "masktypedefs.h"
#include "nearest.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::priority_queue;

using vigra::NumericTraits;
using vigra::triple;
using vigra::stride;
using vigra::StridedImageIterator;
using vigra::Kernel2D;
using vigra::Kernel1D;
using vigra::kernel2d;
using vigra::kernel1d;
using vigra::BorderTreatmentMode;
using vigra::labelImage;
using namespace vigra::functor;

#define BIT_MASK_DIR 0x03
#define BIT_MASK_OPDIR 0x02
#define BIT_MASK_OPEN 0x04
#define BIT_MASK_DIVIDE 0x0F
//#define GRAPHCUT_DBG

namespace enblend{
    
    struct pointHash{
        std::size_t operator()(vigra::Point2D const& p) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, p.x);
            boost::hash_combine(seed, p.y);
            return seed;
        }
    };
    
    class CheckpointPixels{
        public:
            CheckpointPixels(){}
            boost::unordered_set<Point2D, pointHash> top, bottom;
            ~CheckpointPixels(){
                this->clear();
            }
            
            void clear(){
                top.clear();
                bottom.clear();
            }
            
    };
    
    int distab(const Point2D& a, const Point2D& b){
        return std::abs((a.x-b.x))+std::abs((a.y-b.y));
    }
    
    template<class MaskImageIterator, class MaskAccessor, class MaskPixelType,
                        class DestImageIterator, class DestAccessor>
    vector<Point2D>* findExtremePoints(MaskImageIterator mask1_upperleft, MaskImageIterator mask1_lowerright, MaskAccessor ma1,
                            MaskImageIterator mask2_upperleft, MaskAccessor ma2, 
                            DestImageIterator dest_upperleft, DestAccessor da,
                            nearest_neigbor_metric_t& norm, boundary_t& boundary,
                            const Rect2D& iBB){
        
        typedef vigra::NumericTraits<MaskPixelType> MaskPixelTraits;
        typedef typename IMAGETYPE<MaskPixelType>::traverser IteratorType;
        typedef vigra::CrackContourCirculator<IteratorType> Circulator;
		typedef triple<IteratorType, IteratorType, uint> EPointContainer;
        IMAGETYPE<MaskPixelType> nfttmp(mask1_lowerright - mask1_upperleft + Diff2D(2,2), MaskPixelTraits::max()/2),
                                nft(iBB.lowerRight() - iBB.upperLeft() + Diff2D(2,2), MaskPixelTraits::max()/2),
                                overlap(iBB.lowerRight() - iBB.upperLeft() + Diff2D(2,2));
        IteratorType nftiter = nft.upperLeft(), overlapiter = overlap.upperLeft(), previous;
        Circulator *circ, *end; bool offTheBorder = false, inOverlap = false, ready = false;
        vector<Point2D> *list = new vector<Point2D>();
		vector<EPointContainer> entryPointList;
		EPointContainer *max = NULL;
		pair<IteratorType, IteratorType> tmp;
        Point2D tmpPoint;
		uint counter = 0;
        
        nearestFeatureTransform(srcIterRange(mask1_upperleft, mask1_lowerright, ma1),
                            srcIter(mask2_upperleft, ma2),
                            destIter(nfttmp.upperLeft() + Diff2D(1,1)),
                            norm, boundary);
        
        copyImage(srcIterRange(nfttmp.upperLeft() + Diff2D(1,1), nfttmp.lowerRight() - Diff2D(1,1)),
                  destIter(dest_upperleft, da));
	
		copyImage(iBB.apply(srcIterRange(nfttmp.upperLeft()+Diff2D(1,1), nfttmp.lowerRight()-Diff2D(1,1))),
                  destIter(nft.upperLeft() + Diff2D(1,1)));
    
        combineTwoImagesMP(iBB.apply(srcIterRange(mask1_upperleft, mask1_lowerright, ma1)),
                            iBB.apply(srcIter(mask2_upperleft, ma2)),
                            destIter(overlap.upperLeft() + Diff2D(1,1)),
                            ifThenElse(Arg1()&&Arg2(), Param(MaskPixelTraits::max()), Param(MaskPixelTraits::zero())));
        
        #ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(nfttmp), ImageExportInfo("./debug/nfttotal.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(nft), ImageExportInfo("./debug/nft.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(overlap), ImageExportInfo("./debug/overlap.tif").setPixelType("UINT8"));
#endif
        
        circ = new Circulator(nftiter+Diff2D(1,0), vigra::FourNeighborCode::South);
        end = new Circulator(nftiter+Diff2D(1,0), vigra::FourNeighborCode::South);
        
        previous = circ->outerPixel();
        
        do{
            (*circ)++;
            if(nft.accessor()(previous) != nft.accessor()(circ->outerPixel())){
				if(ready){
					entryPointList.push_back(EPointContainer(tmp.first, tmp.second, counter));
					tmp.first = previous;
					tmp.second = circ->outerPixel();
				}
				else{
					ready = true;
					tmp = pair<IteratorType, IteratorType>(previous, circ->outerPixel());
				}
				counter = 0;
            }
			++counter;
            previous = circ->outerPixel();
        } while(*circ != *end);
        
        if(entryPointList.empty())
			return list;
        
        for(typename vector<EPointContainer>::iterator i = entryPointList.begin(); i != entryPointList.end(); ++i){
			if(max == NULL){
				max = &(*i);
				continue;
			}
			if(i->third > max->third)
				max = &(*i);
		}
		
		if(max != NULL){
			delete circ;
			delete end;
			Diff2D dir = max->second - max->first;
			if(dir.x != 0){
				if(max->second.x < max->first.x)
					max->second = max->first;
				circ = new Circulator(max->second, vigra::FourNeighborCode::West);
				end = new Circulator(max->second, vigra::FourNeighborCode::West);
			}
			else if(dir.y != 0){
				if(max->second.y < max->first.y)
					max->second = max->first;
				circ = new Circulator(max->second, vigra::FourNeighborCode::North);
				end = new Circulator(max->second, vigra::FourNeighborCode::North);
			}    
		}
		else return list;
        
        while(circ != end){
            (*circ)++;
            tmpPoint = circ->outerPixel() - nft.upperLeft();
            if(!offTheBorder && !(nft.accessor()(circ->outerPixel()) == MaskPixelTraits::max()/2)){
                offTheBorder = true;
                Point2D tmpPoint2 = Point2D((tmpPoint/*- iBB.upperLeft()*/)*2 - Diff2D(1,1));
                if((tmpPoint.x /*- iBB.upperLeft().y*/)>= 0 && 
                        (tmpPoint.y /*- iBB.upperLeft().y*/)>= 0 &&
                        (tmpPoint.x /*- iBB.upperLeft().x*/)< iBB.lowerRight().x && 
                        (tmpPoint.y /*- iBB.upperLeft().y*/)< iBB.lowerRight().y)
                        list->push_back(tmpPoint2);
                #ifdef GRAPHCUT_DBG
                    cout<<"Start point: "<<tmpPoint2<<endl;
                    #endif
            }
            
            if(offTheBorder){
                
                if(!inOverlap && overlap[tmpPoint] == MaskPixelTraits::max()){
                    inOverlap = true;
                    Point2D tmpPoint2 = Point2D((tmpPoint/*- iBB.upperLeft()*/)*2 - Diff2D(1,1));
                    if(!list->empty() && *(list->begin()) != tmpPoint2 && tmpPoint.x >= 0 && tmpPoint.y >= 0)
                        list->push_back(tmpPoint2);

                    #ifdef GRAPHCUT_DBG
                    cout<<"Entering overlap: "<<tmpPoint2<<endl;
                    #endif
                }
                if(inOverlap && overlap[tmpPoint] == MaskPixelTraits::zero() &&
                        nft.accessor()(circ->outerPixel()) != MaskPixelTraits::max()/2){
                    inOverlap = false;
                    Point2D tmpPoint2 = Point2D((tmpPoint/*- iBB.upperLeft()*/)*2 - Diff2D(1,1));
                    if(tmpPoint.x >= 0 && tmpPoint.y >= 0)
                        list->push_back(tmpPoint2);

                    #ifdef GRAPHCUT_DBG
                    cout<<"Leaving overlap: "<<tmpPoint2<<endl;
                    #endif
                }

                if(nft.accessor()(circ->outerPixel()) == MaskPixelTraits::max()/2){
                    Point2D tmpPoint2 = Point2D((previous - nft.upperLeft()/*- iBB.upperLeft()- Diff2D(0,1)*/)*2 - Diff2D(1,1));
                    if(!list->empty() && *(list->rbegin()) != tmpPoint2 && tmpPoint.x >= 0 && tmpPoint.y >= 0)
                        list->push_back(tmpPoint2);

                    #ifdef GRAPHCUT_DBG
                    cout<<"Endpoint reached: "<<tmpPoint2<<endl;
                    #endif
                    break;
                }
            }
            previous = circ->outerPixel();
        }

        delete circ;
        delete end;
        return list;
    }
    
    template <typename ImageType>
    class CostComparer
    {
        public:
            CostComparer(const ImageType* image):img(image){}
            bool operator()(const Point2D& a, const Point2D& b) const {
                if(a == Point2D(-20,-20))
                    return totalScore >  (*img)[b];
                else if(b == Point2D(-20,-20))
                    return (*img)[a] > totalScore;
                return (*img)[a] > (*img)[b];
            }
            void setTotalScore(long s){totalScore = s;}
        protected:
            const ImageType* img;
            long totalScore;
    };
    
    struct PathEqualityFunctor
    {
        public:
            PathEqualityFunctor(boost::unordered_set<Point2D, pointHash>* a_, 
                    boost::unordered_set<Point2D, pointHash>* b_,
                    Point2D offset_):left(a_), right(b_), offset(offset_){}
            bool operator()(Diff2D a2, Diff2D b2){
                Point2D a(a2), b(b2);
                //add border to detect seams close to border
                a-=Point2D(1,1);b-=Point2D(1,1);
                //a-= offset; b-= offset;
                if(((left->find(a)) != left->end() && (right->find(b)) != right->end()) ||
                   ((right->find(a)) != right->end() && (left->find(b) != left->end())))
                        return false;
                else return true;
            }
        protected:
            boost::unordered_set<Point2D, pointHash>* left;
            boost::unordered_set<Point2D, pointHash>* right;
            Point2D offset;
    };
    template<class MaskPixelType, class DestPixelType>
    struct FinalMaskFunctor
    {
        public:
            
            FinalMaskFunctor(int l, int w):leftMaskLabel(l),leftMaskColor(w){
                if(leftMaskColor == 255)
                    rightMaskColor = 0;
                else rightMaskColor = 255;
            }
            
            DestPixelType operator()(MaskPixelType const& arg1) const{
                
                    if(arg1 == leftMaskLabel)
                        return leftMaskColor;
                    else return rightMaskColor;
                
            }
        protected:
            int leftMaskLabel;
            int leftMaskColor;
            int rightMaskColor;
    };
    template<class MaskPixelType>
    struct CountFunctor
    {
        public:
            
            CountFunctor(int *c, int *c2):color(c), count(c2){}
            
            MaskPixelType operator()(MaskPixelType const& arg1, MaskPixelType const& arg2) const{
                
                    if(arg1 > 0 && arg1 == arg2)
                        (*color)++;
                    if(arg1 == 255)
                        (*count)++;
               return arg1; 
            }
            int *color, *count;
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
    void getNeighbourList(Point2D src, ImageType* img, Point2D* list, 
                Diff2D bounds, CheckpointPixels* bp){
        //MaskPixelType maskVal = NumericTraits<MaskPixelType>::max();
        //return neighbour points from top to left in clockwise order
        bool check = false;
        if(src.y == 1 /*|| (*img)[src.y-2][src.x] == maskVal*/){
            list[0] = Point2D(-1,-1);
            check = true;
        }
        else
            list[0] = src(0,-2);

        if(src.x == bounds.x - 1/*|| (*img)[src.y][src.x+2] == maskVal*/){
            list[1] = Point2D(-1,-1);
            check = true;
        }
        else
            list[1] = src(2,0);

        if(src.y == bounds.y - 1/*|| (*img)[src.y+2][src.x] == maskVal*/){
            list[2] = Point2D(-1,-1);
            check = true;
        }
        else
            list[2] = src(0,2);

        if(src.x == 1 /*|| (*img)[src.y][src.x-2] == maskVal*/){
            list[3] = Point2D(-1,-1);
            check = true;
        }
        else
            list[3] = src(-2,0);
        
        if(bp->bottom.find(src) != bp->bottom.end()){
            list[0] = Point2D(-20,-20);
            list[1] = Point2D(-20,-20);
            list[2] = Point2D(-20,-20);
            list[3] = Point2D(-20,-20);
        }
        
        if(check){
            
		if(bp->bottom.find(src) != bp->bottom.end()
			||bp->bottom.find(src(1,0)) != bp->bottom.end()
			||bp->bottom.find(src(1,1)) != bp->bottom.end()
			||bp->bottom.find(src(0,1)) != bp->bottom.end()){
			if(list[1] == Point2D(-1,-1))
			list[1] = Point2D(-20,-20);
			else if(list[2] == Point2D(-1,-1))
			list[2] = Point2D(-20,-20);
			else if(list[3] == Point2D(-1,-1))
			list[3] = Point2D(-20,-20);
		}
            
        }
    }
    
    template <class ImageType, class MaskPixelType>
    void getNeighbourList(CheckpointPixels* bp, boost::unordered_set<Point2D, pointHash>::iterator* auxList1,
        boost::unordered_set<Point2D, pointHash>::iterator* auxList2){
        
		*auxList1 = bp->top.begin();
		*auxList2 = bp->top.end();
	
    }
    
    template <class ImageType>
    vector<Point2D>* path(Point2D pt, Point2D srcpt, ImageType* img, CheckpointPixels* bp){
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
        }while(bp->top.find(current) == bp->top.end());
        return vec;
    }
    
    template <class ImageType>
    uint getEdgeWeight(int dir, Point2D pt, ImageType* img, bool endpt, Diff2D bounds){
        if(!endpt){
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
        else{
            if(pt.y == 1)
                return (*img)[pt(1,0)];
            else if(pt.x == bounds.x - 2)
                return (*img)[pt(2,1)];
            else if(pt.y == bounds.y - 2)
                return (*img)[pt(1,2)];
            else if(pt.x == 1)
                return (*img)[pt(0,1)];
        }
        return 0;
    }
    
    template <class ImageType, class GradientImageType, class MaskPixelType>
    vector<Point2D>* A_star(Point2D srcpt, Point2D destpt, ImageType* img
                ,GradientImageType* gradientX, GradientImageType* gradientY, Diff2D bounds,
                CheckpointPixels* bp){
        MaskPixelType zeroVal = NumericTraits<MaskPixelType>::zero();
        typedef priority_queue<Point2D, vector<Point2D>, CostComparer<ImageType> > Queue;
        CostComparer<ImageType> costcomp(img);
        Queue *openset = new Queue(costcomp);
        //uint *heur(Point2D, Point2D, ImageType*) = &heuristic_dummy;
        long score = 0, totalScore = 0;
        long count = 0;
        int grada, gradb;
        bool scoreIsBetter, pushflag, destOpen = false;
        Point2D list[4], current, neighbour, destNeighbour;
        boost::unordered_set<Point2D, pointHash>::iterator auxListBegin, auxListEnd;
        openset->push(srcpt);
        
        while(!openset->empty()){
            current = openset->top();
            openset->pop();
            count++;
            if(current == destpt){
			#ifdef GRAPHCUT_DBG
                cout<<"Graphcut completed after visiting "<<count<<" nodes"<<endl;
			#endif
                delete openset;
                return path<ImageType>(destNeighbour, srcpt, img, bp);
            }
            if(current == Point2D(-10,-10))
                getNeighbourList<ImageType, MaskPixelType>(bp, &auxListBegin, &auxListEnd);
            else getNeighbourList<ImageType, MaskPixelType>(current, img, list, bounds, bp);
            if(current != Point2D(-10,-10)){
                for(int i = 0; i < 4; i++){
                    score = 0;
                    scoreIsBetter = false;
                    pushflag = false;
                    neighbour = list[i];
                    if(neighbour == Point2D(-20,-20) || (neighbour != Point2D(-1,-1) && 
                            neighbour != Point2D(-10,-10) && (*img)[neighbour(1,1)] == zeroVal)){
                        if(neighbour == Point2D(-20,-20))
                            score = (*img)[current] + getEdgeWeight(i, current, img, true, bounds);
                        else{ 
                            if(i%2 == 0){
                                grada = std::abs((*gradientY)[(current)/ 2 ]);
                                gradb = std::abs((*gradientY)[(neighbour)/ 2 ]);
                            } else{
                                grada = std::abs((*gradientX)[(current)/ 2 ]);
                                gradb = std::abs((*gradientX)[(neighbour)/ 2 ]);
                            }
                            if((grada + gradb) > 0)
                                score = (*img)[current] + getEdgeWeight(i, current, img, false, bounds) *
                                        (grada + gradb);
                            else score = (*img)[current] + getEdgeWeight(i, current, img, false, bounds);
                        }
                        if(neighbour == Point2D(-20,-20) && !destOpen){
                            pushflag = true;
                            destOpen = true;
                            scoreIsBetter = true;
                        }
                        else if(neighbour != Point2D(-20,-20) && ((*img)[neighbour(1,1)] & BIT_MASK_OPEN) == 0){
                            pushflag = true;
                            (*img)[neighbour(1,1)] += BIT_MASK_OPEN;
                            scoreIsBetter = true;
                        }
                        else if((neighbour == Point2D(-20,-20) && score < totalScore) || (neighbour != Point2D(-20,-20) && score < (*img)[neighbour]))
                            scoreIsBetter = true;

                        if(scoreIsBetter){
                            if(neighbour == Point2D(-20,-20)){
                                totalScore = score;
                                destNeighbour = current;
                                costcomp.setTotalScore(totalScore);
                            } else{
                                (*img)[neighbour(1,1)] &= BIT_MASK_OPEN;
                                (*img)[neighbour(1,1)] += i;
                                (*img)[neighbour(1,1)] ^= BIT_MASK_OPDIR;
                                (*img)[neighbour] = score;
                            }
                            if(pushflag)
                                    openset->push(neighbour);
                        }

                    }
                }
            } else {
                for(boost::unordered_set<Point2D, pointHash>::iterator x = auxListBegin; x != auxListEnd; ++x){
                    score = 0;
                    scoreIsBetter = false;
                    pushflag = false;
                    neighbour = *x;
                    if(neighbour != Point2D(-1,-1) && (*img)[neighbour(1,1)] == zeroVal){
                        score = 0;    
                        if(((*img)[neighbour(1,1)] & BIT_MASK_OPEN) == 0){
                            pushflag = true;
                            (*img)[neighbour(1,1)] += BIT_MASK_OPEN;
                            scoreIsBetter = true;
                        }
                        else if(score < (*img)[neighbour])
                            scoreIsBetter = true;

                        if(scoreIsBetter){
                            (*img)[neighbour(1,1)] &= BIT_MASK_OPEN;
                            (*img)[neighbour] = score;
                            if(pushflag)
                                    openset->push(neighbour);
                        }

                    }
                }
            }
        }
        #ifdef GRAPHCUT_DBG
			cout<<"Graphcut failed after visiting "<<count<<" nodes"<<endl;
		#endif
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
                        boost::unordered_set<Point2D, pointHash>* right,
                        Diff2D bounds, const Rect2D& iBB){
        Point2D previous, current, next, tmp;
        ushort closest = 0;
        Diff2D dim(iBB.lowerRight() - iBB.upperLeft());
        std::vector<Point2D>::iterator previousDual, nextDual;
        for (std::vector<Point2D>::iterator currentDual = cut->begin(),
                                            nextDual = cut->begin();
                                            currentDual != cut->end();
                                            ++currentDual) {
            
            current = convertFromDual(*currentDual);
            next = convertFromDual(*nextDual);
            
            if(currentDual == cut->begin()){
                ++nextDual;
                next = convertFromDual(*nextDual);
                //find closest edge
                int dist;
                Diff2D toLowerRight = dim - current;
                closest = 0;
                dist = current.y;
                if(toLowerRight.x < dist){
                    closest = 1;
                    dist = toLowerRight.x;
                }
                if(toLowerRight.y < dist){
                    closest = 2;
                    dist = toLowerRight.y;
                }
                if(current.x < dist){
                    closest = 3;
                }
                
                switch(closest){
                case 0:
                    right->insert(current(0,0));
                    left->insert(current(1,0));
                    break;
                case 1:
                    left->insert(current(1,1));
                    right->insert(current(1,0));
                    break;
                case 2:
                    left->insert(current(0,1));
                    right->insert(current(1,1));
                    break;
                case 3:
                    right->insert(current(0,1));
                    left->insert(current(0,0));
                    break;
                default:
                    #ifdef GRAPHCUT_DBG
                        cout<<"path dividing error"<<endl;
                    #endif
                    break;
                }
                
                int a;
                if(next.x == current.x){
                    if(next.y < current.y)
                        a = 0;
                    else a = 2;
                } else{
                    if(next.x > current.x)
                        a = 1;
                    else a = 3;
                }
                
                switch(closest){
                    case 0:
                        switch(a){
                        case 2:
                            right->insert(current(0,1));
                            left->insert(current(1,1));
                            break;
                        case 1:
                            right->insert(current(0,1));
                            right->insert(current(1,1));
                            break;
                        case 3:
                            left->insert(current(0,1));
                            left->insert(current(1,1));
                            break;
                        case 0:
				left->insert(current(0,0));
                            right->insert(current(1,0));
                            break;
                        default:
                            #ifdef GRAPHCUT_DBG
                                cout<<"path dividing error"<<endl;
                            #endif
                            break;
                        }
                        break;
                        
                    case 1:
                        switch(a){
                        case 0:
                            left->insert(current(0,0));
                            left->insert(current(0,1));
                            break;
                        case 2:
                            right->insert(current(0,0));
                            right->insert(current(0,1));
                            break;
                        case 3:
                            right->insert(current(0,0));
                            left->insert(current(0,1));
                            break;
                        case 1:
				right->insert(current(1,1));
                            left->insert(current(1,0));
                            break;
                        default:
                            #ifdef GRAPHCUT_DBG
                                cout<<"path dividing error"<<endl;
                            #endif
                            break;
                        }
                        break;
                        
                    case 2:
                        switch(a){
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
                            left->insert(current(1,1));
                            right->insert(current(0,1));
                            break;
                        default:
                            #ifdef GRAPHCUT_DBG
                                cout<<"path dividing error"<<endl;
                            #endif
                            break;
                        }
                        break;
                        
                    case 3:
                        switch(a){
                        case 2:
                            left->insert(current(1,0));
                            left->insert(current(1,1));
                            break;
                        case 0:
                            right->insert(current(1,0));
                            right->insert(current(1,1));
                            break;
                        case 1:
                            right->insert(current(1,1));
                            left->insert(current(1,0));
                            break;
                        case 3:
                            right->insert(current(0,0));
                            left->insert(current(0,1));
                            break;
                        default:
                            #ifdef GRAPHCUT_DBG
                                cout<<"path dividing error"<<endl;
                            #endif
                            break;
                        }
                        break;
                }
                
                switch(closest){
                    case 0:
                        tmp = current;
                        while(tmp.y + iBB.upperLeft().y > -2){
                            tmp.y--;
                            left->insert(tmp(1,0));
                            right->insert(tmp);
                            }
                        break;
                    case 1:
                        tmp = current;
                        while(tmp.x < iBB.lowerRight().x + 2){
                            tmp.x++;
                            left->insert(tmp(0,1));
                            right->insert(tmp);
                            }
                        break;
                    case 2:
                        tmp = current;
                        while(tmp.y < iBB.lowerRight().y + 2){
                            tmp.y++;
                            left->insert(tmp);
                            right->insert(tmp(1,0));
                            }
                        break;
                    case 3:
                        tmp = current;
                        while(tmp.x > -2){
                            tmp.x--;
                            left->insert(tmp);
                            right->insert(tmp(0,1));
                            }
                        break;
                    default:
                        #ifdef GRAPHCUT_DBG
                            cout<<"path dividing error"<<endl;
                        #endif
                        break;

                }
            }
            else{
                
                int a, b;
                
                if(previous.x == current.x){
                    if(previous.y > current.y)
                        a = 0;
                    else a = 2;
                } else{
                    if(previous.x < current.x)
                        a = 1;
                    else a = 3;
                }
                
                a = a << 2;
                
                if(nextDual == cut->end()){
                    int dist;
                    Diff2D toLowerRight = dim - current;
                    closest = 0;
                    dist = current.y;
                    if(toLowerRight.x < dist){
                        closest = 1;
                        dist = toLowerRight.x;
                    }
                    if(toLowerRight.y < dist){
                        closest = 2;
                        dist = toLowerRight.y;
                    }
                    if(current.x < dist){
                        closest = 3;
                    }
                    
                    b = closest;
                } else if(next.x == current.x){
                    if(next.y < current.y)
                        b = 0;
                    else b = 2;
                } else{
                    if(next.x > current.x)
                        b = 1;
                    else b = 3;
                }
                
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
                        cout<<"path dividing error: two-point loop"<<endl;
                    #endif
                    break;
                }
                
                if(*currentDual == cut->back()){
                  
                    switch(closest /*^ BIT_MASK_OPDIR*/){
                        case 2:
                            tmp = current;
                            switch(b){
                                case 1:
                                    left->insert(current);
                                    left->insert(current(1,0));
                                    right->insert(current(0,1));
                                    left->insert(current(1,1));
                                    break;
                                case 3:
                                    right->insert(current);
                                    right->insert(current(1,0));
                                    right->insert(current(0,1));
                                    left->insert(current(1,1));
                                    break;
                            }
                            while(tmp.y < iBB.lowerRight().y + 2){
                                tmp.y++;
                                left->insert(tmp(1,0));
                                right->insert(tmp);
                                }
                            break;
                        case 3:
                            tmp = current;
                            switch(b){
                                case 0:
                                    right->insert(current);
                                    right->insert(current(1,0));
                                    left->insert(current(0,1));
                                    right->insert(current(1,1));
                                    break;
                                case 2:
                                    right->insert(current);
                                    left->insert(current(1,0));
                                    left->insert(current(0,1));
                                    left->insert(current(1,1));
                                    break;
                            }
                            
                            while(tmp.x > -2){
                                tmp.x--;
                                left->insert(tmp(0,1));
                                right->insert(tmp);
                                }
                            break;
                        case 0:
                            tmp = current;
                            switch(b){
                                case 1:
                                    left->insert(current);
                                    right->insert(current(1,0));
                                    right->insert(current(0,1));
                                    right->insert(current(1,1));
                                    break;
                                case 3:
                                    left->insert(current);
                                    right->insert(current(1,0));
                                    left->insert(current(0,1));
                                    left->insert(current(1,1));
                                    break;
                            }
                            while(tmp.y + iBB.upperLeft().y  > -2){
                                tmp.y--;
                                left->insert(tmp);
                                right->insert(tmp(1,0));
                                }
                            break;
                        case 1:
                            tmp = current;
                            switch(b){
                                case 0:
                                    left->insert(current);
                                    left->insert(current(1,0));
                                    left->insert(current(0,1));
                                    right->insert(current(1,1));
                                    break;
                                case 2:
                                    right->insert(current);
                                    left->insert(current(1,0));
                                    right->insert(current(0,1));
                                    right->insert(current(1,1));
                                    break;
                            }
                            while(tmp.x < iBB.lowerRight().x + 2){
                                tmp.x++;
                                left->insert(tmp);
                                right->insert(tmp(0,1));
                                }
                            break;
                        default:
                            #ifdef GRAPHCUT_DBG
                                cout<<"path dividing error"<<endl;
                            #endif
                            break;
                        }
                }
            }
            
            if(nextDual != cut->end()) ++nextDual;
            
            previous = current;
            previousDual = currentDual;
        }
    }
    /* Graph-cut
     * implementation of an algorithm based on an article by:
     * V.Kwatra, A.Schoedl, I.Essa, G.Turk, A.Bobick
     * "Graphcut Textures: Image and Video Synthesis Using Graph Cuts"
     * 
     * The algorithm finds the best seam between two images using a graph-based approach.
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
        typedef float GradientPixelType;
        typedef vigra::NumericTraits<BasePixelType> BasePixelTraits;
        typedef vigra::NumericTraits<GradientPixelType> GradientPixelTraits;
        typedef typename BasePixelTraits::Promote BasePromotePixelType;
        typedef vigra::NumericTraits<BasePromotePixelType> BasePromotePixelTraits;
        typedef typename BasePromotePixelTraits::Promote GraphPixelType;

        const Diff2D size(src1_lowerright.x - src1_upperleft.x,
                          src1_lowerright.y - src1_upperleft.y),
                          masksize(mask1_lowerright.x - mask1_upperleft.x,
                          mask1_lowerright.y - mask1_upperleft.y);

        IMAGETYPE<BasePixelType> tmp(size);
        IMAGETYPE<BasePromotePixelType> graphtmp(size + size + Diff2D(1,1)), tmp2(size);
        IMAGETYPE<GradientPixelType> gradientX(size), gradientY(size);
        IMAGETYPE<GraphPixelType> graph(size + size + Diff2D(1,1));
        IMAGETYPE<MaskPixelType> finalmask(size+Diff2D(2,2));
	
        vector<Point2D> *dualPath = NULL;
        vector<Point2D> totalDualPath;
        vector<Point2D> *list;
        Point2D tmpPoint;
        boost::unordered_set<Point2D, pointHash> a, b;
        CheckpointPixels* bounds = new CheckpointPixels();
        
        
        const Diff2D graphsize(graph.lowerRight().x - graph.upperLeft().x,
                          graph.lowerRight().y - graph.upperLeft().y);
        const Rect2D gBB(Point2D(1,1), Size2D(graphsize)),
                     bBB(Point2D(1,1), Size2D(size));
        Kernel2D<BasePromotePixelType> edgeWeightKernel;
        edgeWeightKernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
                         0, 1, 0,
                         1, 1, 1,
                         0, 1, 0;
        edgeWeightKernel.setBorderTreatment(vigra::BORDER_TREATMENT_CLIP);
        Kernel1D<float> gradientKernel;
        gradientKernel.initSymmetricGradient();
        //gradient images calculation
        transformImage(src1_upperleft, src1_lowerright, sa1,
                tmp2.upperLeft(), tmp2.accessor(),
                MapFunctor<SrcPixelType, BasePixelType>()
                ); 
        transformImage(src2_upperleft, src1_lowerright, sa2,
                tmp.upperLeft(), tmp.accessor(),
                MapFunctor<SrcPixelType, BasePixelType>()
                ); 
        combineTwoImagesMP(srcImageRange(tmp2), srcImage(tmp), destImage(tmp2), Arg1()+Arg2());
        #ifdef GRAPHCUT_DBG
                exportImage(srcImageRange(tmp), ImageExportInfo("./debug/sum.tif").setPixelType("UINT8"));
                exportImage(srcImageRange(tmp2), ImageExportInfo("./debug/sum2.tif").setPixelType("UINT8"));
        #endif
        //computing image gradient for a better cost function
        vigra::separableConvolveX(srcImageRange(tmp2), destImage(gradientX), kernel1d(gradientKernel));
        vigra::separableConvolveY(srcImageRange(tmp2), destImage(gradientY), kernel1d(gradientKernel));

        //difference image calculation
        combineTwoImagesMP(src1_upperleft, src1_lowerright, sa1, 
                src2_upperleft, sa2,
                tmp.upperLeft(), tmp.accessor(),
                PixelDifferenceFunctor<SrcPixelType, BasePixelType>()
                );     
        
        //masking overlap region borders
        combineThreeImagesMP(iBB.apply(srcIterRange(mask1_upperleft, mask1_lowerright, ma1)), 
                iBB.apply(srcIter(mask2_upperleft, ma2)),
                srcIter(tmp.upperLeft(), tmp.accessor()),
                destIter(tmp.upperLeft(), tmp.accessor()),
                ifThenElse(!(Arg1() & Arg2()), Param(BasePixelTraits::max()), Arg3())
                );
        
        
        //look for possible start and end points     
        list = findExtremePoints<MaskImageIterator, MaskAccessor, BasePixelType>
                (mask1_upperleft, mask1_lowerright, ma1,
                            mask2_upperleft, ma2, dest_upperleft, da,
                            norm, boundary, iBB);
				
		//in case something goes wrong with point finding, returns regular nft image
		if(list->empty())
			return;
        
        //copying to a grid
        copyImage(srcImageRange(tmp), stride(2,2,gBB.apply(destImage(graphtmp))));
                     
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(tmp), ImageExportInfo("./debug/diff.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(graphtmp), ImageExportInfo("./debug/diff2.tif").setPixelType("UINT8"));
        exportImage(mask1_upperleft, mask1_lowerright, ma1, ImageExportInfo("./debug/mask1.tif").setPixelType("UINT8"));
        exportImage(mask2_upperleft, mask1_lowerright, ma2, ImageExportInfo("./debug/mask2.tif").setPixelType("UINT8"));
        exportImage(src1_upperleft, src1_lowerright, sa1, ImageExportInfo("./debug/src1.tif").setPixelType("UINT8"));
        exportImage(src2_upperleft, src1_lowerright, sa2, ImageExportInfo("./debug/src2.tif").setPixelType("UINT8"));
   
        exportImage(srcImageRange(gradientX), ImageExportInfo("./debug/gradx.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(gradientY), ImageExportInfo("./debug/grady.tif").setPixelType("UINT8"));
#endif
        //calculating differences between pixels that are adjacent in the original image
        convolveImage(srcImageRange(graphtmp), destImage(graph), kernel2d(edgeWeightKernel));
        copyImage(srcImageRange(graph), destImage(graphtmp));
        
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(graph), ImageExportInfo("./debug/graph.tif").setPixelType("UINT8"));
#endif
        //find optimal cut in dual graph
        for(vector<Point2D>::iterator i = list->begin(); i != list->end(); ++i){
            
            if(i == list->begin()){
                tmpPoint = *i;
                continue;
            }
            
            bounds->clear();
            bounds->top.insert(tmpPoint);
            bounds->bottom.insert(*i);
            
             #ifdef GRAPHCUT_DBG
                    cout<<"Running graph-cut: "<<tmpPoint<<":"<<*i<<endl;
             #endif
        
            dualPath = A_star<IMAGETYPE<GraphPixelType>, IMAGETYPE<GradientPixelType>, BasePixelType>(
                Point2D(-10,-10), 
                Point2D(-20,-20),
                &graphtmp, &gradientX, &gradientY, graphsize - Diff2D(1,1),
                bounds);
            
            for(vector<Point2D>::reverse_iterator j = dualPath->rbegin(); j < dualPath->rend(); j++){
                if((j == dualPath->rbegin() && totalDualPath.empty()) || j != dualPath->rbegin())
                    totalDualPath.push_back(*j);
            }
            copyImage(srcImageRange(graph), destImage(graphtmp));
            tmpPoint = *i;
        }
        
        dividePath<IMAGETYPE<GraphPixelType> >(&graph, &totalDualPath, &a, &b, graphsize - Diff2D(1,1), iBB);
        
        //labels areas that belong to left/right images
        //adds a 1-pixel border to catch any area that is cut off by the seam
        labelImage(srcIterRange(Diff2D(), Diff2D() + size + Diff2D(2,2)), 
                        destImage(finalmask), 
                        false, PathEqualityFunctor(&a, &b, iBB.upperLeft()));
        
        #ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(finalmask), ImageExportInfo("./debug/labels.tif").setPixelType("UINT8"));
#endif
        
        int colorSum=0, count=0;
        CountFunctor<MaskPixelType> c(&colorSum, &count); 
        
        transformImage(srcIterRange(finalmask.upperLeft()+Diff2D(1,1), 
                        finalmask.lowerRight()-Diff2D(1,1)), 
                        destIter(finalmask.upperLeft()+Diff2D(1,1)),
                        ifThenElse(Arg1() == Param(1), Param(BasePixelTraits::max()), 
                        Param(BasePixelTraits::zero())));
        inspectTwoImages(srcIterRange(finalmask.upperLeft()+Diff2D(1,1), 
                        finalmask.lowerRight()-Diff2D(1,1)),
                        srcIter(dest_upperleft + iBB.upperLeft(), da),
                        c);
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(finalmask), ImageExportInfo("./debug/labels1.tif").setPixelType("UINT8"));
#endif
        
        //if colorSum < half the pixels labeled as "1" in finalmask should be black
        if(colorSum < count / 2){
            transformImage(srcIterRange(finalmask.upperLeft()+Diff2D(1,1), 
                        finalmask.lowerRight()-Diff2D(1,1)), 
                        destIter(finalmask.upperLeft()+Diff2D(1,1)),
                        ifThenElse(Arg1() == Param(0), Param(BasePixelTraits::max()), 
                        Param(BasePixelTraits::zero())));
        } else{
            transformImage(srcIterRange(finalmask.upperLeft()+Diff2D(1,1), 
                        finalmask.lowerRight()-Diff2D(1,1)), 
                        destIter(finalmask.upperLeft()+Diff2D(1,1)),
                        ifThenElse(Arg1() == Param(0), Param(BasePixelTraits::zero()), 
                        Param(BasePixelTraits::max())));
        }
        
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(finalmask), ImageExportInfo("./debug/labels2.tif").setPixelType("UINT8"));
#endif
	combineTwoImagesMP(iBB.apply(srcIterRange(mask1_upperleft, mask1_lowerright, ma1)),
                            iBB.apply(srcIter(mask2_upperleft, ma2)),
                            destImage(gradientX),
                            ifThenElse(Arg1()&&Arg2(), Param(MaskPixelTraits::max()), Param(MaskPixelTraits::zero())));
	
	copyImageIf(srcIterRange(finalmask.upperLeft()+Diff2D(1,1), 
                        finalmask.lowerRight()-Diff2D(1,1)),
			srcImage(gradientX),
                        destIter(dest_upperleft + iBB.upperLeft(), da));
        
        delete list;
        delete dualPath;
    }
    

}/* namespace enblend */
#endif	/* GRAPHCUT_H */

