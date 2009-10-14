/*
 * Copyright (C) 2004-2009 Andrew Mihal
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
#ifndef __MASK_H__
#define __MASK_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <functional>
#include <numeric>
#ifdef HAVE_EXT_SLIST
#include <ext/slist>
#else
#include <slist>
#endif

#include "common.h"
#include "anneal.h"
#include "nearest.h"
#include "path.h"

#include "vigra/contourcirculator.hxx"
#include "vigra/error.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/impex.hxx"
#include "vigra/initimage.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/stdcachedfileimage.hxx"
#include "vigra_ext/impexalpha.hxx"
#include "vigra_ext/XMIWrapper.h"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/algorithm.hpp>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/if.hpp>

using std::make_pair;
using std::vector;
#ifdef HAVE_EXT_SLIST
using __gnu_cxx::slist;
#else
using std::slist;
#endif

using vigra::combineThreeImages;
using vigra::combineTwoImages;
using vigra::CrackContourCirculator;
using vigra::exportImage;
using vigra::ImageExportInfo;
using vigra::ImageImportInfo;
using vigra::importImageAlpha;
using vigra::initImageIf;
using vigra::LinearIntensityTransform;
using vigra::linearRangeMapping;
using vigra::NumericTraits;
using vigra::Point2D;
using vigra::Rect2D;
using vigra::RGBToGrayAccessor;
using vigra::Size2D;
using vigra::transformImage;
using vigra::transformImageIf;

using vigra::functor::Arg1;
using vigra::functor::Arg2;
using vigra::functor::Arg3;
using vigra::functor::ifThenElse;
using vigra::functor::Param;
using vigra::functor::UnaryFunctor;

using vigra_ext::copyPaintedSetToImage;

using boost::lambda::_1;
using boost::lambda::_2;
using boost::lambda::if_then_else_return;
using boost::lambda::call_begin;
using boost::lambda::call_end;
using boost::lambda::constant;
using boost::lambda::protect;
using boost::lambda::ret;

namespace enblend {

/** Data structures for vector masks */
typedef pair<bool, Point2D> SegmentPoint;
typedef slist<SegmentPoint> Segment;
typedef vector<Segment*> Contour;
typedef vector<Contour*> ContourVector;


void
dump_segment(const Segment& segment,
             const std::string& prefix = "", std::ostream& out = std::cout)
{
    const unsigned points_per_line = 5U;
    unsigned n = 1U;
    const size_t size = segment.size();

    out << prefix << "{segment with " << size << " point(s):\n" << prefix;
    for (Segment::const_iterator i = segment.begin(); i != segment.end(); ++i, ++n)
    {
        out << ' ' << i->second << (i->first ? 'f' : ' ');
        if (n % points_per_line == 0U && n != size)
        {
            out << '\n' << prefix;
        }
    }
    out << "\n" << prefix << "}\n";
}


void
dump_contour(const Contour& contour,
             const std::string& prefix = "", std::ostream& out = std::cout)
{
    out << prefix << "{contour with " << contour.size() << " segment(s):\n";
    for (Contour::const_iterator i = contour.begin(); i != contour.end(); ++i)
    {
        dump_segment(**i, prefix + "    ", out);
    }
    out << prefix << "}\n";
}


void
dump_contourvector(const ContourVector& contourvector,
                   const std::string& prefix = "", std::ostream& out = std::cout)
{
    out << prefix << "{contourvector with " << contourvector.size() << " contour(s):\n";
    for (ContourVector::const_iterator i = contourvector.begin(); i != contourvector.end(); ++i)
    {
        dump_contour(**i, prefix + "    ", out);
    }
    out << prefix << "}\n";
}


// Draw a line from begin to end in image.  Use value as pixel color.
//
// The algorithm is not implemented efficiently.  However, the
// function is only called to visualize the initial seam line, which
// we do not consider a performance-critical job.
template <typename ImageType>
void
visualizeLine(ImageType& image,
              const Point2D& begin, const Point2D& end,
              typename ImageType::PixelType value)
{
    typedef typename ImageType::Iterator Iterator;

    const vigra::Diff2D difference(end - begin);
    const int stepX = sign(difference.x);
    const int stepY = sign(difference.y);
    Iterator point(image.upperLeft() + begin);
    const Iterator stop(image.upperLeft() + end);
    double error = 0.0;

    //std::cout << "+ [" << begin << " .. " << end << "]\n";

    const vigra::Size2D size(image.size());
    // Exit immediately if both start point and end point are outside of the image.
    if (!(begin.x >= 0 && begin.x < size.x && begin.y >= 0 && begin.y < size.y &&
          end.x >= 0 && end.x < size.x && end.y >= 0 && end.y < size.y)) {
        return;
    }

    if (abs(difference.x) >= abs(difference.y)) {
        const double delta_error =
            difference.x == 0 ?
            0.0 :
            fabs(static_cast<double>(difference.y) / static_cast<double>(difference.x));
        while (true) {
            if (point.x >= image.upperLeft().x && point.x < image.lowerRight().x &&
                point.y >= image.upperLeft().y && point.y < image.lowerRight().y) {
                image.accessor().set(value, point);
            }
            if (point.x == stop.x) {
                break;
            }
            error += delta_error;
            if (fabs(error) >= 0.5) {
                point.y += stepY;
                error -= 1.0;
            }
            point.x += stepX;
        }
    } else {
        const double delta_error =
            difference.y == 0 ?
            0.0 :
            fabs(static_cast<double>(difference.x) / static_cast<double>(difference.y));
        while (true) {
            if (point.x >= image.upperLeft().x && point.x < image.lowerRight().x &&
                point.y >= image.upperLeft().y && point.y < image.lowerRight().y) {
                image.accessor().set(value, point);
            }
            if (point.y == stop.y) {
                break;
            }
            error += delta_error;
            if (fabs(error) >= 0.5) {
                point.x += stepX;
                error -= 1.0;
            }
            point.y += stepY;
        }
    }
}


template <typename PixelType, typename ResultType>
class PixelDifferenceFunctor
{
    typedef typename EnblendNumericTraits<PixelType>::ImagePixelComponentType PixelComponentType;
    typedef typename EnblendNumericTraits<ResultType>::ImagePixelComponentType ResultPixelComponentType;
    typedef LinearIntensityTransform<ResultType> RangeMapper;

public:
    PixelDifferenceFunctor() :
        rm(linearRangeMapping(NumericTraits<PixelComponentType>::min(),
                              NumericTraits<PixelComponentType>::max(),
                              ResultType(NumericTraits<ResultPixelComponentType>::min()),
                              ResultType(NumericTraits<ResultPixelComponentType>::max()))) {}

    ResultType operator()(const PixelType& a, const PixelType& b) const {
        typedef typename NumericTraits<PixelType>::isScalar src_is_scalar;
        return diff(a, b, src_is_scalar());
    }

protected:
    ResultType diff(const PixelType& a, const PixelType& b, VigraFalseType) const {
        PixelComponentType aLum = a.luminance();
        PixelComponentType bLum = b.luminance();
        PixelComponentType aHue = a.hue();
        PixelComponentType bHue = b.hue();
        PixelComponentType lumDiff = (aLum > bLum) ? (aLum - bLum) : (bLum - aLum);
        PixelComponentType hueDiff = (aHue > bHue) ? (aHue - bHue) : (bHue - aHue);
        if (hueDiff > (NumericTraits<PixelComponentType>::max() / 2)) {
            hueDiff = NumericTraits<PixelComponentType>::max() - hueDiff;
        }
        return rm(std::max(hueDiff, lumDiff));
    }

    ResultType diff(const PixelType& a, const PixelType& b, VigraTrueType) const {
        typedef typename NumericTraits<PixelType>::isSigned src_is_signed;
        return scalar_diff(a, b, src_is_signed());
    }

    ResultType scalar_diff(const PixelType& a, const PixelType& b, VigraTrueType) const {
        return rm(std::abs(a - b));
    }

    // This appears necessary because NumericTraits<unsigned int>::Promote
    // is an unsigned int instead of an int.
    ResultType scalar_diff(const PixelType& a, const PixelType& b, VigraFalseType) const {
        return rm(std::abs(static_cast<int>(a) - static_cast<int>(b)));
    }

    RangeMapper rm;
};


template <typename ValueType, typename AccessorType>
class XorAccessor
{
public:
    typedef ValueType value_type;
    typedef AccessorType accessor_type;

    XorAccessor(AccessorType a) : acc(a) {}

    template <class Iterator>
    ValueType operator()(const Iterator& i) const {return acc(i);}

    ValueType operator()(const ValueType* i) const {return acc(i);}

    template <class Iterator, class Difference>
    ValueType operator()(const Iterator& i, Difference d) const {return acc(i, d);}

    template <class Value, class Iterator>
    void set(const Value& v, const Iterator& i) const {acc.set(v ^ acc(i), i);}

    template <class Value, class Iterator>
    void set(const Value& v, Iterator& i) const {acc.set(v ^ acc(i), i);}

    template <class Value, class Iterator, class Difference>
    void set(const Value& v, const Iterator& i, const Difference& d) const {acc.set(v ^ acc(i, d), i, d);}

private:
    AccessorType acc;
};


template <typename MaskType>
void fillContour(MaskType* mask, const Contour& contour, Diff2D offset)
{
    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::Accessor MaskAccessor;
    typedef NumericTraits<MaskPixelType> MaskPixelTraits;

    const size_t totalPoints =
        std::accumulate(contour.begin(), contour.end(),
                        0U, ret<size_t>(_1 + bind(&Segment::size, _2)));

    if (totalPoints == 0U) {
        return;
    }

    miPixel pixels[2] = {MaskPixelTraits::max(), MaskPixelTraits::max()};
    miGC* pGC = miNewGC(2, pixels);
    miPoint* const points = new miPoint[totalPoints];
    miPoint* p = points;

    for (Contour::const_iterator currentSegment = contour.begin();
         currentSegment != contour.end();
         ++currentSegment) {
        for (Segment::iterator vertexIterator = (*currentSegment)->begin();
             vertexIterator != (*currentSegment)->end();
             ++vertexIterator, ++p) {
            p->x = vertexIterator->second.x;
            p->y = vertexIterator->second.y;
        }
    }

    miPaintedSet* paintedSet = miNewPaintedSet();
    miFillPolygon(paintedSet, pGC,
                  MI_SHAPE_GENERAL, MI_COORD_MODE_ORIGIN,
                  totalPoints, points);

    delete [] points;

    copyPaintedSetToImage(mask->upperLeft(), mask->lowerRight(),
                          XorAccessor<MaskPixelType, MaskAccessor>(mask->accessor()),
                          paintedSet, offset);

    miDeletePaintedSet(paintedSet);
    miDeleteGC(pGC);
}


template <typename MaskType>
void maskBounds(MaskType* mask, const Rect2D& uBB, Rect2D& mBB)
{
    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;
    typedef typename MaskType::Accessor MaskAccessor;

    // Find the bounding box of the mask transition line and put it in mBB.
    // mBB starts out as empty rectangle.
    mBB = Rect2D(Point2D(mask->size()), Point2D(0, 0));

    MaskIteratorType myPrev = mask->upperLeft();
    MaskIteratorType my = mask->upperLeft() + Diff2D(0, 1);
    MaskIteratorType mend = mask->lowerRight();
    MaskIteratorType mxLeft = myPrev;
    MaskIteratorType mx = myPrev + Diff2D(1, 0);

    for (int x = 1; mx.x < mend.x; ++x, ++mx.x, ++mxLeft.x) {
        if (*mxLeft != *mx) {
            mBB |= Rect2D(x - 1, 0, x + 1, 1);
        }
    }
    for (int y = 1; my.y < mend.y; ++y, ++my.y, ++myPrev.y) {
        mxLeft = my;
        mx = my + Diff2D(1, 0);
        MaskIteratorType mxUpLeft = myPrev;
        MaskIteratorType mxUp = myPrev + Diff2D(1, 0);

        if (*mxUpLeft != *mxLeft) {
            // Transition line is between mxUpLeft and mxLeft.
            mBB |= Rect2D(0, y - 1, 1, y + 1);
        }

        for (int x = 1; mx.x < mend.x; ++x, ++mx.x, ++mxLeft.x, ++mxUp.x) {
            if (*mxLeft != *mx) {
                mBB |= Rect2D(x - 1, y, x + 1, y + 1);
            }
            if (*mxUp != *mx) {
                mBB |= Rect2D(x, y - 1, x + 1, y + 1);
            }
        }
    }

    // Check that mBB is well-defined.
    if (mBB.isEmpty()) {
        // No transition pixels were found in the mask at all.  This
        // means that one image has no contribution.
        if (*(mask->upperLeft()) == NumericTraits<MaskPixelType>::zero()) {
            // If the mask is entirely black, then inspectOverlap
            // should have caught this.  It should have said that the
            // white image is redundant.
            cerr << command
                 << ": mask is entirely black, but white image was not identified as redundant"
                 << endl;
            exit(1);
        } else {
            // If the mask is entirely white, then the black image
            // would have been identified as redundant if black and
            // white were swapped.  Set mBB to the full size of the
            // mask.
            mBB = uBB;
            // Explain why the black image disappears completely.
            cerr << command
                 << ": warning: previous images are completely overlapped by the current images"
                 << endl;
        }
    } else {
        // mBB is defined relative to inputUnion origin
        //cerr << command << ": info: mBB relative to mask: " << mBB << endl;
        mBB.moveBy(uBB.upperLeft());
    }

    if (Verbose >= VERBOSE_ROIBB_SIZE_MESSAGES) {
        cerr << command
             << ": info: mask transition line bounding box: "
             << mBB
             << endl;
    }
}


/** Vectorize the seam line defined in nftOutputImage into the contour
 *  rawSegments. */
template <typename MaskType, typename AlphaType>
void vectorizeSeamLine(Contour& rawSegments,
                       const AlphaType* const whiteAlpha, const AlphaType* const blackAlpha,
                       const Rect2D& uBB,
                       int nftStride, MaskType* nftOutputImage)
{
    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;

    const double diagonalLength =
        hypot(static_cast<double>(nftOutputImage->width()),
              static_cast<double>(nftOutputImage->height()));
    int vectorizeDistance =
        MaskVectorizeDistance.isPercentage ?
        static_cast<int>(ceil(MaskVectorizeDistance.value / 100.0 * diagonalLength)) :
        MaskVectorizeDistance.value;
    if (vectorizeDistance < minimumVectorizeDistance) {
        cerr << command
             << ": warning: mask vectorization distance "
             << vectorizeDistance
             << " ("
             << 100.0 * vectorizeDistance / diagonalLength
             << "% of diagonal) is smaller\n"
             << command
             << ": warning:   than minimum of " << minimumVectorizeDistance
             << "; will use " << minimumVectorizeDistance << " ("
             << 100.0 * minimumVectorizeDistance / diagonalLength
             << "% of diagonal)"
             << endl;
        vectorizeDistance = minimumVectorizeDistance;
    }

    const Rect2D border(1, 1, nftOutputImage->width() - 1, nftOutputImage->height() - 1);

    MaskIteratorType my = nftOutputImage->upperLeft() + Diff2D(1, 1);
    MaskIteratorType mend = nftOutputImage->lowerRight() + Diff2D(-1, -1);
    for (int y = 1; my.y < mend.y; ++y, ++my.y) {
        MaskIteratorType mx = my;
        MaskPixelType lastColor = NumericTraits<MaskPixelType>::zero();

        for (int x = 1; mx.x < mend.x; ++x, ++mx.x) {
            if (*mx == NumericTraits<MaskPixelType>::max()
                && lastColor == NumericTraits<MaskPixelType>::zero()) {
                // Found the corner of a previously unvisited white region.
                // Create a Segment to hold the border of this region.
                vector<Point2D> excessPoints;
                Segment* snake = new Segment();
                rawSegments.push_back(snake);

                // Walk around border of white region.
                CrackContourCirculator<MaskIteratorType> crack(mx);
                CrackContourCirculator<MaskIteratorType> crackEnd(crack);
                bool lastPointFrozen = false;
                int distanceLastPoint = 0;
                do {
                    Point2D currentPoint = *crack + Diff2D(x, y);
                    crack++;
                    Point2D nextPoint = *crack + Diff2D(x, y);

                    // See if currentPoint lies on border.
                    if (currentPoint.x == border.left()
                        || currentPoint.x == border.right()
                        || currentPoint.y == border.top()
                        || currentPoint.y == border.bottom()) {
                        // See if currentPoint is in a corner.
                        if ((currentPoint.x == border.left() && currentPoint.y == border.top())
                            || (currentPoint.x == border.left() && currentPoint.y == border.bottom())
                            || (currentPoint.x == border.right() && currentPoint.y == border.top())
                            || (currentPoint.x == border.right() && currentPoint.y == border.bottom())) {
                            snake->push_front(make_pair(false, currentPoint));
                            distanceLastPoint = 0;
                        } else if (!lastPointFrozen
                                   || (nextPoint.x != border.left()
                                       && nextPoint.x != border.right()
                                       && nextPoint.y != border.top()
                                       && nextPoint.y != border.bottom())) {
                            snake->push_front(make_pair(false, currentPoint));
                            distanceLastPoint = 0;
                        } else {
                            excessPoints.push_back(currentPoint);
                        }
                        lastPointFrozen = true;
                    } else {
                        // Current point is not frozen.
                        if (distanceLastPoint % vectorizeDistance == 0) {
                            snake->push_front(make_pair(true, currentPoint));
                            distanceLastPoint = 0;
                        } else {
                            excessPoints.push_back(currentPoint);
                        }
                        lastPointFrozen = false;
                    }
                    distanceLastPoint++;
                } while (crack != crackEnd);

                // Paint the border so this region will not be found again
                for (Segment::iterator vertexIterator = snake->begin();
                     vertexIterator != snake->end(); ++vertexIterator) {
                    (*nftOutputImage)[vertexIterator->second] = NumericTraits<MaskPixelType>::one();

                    // While we're at it, convert vertices to uBB-relative coordinates.
                    vertexIterator->second =
                        nftStride * (vertexIterator->second + Diff2D(-1, -1));

                    // While we're at it, mark vertices outside the union region as not moveable.
                    if (vertexIterator->first
                        && (*whiteAlpha)[vertexIterator->second + uBB.upperLeft()] == NumericTraits<MaskPixelType>::zero()
                        && (*blackAlpha)[vertexIterator->second + uBB.upperLeft()] == NumericTraits<MaskPixelType>::zero()) {
                        vertexIterator->first = false;
                    }
                }
                for (vector<Point2D>::iterator vertexIterator = excessPoints.begin();
                     vertexIterator != excessPoints.end(); ++vertexIterator) {
                    // These are points on the border of the white
                    // region that are not in the snake.  Recolor them
                    // so that this white region will not be found
                    // again.
                    (*nftOutputImage)[*vertexIterator] = NumericTraits<MaskPixelType>::one();
                }
            }

            lastColor = *mx;
        }
    }
}


/** Convert rawContours snakes into segments with unbroken runs of
 *  moveable vertices. */
void reorderSnakesToMovableRuns(ContourVector& contours, const Contour& rawSegments)
{
    for (Contour::const_iterator segments = rawSegments.begin();
         segments != rawSegments.end();
         ++segments) {
        Segment* snake = *segments;

        // Snake becomes multiple separate segments in one contour
        Contour* currentContour = new Contour();
        contours.push_back(currentContour);

        // Check if snake is a closed contour
        bool closedContour = true;
        Segment::iterator vertexIterator = snake->begin();
        for (Segment::iterator vertexIterator = snake->begin();
             vertexIterator != snake->end();
             ++vertexIterator) {
            if (!vertexIterator->first) {
                closedContour = false;
                break;
            }
        }

        // Closed contours consist of only moveable vertices.
        if (closedContour) {
            currentContour->push_back(snake);
            continue;
        }

        if (snake->front().first) {
            // First vertex is moveable.  Rotate list so that first
            // vertex is nonmoveable.
            Segment::iterator firstNonmoveableVertex = snake->begin();
            while (firstNonmoveableVertex->first) {
                ++firstNonmoveableVertex;
            }

            // Copy initial run on moveable vertices and first
            // non-moveable vertex to end of list.
            Segment::iterator firstNonmoveableSuccessor = firstNonmoveableVertex;
            ++firstNonmoveableSuccessor;
            if (EXPECT_RESULT(firstNonmoveableSuccessor == snake->end(), false)) {
                snake->insert_after(firstNonmoveableVertex,
                                    snake->begin(), firstNonmoveableVertex);
            } else {
                snake->insert(snake->end(), // append at the end
                              snake->begin(), firstNonmoveableSuccessor);
            }

            // Erase initial run of moveable vertices.
            snake->erase(snake->begin(), firstNonmoveableVertex);
        }

        // Find last moveable vertex.
        Segment::iterator lastMoveableVertex = snake->begin();
        for (Segment::iterator vertexIterator = snake->begin();
             vertexIterator != snake->end();
             ++vertexIterator) {
            if (vertexIterator->first) {
                lastMoveableVertex = vertexIterator;
            }
        }

        Segment* currentSegment = NULL;
        bool insideMoveableSegment = false;
        bool passedLastMoveableVertex = false;
        Segment::iterator lastNonmoveableVertex = snake->begin();
        for (Segment::iterator vertexIterator = snake->begin();
             vertexIterator != snake->end();
             ++vertexIterator) {
            // Create a new segment if necessary.
            if (currentSegment == NULL) {
                currentSegment = new Segment();
                currentContour->push_back(currentSegment);
            }

            // Keep track of when we visit the last moveable vertex.
            // Don't create new segments after this point.
            // Add all remaining nonmoveable vertices to current segment.
            if (vertexIterator == lastMoveableVertex) {
                passedLastMoveableVertex = true;
            }

            // Keep track of last nonmoveable vertex.
            if (!vertexIterator->first) {
                lastNonmoveableVertex = vertexIterator;
            }

            // All segments must begin with a nonmoveable vertex.
            // If only one nonmoveable vertex separates two runs of moveable vertices,
            // that vertex is copied into the beginning of the current segment.
            // It was previously added at the end of the last segment.
            if (vertexIterator->first && currentSegment->empty()) {
                currentSegment->push_front(*lastNonmoveableVertex);
            }

            // Add the current vertex to the current segment.
            currentSegment->push_front(*vertexIterator);

            if (!insideMoveableSegment && vertexIterator->first) {
                // Beginning a new moveable segment.
                insideMoveableSegment = true;
            } else if (insideMoveableSegment
                       && !vertexIterator->first
                       && !passedLastMoveableVertex) {
                // End of currentSegment.
                insideMoveableSegment = false;
                // Correct for the push_fronts we've been doing
                currentSegment->reverse();
                // Cause a new segment to be generated on next vertex.
                currentSegment = NULL;
            }
        }

        // Reverse the final segment.
        if (currentSegment != NULL) {
            currentSegment->reverse();
        }

        delete snake;
    }
}


/** Find extent of moveable snake vertices, and vertices bordering
 *  moveable vertices vertex bounding box. */
Rect2D vertexBoundingBox(ContourVector& contours)
{
    Rect2D box;
    bool initializedVBB = false;

    for (ContourVector::iterator currentContour = contours.begin();
         currentContour != contours.end();
         ++currentContour) {
        for (Contour::iterator currentSegment = (*currentContour)->begin();
             currentSegment != (*currentContour)->end();
             ++currentSegment) {
            Segment::iterator lastVertex = (*currentSegment)->begin();
            bool foundFirstMoveableVertex = false;
            for (Segment::iterator vertexIterator = (*currentSegment)->begin();
                 vertexIterator != (*currentSegment)->end();
                 ++vertexIterator) {
                if (vertexIterator->first) {
                    if (!initializedVBB) {
                        box = Rect2D(vertexIterator->second, Size2D(1, 1));
                        initializedVBB = true;
                    } else {
                        box |= vertexIterator->second;
                    }

                    if (!foundFirstMoveableVertex) {
                        box |= lastVertex->second;
                    }

                    foundFirstMoveableVertex = true;
                } else if (foundFirstMoveableVertex) {
                    // First nonmoveable vertex at end of run.
                    box |= vertexIterator->second;
                    break;
                }

                lastVertex = vertexIterator;
            }
        }
    }

    return box;
}


/** Calculate a blending mask between whiteImage and blackImage.
 */
template <typename ImageType, typename AlphaType, typename MaskType>
MaskType* createMask(const ImageType* const white,
                     const ImageType* const black,
                     const AlphaType* const whiteAlpha,
                     const AlphaType* const blackAlpha,
                     const Rect2D& uBB,
                     const Rect2D& iBB,
                     bool wraparound,
                     unsigned numberOfImages,
                     list<char*>::const_iterator inputFileNameIterator,
                     unsigned m)
{
    typedef typename ImageType::PixelType ImagePixelType;
    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;
    typedef typename MaskType::Accessor MaskAccessor;

    if (LoadMasks) {
        // Read mask from a file instead of calculating it.
        MaskType* mask = new MaskType(uBB.size());
        const std::string maskFilename =
            enblend::expandFilenameTemplate(LoadMaskTemplate,
                                            numberOfImages,
                                            *inputFileNameIterator,
                                            OutputFileName,
                                            m);
        ImageImportInfo maskInfo(maskFilename.c_str());
        if (Verbose >= VERBOSE_MASK_MESSAGES) {
            cerr << command
                 << ": info: loading mask \"" << maskFilename << "\"" << endl;
        }
        if (maskInfo.width() != uBB.width() || maskInfo.height() != uBB.height()) {
            cerr << command
                 << ": warning: mask in \"" << maskFilename << "\" has size "
                 << "(" << maskInfo.width() << "x" << maskInfo.height() << "),\n"
                 << command
                 << ": warning:     but image union has size " << uBB.size() << ";\n"
                 << command
                 << ": warning:     make sure this is the right mask for the given images"
                 << endl;
        }
        importImage(maskInfo, destImage(*mask));
        return mask;
    }

    // Start by using the nearest feature transform to generate a mask.
    Size2D nftInputSize;
    Rect2D nftInputBB;
    int nftStride;
    if (CoarseMask) {
        // Do NFT at 1/CoarsenessFactor scale.
        // uBB rounded up to multiple of CoarsenessFactor pixels in each direction
        nftInputSize = Size2D((uBB.width() + CoarsenessFactor - 1) / CoarsenessFactor,
                              (uBB.height() + CoarsenessFactor - 1) / CoarsenessFactor);
        nftInputBB = Rect2D(Size2D(uBB.width() / CoarsenessFactor,
                                   uBB.height() / CoarsenessFactor));
        nftStride = CoarsenessFactor;
    } else {
        // Do NFT at 1/1 scale.
        nftInputSize = uBB.size();
        nftInputBB = Rect2D(nftInputSize);
        nftStride = 1;
    }

    Size2D nftOutputSize;
    Diff2D nftOutputOffset;
    if (!CoarseMask && !OptimizeMask) {
        // We are not going to vectorize the mask.
        nftOutputSize = nftInputSize;
        nftOutputOffset = Diff2D(0, 0);
    } else {
        // Add 1-pixel border all around the image for the vectorization algorithm.
        nftOutputSize = nftInputSize + Diff2D(2, 2);
        nftOutputOffset = Diff2D(1, 1);
    }

    // mem usage before: 0
    // mem usage after: CoarseMask: 1/8 * uBB * MaskType
    //                  !CoarseMask: uBB * MaskType
    MaskType* nftOutputImage = new MaskType(nftOutputSize);

    nearestFeatureTransform(stride(nftStride, nftStride, uBB.apply(srcImageRange(*whiteAlpha))),
                            stride(nftStride, nftStride, uBB.apply(srcImage(*blackAlpha))),
                            destIter(nftOutputImage->upperLeft() + nftOutputOffset),
                            ManhattanDistance,
                            wraparound ? HorizontalStrip : OpenBoundaries);

#ifdef DEBUG_NEAREST_FEATURE_TRANSFORM
    {
        typedef pair<const char*, const MaskType*> ImagePair;

        const ImagePair nft[] = {
            std::make_pair("blackmask", blackAlpha),
            std::make_pair("whitemask", whiteAlpha),
            //std::make_pair("nft-input", nftInputImage),
            std::make_pair("nft-output", nftOutputImage)
        };

        for (size_t i = 0; i < sizeof(nft) / sizeof(ImagePair); ++i) {
            const std::string nftMaskTemplate(command + "-" + nft[i].first + "-%n.tif");
            const std::string nftMaskFilename =
                enblend::expandFilenameTemplate(nftMaskTemplate,
                                                numberOfImages,
                                                *inputFileNameIterator,
                                                OutputFileName,
                                                m);
            if (Verbose >= VERBOSE_NFT_MESSAGES) {
                cerr << command
                     << ": info: saving nearest-feature-transform image \""
                     << nftMaskFilename << "\"" << endl;
            }
            ImageExportInfo nftMaskInfo(nftMaskFilename.c_str());
            nftMaskInfo.setCompression(MASK_COMPRESSION);
            exportImage(srcImageRange(*nft[i].second), nftMaskInfo);
        }
    }
#endif

    // mem usage before: CoarseMask: 2/8 * uBB * MaskType
    //                   !CoarseMask: 2 * uBB * MaskType
    // mem usage after: CoarseMask: 1/8 * uBB * MaskType
    //                  !CoarseMask: uBB * MaskType

    if (!CoarseMask && !OptimizeMask) {
        // nftOutputImage is the final mask in this case.
        return nftOutputImage;
    }

    // Vectorize the seam lines found in nftOutputImage.
    Contour rawSegments;
    vectorizeSeamLine(rawSegments,
                      whiteAlpha, blackAlpha,
                      uBB,
                      nftStride, nftOutputImage);
    delete nftOutputImage;

#ifdef DEBUG_SEAM_LINE
    std::cout << "+ createMask: rawSegments\n";
    dump_contour(rawSegments, "+ createMask: ");
#endif

    // mem usage after: 0

    if (!OptimizeMask) {
        // Simply fill contours to get final unoptimized mask.
        MaskType* mask = new MaskType(uBB.size());
        fillContour(mask, rawSegments, Diff2D(0, 0));
        // delete all segments in rawSegments
        std::for_each(rawSegments.begin(), rawSegments.end(), bind(delete_ptr(), _1));
        return mask;
    }

    ContourVector contours;
    reorderSnakesToMovableRuns(contours, rawSegments);
    rawSegments.clear();

#ifdef DEBUG_SEAM_LINE
    std::cout << "+ createMask: contours\n";
    dump_contourvector(contours, "+ createMask: ");
#endif

    {
        const size_t totalSegments =
            std::accumulate(contours.begin(), contours.end(),
                            0U, ret<size_t>(_1 + bind(&Contour::size, _2)));

        if (Verbose >= VERBOSE_MASK_MESSAGES) {
            cerr << command << ": info: optimizing ";
            if (totalSegments == 1U) {
                cerr << "1 distinct seam";
            } else {
                cerr << totalSegments << " distinct seams";
            }
            cerr << endl;
        }
        if (totalSegments == 0U) {
            cerr << command << ": warning: failed to detect any seam" << endl;
        }
    }

    Rect2D vBB = vertexBoundingBox(contours);
    vBB.moveBy(uBB.upperLeft()); // move vBB to be root-relative

    // Make sure that vBB is bigger than iBB by one pixel in each
    // direction.  This will create a max-cost border to keep the seam
    // line from leaving the intersection region.
    Rect2D iBBPlus = iBB;
    iBBPlus.addBorder(1);
    vBB |= iBBPlus;

    // Vertex-Union bounding box: portion of uBB inside vBB
    const Rect2D uvBB = vBB & uBB;

    // Offset between vBB and uvBB
    const Diff2D uvBBOffset = uvBB.upperLeft() - vBB.upperLeft();

    Size2D mismatchImageSize;
    int mismatchImageStride;
    Diff2D uvBBStrideOffset;

    if (CoarseMask) {
        // Prepare to stride by two over uvBB to create cost image.
        // Push ul corner of vBB so that there is an even number of
        // pixels between vBB and uvBB.
        if (uvBBOffset.x % 2) {
            vBB.setUpperLeft(vBB.upperLeft() + Diff2D(-1, 0));
        }
        if (uvBBOffset.y % 2) {
            vBB.setUpperLeft(vBB.upperLeft() + Diff2D(0, -1));
        }
        uvBBStrideOffset = (uvBB.upperLeft() - vBB.upperLeft()) / 2;
        mismatchImageStride = 2;
        mismatchImageSize = (vBB.size() + Diff2D(1, 1)) / 2;
    } else {
        uvBBStrideOffset = uvBBOffset;
        mismatchImageStride = 1;
        mismatchImageSize = vBB.size();
    }

    typedef UInt8 MismatchImagePixelType;
    typedef BasicImage<MismatchImagePixelType> MismatchImageType;
    typedef BasicImage<RGBValue<MismatchImagePixelType> > VisualizeImageType;
    MismatchImageType mismatchImage(mismatchImageSize,
                                    NumericTraits<MismatchImagePixelType>::max());

    // Visualization of optimization output
    VisualizeImageType* visualizeImage = NULL;
    if (VisualizeSeam) {
        visualizeImage = new VisualizeImageType(mismatchImageSize);
    }

    // mem usage after: Visualize && CoarseMask: iBB * UInt8
    //                  Visualize && !CoarseMask: 2 * iBB * UInt8
    //                  !Visualize && CoarseMask: 1/2 * iBB * UInt8
    //                  !Visualize && !CoarseMask: iBB * UInt8

    // Calculate mismatch image
    combineTwoImagesMP(stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImageRange(*white))),
                       stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImage(*black))),
                       destIter(mismatchImage.upperLeft() + uvBBStrideOffset),
                       PixelDifferenceFunctor<ImagePixelType, MismatchImagePixelType>());

    if (DifferenceBlurRadius > 0.0) {
        gaussianSmoothing(srcImageRange(mismatchImage),
                          destImage(mismatchImage),
                          DifferenceBlurRadius);
    }

    if (visualizeImage) {
        // Dump cost image into visualize image.
        copyImage(srcImageRange(mismatchImage), destImage(*visualizeImage));

        const Diff2D offset = Diff2D(vBB.upperLeft()) - Diff2D(uBB.upperLeft());
        // Draw the initial seam line as a reference.
        for (ContourVector::const_iterator v = contours.begin(); v != contours.end(); ++v) {
            for (Contour::const_iterator c = (*v)->begin(); c != (*v)->end(); ++c) {
                for (Segment::const_iterator s = (*c)->begin(); s != (*c)->end(); ++s) {
                    Segment::const_iterator next = s;
                    ++next;
                    if (next != (*c)->end()) {
                        visualizeLine(*visualizeImage,
                                      (s->second - offset) / mismatchImageStride,
                                      (next->second - offset) / mismatchImageStride,
                                      VISUALIZE_INITIAL_PATH);
                    }
                }
            }
        }
    }

    // Areas other than intersection region have maximum cost.
    combineThreeImagesMP(stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImageRange(*whiteAlpha))),
                         stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImage(*blackAlpha))),
                         srcIter(mismatchImage.upperLeft() + uvBBStrideOffset),
                         destIter(mismatchImage.upperLeft() + uvBBStrideOffset),
                         ifThenElse(Arg1() & Arg2(),
                                    Arg3(),
                                    Param(NumericTraits<MismatchImagePixelType>::max())));

    // Strategy 1: Use GDA to optimize placement of snake vertices
    int segmentNumber;
    for (ContourVector::iterator currentContour = contours.begin();
         currentContour != contours.end();
         ++currentContour) {
        segmentNumber = 0;
        for (Contour::iterator currentSegment = (*currentContour)->begin();
             currentSegment != (*currentContour)->end();
             ++currentSegment, ++segmentNumber) {
            Segment* snake = *currentSegment;

            if (Verbose >= VERBOSE_MASK_MESSAGES) {
                cerr << command
                     << ": info: strategy 1, s"
                     << segmentNumber << ":";
                cerr.flush();
            }

            if (snake->empty()) {
                cerr << endl
                     << command
                     << ": warning: seam s"
                     << segmentNumber - 1
                     << " is a tiny closed contour and was removed before optimization"
                     << endl;
                continue;
            }

            // Move snake points to mismatchImage-relative coordinates
            for (Segment::iterator vertexIterator = snake->begin();
                 vertexIterator != snake->end();
                 ++vertexIterator) {
                vertexIterator->second =
                    (vertexIterator->second + uBB.upperLeft() - vBB.upperLeft())
                    / mismatchImageStride;
            }

            annealSnake(&mismatchImage, OptimizerWeights,
                        snake, visualizeImage);

            // Post-process annealed vertices
            Segment::iterator lastVertex = snake->previous(snake->end());
            for (Segment::iterator vertexIterator = snake->begin();
                 vertexIterator != snake->end();) {
                if (vertexIterator->first
                    && mismatchImage[vertexIterator->second] == NumericTraits<MismatchImagePixelType>::max()) {
                    // Vertex is still in max-cost region. Delete it.
                    if (vertexIterator == snake->begin()) {
                        snake->pop_front();
                        vertexIterator = snake->begin();
                    } else {
                        vertexIterator = snake->erase_after(lastVertex);
                    }

                    bool needsBreak = false;
                    if (vertexIterator == snake->end()) {
                        vertexIterator = snake->begin();
                        needsBreak = true;
                    }

                    // vertexIterator now points to next entry.

                    // It is conceivable but very unlikely that every vertex in a closed contour
                    // ended up in the max-cost region after annealing.
                    if (snake->empty()) {
                        break;
                    }

                    if (!(lastVertex->first || vertexIterator->first)) {
                        // We deleted an entire range of moveable points between two nonmoveable points.
                        // insert dummy point after lastVertex so dijkstra can work over this range.
                        if (vertexIterator == snake->begin()) {
                            snake->push_front(make_pair(true, vertexIterator->second));
                            lastVertex = snake->begin();
                        } else {
                            lastVertex = snake->insert_after(lastVertex,
                                                             make_pair(true, vertexIterator->second));
                        }
                    }

                    if (needsBreak) {
                        break;
                    }
                }
                else {
                    lastVertex = vertexIterator;
                    ++vertexIterator;
                }
            }

            if (Verbose >= VERBOSE_MASK_MESSAGES) {
                cerr << endl;
            }

            // Print an explanation if every vertex in a closed contour ended up in the
            // max-cost region after annealing.
            // FIXME explain how to fix this problem in the error message!
            if (snake->empty()) {
                cerr << endl
                     << command
                     << ": seam s"
                     << segmentNumber - 1
                     << " is a tiny closed contour and was removed after optimization"
                     << endl;
            }
        }
    }

    if (Verbose >= VERBOSE_MASK_MESSAGES) {
        cerr << command
             << ": info: strategy 2:";
        cerr.flush();
    }

    // Adjust cost image for the shortest path algorithm.
    // Areas outside the union region have epsilon cost.
    combineThreeImagesMP(stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImageRange(*whiteAlpha))),
                         stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImage(*blackAlpha))),
                         srcIter(mismatchImage.upperLeft() + uvBBStrideOffset),
                         destIter(mismatchImage.upperLeft() + uvBBStrideOffset),
                         ifThenElse(!(Arg1() || Arg2()),
                                    Param(NumericTraits<MismatchImagePixelType>::one()),
                                    Arg3()));

    if (visualizeImage) {
        combineThreeImagesMP(stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImageRange(*whiteAlpha))),
                             stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImage(*blackAlpha))),
                             srcIter(visualizeImage->upperLeft() + uvBBStrideOffset),
                             destIter(visualizeImage->upperLeft() + uvBBStrideOffset),
                             ifThenElse(Arg1() ^ Arg2(),
                                        Param(VISUALIZE_NO_OVERLAP_VALUE),
                                        Arg3()));
    }

    Rect2D withinMismatchImage(mismatchImageSize);

    // Use Dijkstra to route between moveable snake vertices over mismatchImage.
    for (ContourVector::iterator currentContour = contours.begin();
         currentContour != contours.end();
         ++currentContour) {
        segmentNumber = 0;
        for (Contour::iterator currentSegment = (*currentContour)->begin();
             currentSegment != (*currentContour)->end();
             ++currentSegment, ++segmentNumber) {
            Segment* snake = *currentSegment;

            if (snake->empty()) {
                continue;
            }

            if (Verbose >= VERBOSE_MASK_MESSAGES) {
                cerr << " s" << segmentNumber;
                cerr.flush();
            }

            for (Segment::iterator currentVertex = snake->begin(); ; ) {
                Segment::iterator nextVertex = currentVertex;
                ++nextVertex;
                if (nextVertex == snake->end()) {
                    nextVertex = snake->begin();
                }

                if (currentVertex->first || nextVertex->first) {
                    // Find shortest path between these points
                    Point2D currentPoint = currentVertex->second;
                    Point2D nextPoint = nextVertex->second;

                    Rect2D pointSurround(currentPoint, Size2D(1, 1));
                    pointSurround |= Rect2D(nextPoint, Size2D(1, 1));
                    pointSurround.addBorder(DijkstraRadius);
                    pointSurround &= withinMismatchImage;

                    // Make BasicImage to hold pointSurround portion of mismatchImage.
                    // min cost path needs inexpensive random access to cost image.
                    BasicImage<MismatchImagePixelType> mismatchROIImage(pointSurround.size());
                    copyImage(pointSurround.apply(srcImageRange(mismatchImage)),
                              destImage(mismatchROIImage));

                    vector<Point2D>* shortPath =
                        minCostPath(srcImageRange(mismatchROIImage),
                                    Point2D(nextPoint - pointSurround.upperLeft()),
                                    Point2D(currentPoint - pointSurround.upperLeft()));

                    for (vector<Point2D>::iterator shortPathPoint = shortPath->begin();
                         shortPathPoint != shortPath->end();
                         ++shortPathPoint) {
                        snake->insert_after(currentVertex,
                                            make_pair(false, *shortPathPoint + pointSurround.upperLeft()));

                        if (visualizeImage) {
                            (*visualizeImage)[*shortPathPoint + pointSurround.upperLeft()] =
                                VISUALIZE_SHORT_PATH_VALUE;
                        }
                    }

                    delete shortPath;

                    if (visualizeImage) {
                        (*visualizeImage)[currentPoint] =
                            currentVertex->first ?
                            VISUALIZE_FIRST_VERTEX_VALUE :
                            VISUALIZE_NEXT_VERTEX_VALUE;
                        (*visualizeImage)[nextPoint] =
                            nextVertex->first ?
                            VISUALIZE_FIRST_VERTEX_VALUE :
                            VISUALIZE_NEXT_VERTEX_VALUE;
                    }
                }

                currentVertex = nextVertex;
                if (nextVertex == snake->begin()) {
                    break;
                }
            }

            // Move snake vertices from mismatchImage-relative
            // coordinates to uBB-relative coordinates.
            for (Segment::iterator currentVertex = snake->begin();
                 currentVertex != snake->end();
                 ++currentVertex) {
                currentVertex->second =
                    currentVertex->second * mismatchImageStride
                    + vBB.upperLeft() - uBB.upperLeft();
            }
        }
    }

    if (Verbose >= VERBOSE_MASK_MESSAGES) {
        cerr << endl;
    }

    if (visualizeImage) {
        const std::string visualizeFilename =
            enblend::expandFilenameTemplate(VisualizeTemplate,
                                            numberOfImages,
                                            *inputFileNameIterator,
                                            OutputFileName,
                                            m);
        if (visualizeFilename == *inputFileNameIterator) {
            cerr << command
                 << ": will not overwrite input image \""
                 << *inputFileNameIterator
                 << "\" with seam-visualization image"
                 << endl;
            exit(1);
        } else if (visualizeFilename == OutputFileName) {
            cerr << command
                 << ": will not overwrite output image \""
                 << OutputFileName
                 << "\" with seam-visualization image"
                 << endl;
            exit(1);
        } else {
            if (Verbose >= VERBOSE_MASK_MESSAGES) {
                cerr << command
                     << ": info: saving seam visualization \""
                     << visualizeFilename << "\"" << endl;
            }
            ImageExportInfo visualizeInfo(visualizeFilename.c_str());
            visualizeInfo.setCompression(MASK_COMPRESSION);
            exportImage(srcImageRange(*visualizeImage), visualizeInfo);
        }

        delete visualizeImage;
    }

#ifdef DEBUG_SEAM_LINE
    std::cout << "+ createMask: contours of final optimized mask\n";
    dump_contourvector(contours, "+ createMask: ");
#endif

    // Fill contours to get final optimized mask.
    MaskType* mask = new MaskType(uBB.size());
    std::for_each(contours.begin(), contours.end(),
                  bind(fillContour<MaskType>, mask, *_1, Diff2D(0, 0)));

    // Clean up contours
    std::for_each(contours.begin(), contours.end(),
                  bind(boost::lambda::ll::for_each(), bind(call_begin(), (*_1)), bind(call_end(), (*_1)),
                       protect(bind(delete_ptr(), _1))));

    std::for_each(contours.begin(), contours.end(), bind(delete_ptr(), _1));

    return mask;
}

} // namespace enblend

#endif /* __MASK_H__ */

// Local Variables:
// mode: c++
// End:
