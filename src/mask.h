/*
 * Copyright (C) 2004-2011 Andrew Mihal
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
#include "postoptimizer.h"
#include "graphcut.h"
#include "maskcommon.h"

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

void
dump_segment(const Segment& segment,
             const std::string& prefix = "", std::ostream& out = std::cout)
{
    const unsigned points_per_line = 5U;
    unsigned n = 1U;
    const size_t size = segment.size();

    out <<
        prefix << "{segment with " << size <<
        " point(s): // suffix 'f' means frozen point\n" << prefix;
    for (Segment::const_iterator i = segment.begin(); i != segment.end(); ++i, ++n)
    {
        // Append an 'f' to non-movable ("frozen") segment points.
        out << ' ' << i->second << (i->first ? ' ' : 'f');
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


template <typename ImageType>
void
visualizePoint(ImageType& image,
               const Point2D& location,
               typename ImageType::PixelType value)
{
    typedef typename ImageType::Iterator Iterator;

    Iterator point(image.upperLeft() + location);

    if (point.x >= image.upperLeft().x && point.x < image.lowerRight().x &&
        point.y >= image.upperLeft().y && point.y < image.lowerRight().y) {
        image.accessor().set(value, point);
    }
}


template <typename T>
inline T
round_to(double x)
{
    return static_cast<T>(x >= 0.0 ? x + 0.5 : x - 0.5);
}


template <typename ImageType>
void
visualizePoint(ImageType& image,
               const Point2D& location,
               typename ImageType::PixelType value,
               marker_t marker,
               int radius)
{
    if (radius <= 0)
    {
        return;
    }

    const int r_sqrt2 = round_to<int>(static_cast<double>(radius) / 1.414213562373095);

    switch (marker)
    {
    case NO_MARKER: return;

    case DOT_MARKER:
        visualizePoint(image, location, value);
        break;

    case PLUS_MARKER:
        for (int i = -radius; i <= radius; ++i)
        {
            visualizePoint(image, location + Point2D(i, 0), value);
            visualizePoint(image, location + Point2D(0, i), value);
        }
        break;

    case CROSS_MARKER:
        for (int i = -r_sqrt2; i <= r_sqrt2; ++i)
        {
            visualizePoint(image, location + Point2D(i, i), value);
            visualizePoint(image, location + Point2D(-i, i), value);
        }
        break;

    case HOLLOW_SQUARE_MARKER:
        for (int i = -r_sqrt2; i <= r_sqrt2; ++i)
        {
            visualizePoint(image, location + Point2D(-r_sqrt2, i), value);
            visualizePoint(image, location + Point2D(r_sqrt2, i), value);
            visualizePoint(image, location + Point2D(i, -r_sqrt2), value);
            visualizePoint(image, location + Point2D(i, r_sqrt2), value);
        }
        break;

    case HOLLOW_DIAMOND_MARKER:
        for (int i = 0; i <= radius; ++i)
        {
            visualizePoint(image, location + Point2D(i - radius, i), value);
            visualizePoint(image, location + Point2D(i, i - radius), value);
            visualizePoint(image, location + Point2D(-i + radius, i), value);
            visualizePoint(image, location + Point2D(-i, i - radius), value);
        }
        break;

    default:
        visualizePoint(image, location, value);
    }
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
                       int nftStride, MaskType* nftOutputImage, int vectorizeDistance = 0)
{
    typedef typename MaskType::PixelType MaskPixelType;
    typedef typename MaskType::traverser MaskIteratorType;

    const double diagonalLength =
        hypot(static_cast<double>(nftOutputImage->width()),
              static_cast<double>(nftOutputImage->height()));
    if(vectorizeDistance == 0)
        vectorizeDistance =
            MaskVectorizeDistance.is_percentage() ?
            static_cast<int>(ceil(MaskVectorizeDistance.value() / 100.0 * diagonalLength)) :
            MaskVectorizeDistance.value();
    
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
                     FileNameList::const_iterator inputFileNameIterator,
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
            const bool too_small = maskInfo.width() < uBB.width() || maskInfo.height() < uBB.height();

            std::string category;
            if (too_small) {
                category = "warning: ";
            }

            cerr << command
                 << ": " << category << "mask in \"" << maskFilename << "\" has size "
                 << "(" << maskInfo.width() << "x" << maskInfo.height() << "),\n"
                 << command
                 << ": " << category << "    but image union has size " << uBB.size() << ";\n"
                 << command
                 << ": " << category << "    make sure this is the right mask for the given images"
                 << endl;

            if (!too_small) {
                // Mask is too large, loading it would cause a segmentation fault.
                exit(1);
            }
        }
        importImage(maskInfo, destImage(*mask));
        return mask;
    }

    // Start by using the nearest feature transform to generate a mask.
    Size2D mainInputSize, mainInputBBSize;
    Rect2D mainInputBB;
#define SKIP_OPTIMIZER
    bool graphCutDebug = true;
    
    int mainStride;
    if (CoarseMask && !graphCutDebug) {
        // Do NFT at 1/CoarsenessFactor scale.
        // uBB rounded up to multiple of CoarsenessFactor pixels in each direction
        mainInputSize = Size2D((uBB.width() + CoarsenessFactor - 1) / CoarsenessFactor,
                              (uBB.height() + CoarsenessFactor - 1) / CoarsenessFactor);
        mainInputBBSize = Size2D((iBB.width() + CoarsenessFactor - 1) / CoarsenessFactor,
                              (iBB.height() + CoarsenessFactor - 1) / CoarsenessFactor);
        //mainInputBB = Rect2D(Size2D(uBB.width() / CoarsenessFactor,
          //                         uBB.height() / CoarsenessFactor));
        mainInputBB = Rect2D(Point2D((iBB.upperLeft().x - uBB.upperLeft().x + CoarsenessFactor - 1)/ CoarsenessFactor, 
                            (iBB.upperLeft().y - uBB.upperLeft().y + CoarsenessFactor - 1)/CoarsenessFactor),
                            mainInputBBSize);
        mainStride = CoarsenessFactor;
    } else {
        // Do NFT at 1/1 scale.
        mainInputSize = uBB.size();
        mainInputBB = Rect2D(iBB);
        if(mainInputBB.upperLeft().x >= uBB.upperLeft().x)
                mainInputBB.moveBy(-uBB.upperLeft());
        else mainInputBB.moveBy(uBB.upperLeft());
        mainStride = 1;
    }

    Size2D mainOutputSize;
    Diff2D mainOutputOffset;
    if (!CoarseMask && !OptimizeMask) {
        // We are not going to vectorize the mask.
        mainOutputSize = mainInputSize;
        mainOutputOffset = Diff2D(0, 0);
    } else {
        // Add 1-pixel border all around the image for the vectorization algorithm.
        mainOutputSize = mainInputSize + Diff2D(2, 2);
        mainOutputOffset = Diff2D(1, 1);
    }

    // mem usage before: 0
    // mem usage after: CoarseMask: 1/8 * uBB * MaskType
    //                  !CoarseMask: uBB * MaskType
    MaskType* mainOutputImage = new MaskType(mainOutputSize);
    
    #ifdef GRAPHCUT_DBG
        cout<<iBB<<endl<<mainInputBB<<endl;
#endif

    graphCut(stride(mainStride, mainStride, iBB.apply(srcImageRange(*white))),
                            stride(mainStride, mainStride, iBB.apply(srcImage(*black))),
                            destIter(mainOutputImage->upperLeft() + mainOutputOffset),
                            stride(mainStride, mainStride, uBB.apply(srcImageRange(*whiteAlpha))),
                            stride(mainStride, mainStride, uBB.apply(srcImage(*blackAlpha))),
                            ManhattanDistance,
                            wraparound ? HorizontalStrip : OpenBoundaries, 
                            mainInputBB);
    
 
    
#ifdef GRAPHCUT_DBG
        exportImage(srcImageRange(*mainOutputImage), ImageExportInfo("./debug/output.tif").setPixelType("UINT8"));
#endif

    /*nearestFeatureTransform(stride(mainStride, mainStride, uBB.apply(srcImageRange(*whiteAlpha))),
                            stride(mainStride, mainStride, uBB.apply(srcImage(*blackAlpha))),
                            destIter(mainOutputImage->upperLeft() + mainOutputOffset),
                            ManhattanDistance,
                            wraparound ? HorizontalStrip : OpenBoundaries);*/

#ifdef DEBUG_NEAREST_FEATURE_TRANSFORM
    {
        typedef pair<const char*, const MaskType*> ImagePair;

        const ImagePair nft[] = {
            std::make_pair("blackmask", blackAlpha),
            std::make_pair("whitemask", whiteAlpha),
            //std::make_pair("nft-input", nftInputImage),
            std::make_pair("nft-output", mainOutputImage)
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
        return mainOutputImage;
    }

    // Vectorize the seam lines found in nftOutputImage.
    Contour rawSegments;
    if(graphCutDebug)
    vectorizeSeamLine(rawSegments,
                      whiteAlpha, blackAlpha,
                      uBB,
                      mainStride, mainOutputImage, 4);
    else 
        vectorizeSeamLine(rawSegments,
                      whiteAlpha, blackAlpha,
                      uBB,
                      mainStride, mainOutputImage);
    delete mainOutputImage;

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

        // Color the parts of the visualize image where the two images
        // to be blended do not overlap.
        combineThreeImagesMP(stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImageRange(*whiteAlpha))),
                             stride(mismatchImageStride, mismatchImageStride, uvBB.apply(srcImage(*blackAlpha))),
                             srcIter(visualizeImage->upperLeft() + uvBBStrideOffset),
                             destIter(visualizeImage->upperLeft() + uvBBStrideOffset),
                             ifThenElse(Arg1() & Arg2(),
                                        Arg3(),
                                        Param(VISUALIZE_NO_OVERLAP_VALUE)));

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
                    visualizePoint(*visualizeImage,
                                   (s->second - offset) / mismatchImageStride,
                                   s->first ? VISUALIZE_MOVABLE_POINT : VISUALIZE_FROZEN_POINT,
                                   s->first ? MARK_MOVABLE_POINT : MARK_FROZEN_POINT,
                                   2);
                }
            }
        }
    }

#ifndef SKIP_OPTIMIZER
    
    int segmentNumber;
    // Move snake points to mismatchImage-relative coordinates
    for (ContourVector::iterator currentContour = contours.begin();
             currentContour != contours.end();
             ++currentContour) {
            segmentNumber = 0;
            for (Contour::iterator currentSegment = (*currentContour)->begin();
                 currentSegment != (*currentContour)->end();
                 ++currentSegment, ++segmentNumber) {
                Segment* snake = *currentSegment;
                for (Segment::iterator vertexIterator = snake->begin();
                     vertexIterator != snake->end();
                     ++vertexIterator) {
                    vertexIterator->second =
                        (vertexIterator->second + uBB.upperLeft() - vBB.upperLeft()) /
                        mismatchImageStride;
                }
            }
    }
    
    vector<double> *params = new(vector<double>);
    
    OptimizerChain<MismatchImagePixelType, MismatchImageType, VisualizeImageType, AlphaType>
        *defaultOptimizerChain = new OptimizerChain<MismatchImagePixelType, MismatchImageType, VisualizeImageType, AlphaType>
                              (&mismatchImage, visualizeImage,
                              &mismatchImageSize, &mismatchImageStride,
                              &uvBBStrideOffset, &contours, &uBB, &vBB, params,
                              whiteAlpha, blackAlpha, &uvBB);
    
        // Add Strategy 1: Use GDA to optimize placement of snake vertices
    defaultOptimizerChain->addOptimizer("anneal");
    
        // Add Strategy 2: Use Dijkstra shortest path algorithm between snake vertices
    defaultOptimizerChain->addOptimizer("dijkstra");
   
        // Fire optimizer chain (runs every optimizer on the list in sequence)
    defaultOptimizerChain->runOptimizerChain();
    
    // Move snake vertices from mismatchImage-relative
    // coordinates to uBB-relative coordinates.    
    for (ContourVector::iterator currentContour = contours.begin();
             currentContour != contours.end();
             ++currentContour) {
            segmentNumber = 0;
            for (Contour::iterator currentSegment = (*currentContour)->begin();
                 currentSegment != (*currentContour)->end();
                 ++currentSegment, ++segmentNumber) {
                Segment* snake = *currentSegment;
                for (Segment::iterator currentVertex = snake->begin();
                     currentVertex != snake->end();
                     ++currentVertex) {
                    currentVertex->second =
                        currentVertex->second * mismatchImageStride +
                        vBB.upperLeft() - uBB.upperLeft();
                }
            }
    }
    
    delete params;
    delete defaultOptimizerChain;
    
#endif // !SKIP_OPTIMIZER

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
