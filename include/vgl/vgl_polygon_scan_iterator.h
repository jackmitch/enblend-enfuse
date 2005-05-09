// This is core/vgl/vgl_polygon_scan_iterator.h
#ifndef vgl_polygon_scan_iterator_h
#define vgl_polygon_scan_iterator_h
//:
// \file
// \author Adapted from FillPolygon by J.L. Mundy
//
// \verbatim
//  Modifications
//   Binary IO added and documentation tidied up NPC, 20/03/01
//   Feb.2002 - Peter Vanroose - brief doxygen comment placed on single line
//   Nov.2003 - Peter Vanroose - made vgl_polygon_scan_iterator a templated class
//   Apr.2004 - Peter Vanroose - corrected an earlier with_boundary fix in next()
// \endverbatim

#include <vector>
#include "vigra/diff2d.hxx"
#include "vigra/iteratortags.hxx"
#include "vigra_ext/ROI.h"

using std::vector;
using vigra::Diff2D;
using vigra::Point2D;

typedef vigra_ext::ROI<Diff2D> vgl_roi;
typedef vector<Point2D> vgl_polygon_sheet;
typedef vector< vgl_polygon_sheet > vgl_polygon;

//: Fill a polygonal face with interior scan lines
//  This class provides an iterator-style interface to polygon scan
//  conversion.  There are convenient constructors from vgl_polygon, and_
//  lists of floats.  An auxiliary clipping window can be specified by the
//  constructor argument, vgl_box_2d<T> win.
//
// Concave Polygon Scan Conversion
// by Paul Heckbert
// from "Graphics Gems", Academic Press, 1990
//
// Scan convert nvert-sided concave non-simple polygon
// with vertices at  (point[i].x, point[i].y)
// for i in [0..nvert-1] within the window win by
// calling spanproc for each visible span of pixels.
// Polygon can be clockwise or counterclockwise.
// Algorithm does uniform point sampling at pixel centers.
// Inside-outside test done by Jordan's rule: a point is considered inside if
// an emanating ray intersects the polygon an odd number of times.
//
// Note: The span limits, startx and endx,  are closed intervals.
// That is, you can use the endpoints of the span as valid interior points.
// Also, the initial and final y scan lines returned by the iterator
// are interior to the polygon.  The constructor argument, win, is a clipping
// window that is intersected with the polygonal region to determine the actual
// scanned area.
//
// Example usage:
// \code
//  vgl_polygon_scan_iterator<float> psi(mypoints);
//  psi.set_include_boundary(true); // optional flag, default is true
//  for (psi.reset(); psi.next(); ) {
//    int y = psi.scany();
//    for (int x = psi.startx(); x <= psi.endx(); ++x)
//         ....
//  }
// \endcode

template <class IMAGE_ITERATOR>
class vgl_polygon_scan_iterator : private IMAGE_ITERATOR
{
public:
    typedef typename IMAGE_ITERATOR::value_type value_type;
    typedef typename IMAGE_ITERATOR::reference reference;
    typedef typename IMAGE_ITERATOR::pointer pointer;
    typedef forward_circulator_tag iterator_category;

private:
    int boundp;       //!< boolean indicating if boundary should be included or not
    int xl;           //!< left bound of current span
    double fxl;            //!< left bound of current span (floating point value)
    int xr;           //!< right bound of current span
    double fxr;            //!< right bound of current span (floating point value)
    int k;            //!< current index of vertices ordered by increasing y
    int y0;           //!< bottommost scan line
    int y1;           //!< topmost scan line
    int y;            //!< current scan line
    double fy;             //!< floating point value of current scan line (i.e. T(y))
    int curcrossedge; //!< crossedge marking start of next scan segment
    vgl_roi win;       //!< clipping window
    bool have_window;

    vgl_polygon poly_; //!< the polygon

    int vigra_x;
    int vigra_y;

public:
    // Constructors/Destructor---------------------------------------------------

    //: Construct with a polygon and bool indicating whether boundary included
    vgl_polygon_scan_iterator(const IMAGE_ITERATOR &iterator,
            vgl_polygon const& face,
            bool boundaryp = true);

    //: Construct with a polygon, bool indicating whether boundary included and window (area visible)
    vgl_polygon_scan_iterator(const IMAGE_ITERATOR &iterator,
            vgl_polygon const& face,
            bool boundaryp,
            vgl_roi const& window);

    //: Destructor
    ~vgl_polygon_scan_iterator();

    //Functions----------------------------------------------------------

    // VIGRA image iterator adaptor functions
    // Move forward pre-increment
    vgl_polygon_scan_iterator & operator++() {
        int nextX = vigra_x + 1;
        int nextY = vigra_y;
        if (nextX > endx()) {
            //int beforeNextY = nextY;
            if (!next()) {
                reset();
                next();
            }
            nextY = scany();
            //if (nextY != (beforeNextY + 1)) cout << "called next() ybefore=" << beforeNextY << " yafter=" << nextY << endl;
            nextX = startx();
        }
        //cout << "vigra (" << vigra_x << ", " << vigra_y << ") "
        //     << "next (" << nextX << ", " << nextY << ")" << endl;
        Diff2D delta(nextX - vigra_x, nextY - vigra_y);
        vigra_x = nextX;
        vigra_y = nextY;
        IMAGE_ITERATOR::operator+=(delta);
        return *this;
    }

    // Move forward post-increment
    vgl_polygon_scan_iterator operator++(int) {
        vgl_polygon_scan_iterator ret(*this);
        ++(*this);
        return ret;
    }

    // equality
    bool operator==(vgl_polygon_scan_iterator const & o) const {
        return IMAGE_ITERATOR::operator==(o);
    }

    // inequality
    bool operator!=(vgl_polygon_scan_iterator const & o) const {
        return IMAGE_ITERATOR::operator!=(o);
    }

    reference operator*() const {
        return IMAGE_ITERATOR::operator*();
    }

    pointer operator->() const {
        return IMAGE_ITERATOR::operator->();
    }

private:
    //: Resets iterator to first segment of first scan line
    void reset();

    //: Moves iterator to next segment
    bool next();

    //: Returns current scan line
    inline int scany() const { return y-1; }

    //: Returns start of current span
    inline int startx() const { return xl; }

    //: Returns end of current span
    inline int endx() const { return xr; }

    //: Returns start of current span (floating point value)
    inline double fstartx() const { return fxl; }

    //: Returns end of current span (floating point value)
    inline double fendx() const { return fxr; }

    //: Returns current scan line (floating point value)
    inline double fscany() const { return fy; }

    //: Vertex index - uniquely identifies a vertex in the array chains
    struct vertind {
        unsigned int chainnum; //!< which chain the vertex is part of
        unsigned int vertnum;  //!< which vertex in the chain
    };

    //: Describes an edge crossing the current scan line
    struct crossedge {
        double x;       //!< x coord of edge's intersection with current scanline
        double dx;      //!< change in x with respect to y
        vertind v; //!< edge goes from vertex v.vertnum to v.vertnum + 1
    };

// Internals ---------------------------------------------------------------

private:

    vertind * yverts;       //!< array of all vertices ordered by y coordinate
    crossedge * crossedges; //!< array of edges crossing current scan line
    int numcrossedges;      //!< number of edges currently crossing scan line
    int numverts;           //!< total number of vertices comprising face

    // Returns x coord of vertex v
    inline double get_x(vertind v) const { return (double)poly_[v.chainnum][v.vertnum].x; }

    // Returns y coord of vertex v
    inline double get_y(vertind v) const { return (double)poly_[v.chainnum][v.vertnum].y; }

    // Returns vertex v
    inline Point2D get_pt( vertind v ) { return poly_[v.chainnum][v.vertnum]; }

    // assumes poly_, win, have_window, boundp are set
    void init();

    // Deletes edge (v,get_next_vert(v)) from crossedges array
    void delete_edge( vertind v );

    // Inserts edge (v,get_next_vert(v)) into crossedges array
    void insert_edge( vertind v );

    // Returns next vertex on chain
    void get_next_vert( vertind v, vertind & next );

    // Returns prev vertex on chain
    void get_prev_vert( vertind v, vertind & prev );

    // For debugging purposes
    void display_chains();
    void display_crossedges();

    inline int irnd(double x) const { return (int)floor(x + 0.5); };
};

#include "vgl/vgl_polygon_scan_iterator.txx"

#endif // vgl_polygon_scan_iterator_h
