/*
 *  Copyright (C) 2004 Andrew Mihal
 *
 *  This software is an extension of the VIGRA computer vision library.
 *  ( Version 1.2.0, Aug 07 2003 )
 *  You may use, modify, and distribute this software according
 *  to the terms stated in the LICENSE file included in
 *  the VIGRA distribution.
 *
 *  VIGRA is Copyright 1998-2002 by Ullrich Koethe
 *  Cognitive Systems Group, University of Hamburg, Germany
 *
 *  The VIGRA Website is
 *      http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/
 *  Please direct questions, bug reports, and contributions to
 *      koethe@informatik.uni-hamburg.de
 *
 *  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR
 *  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */
#ifndef VIGRA_EXT_CACHEDFILEIMAGE_HXX
#define VIGRA_EXT_CACHEDFILEIMAGE_HXX

#define CACHED_FILE_IMAGE_BLOCKSIZE 1000000
#define CACHED_FILE_IMAGE_CACHELINES 10

#include <errno.h>
#include <iostream>
#include <list>
#include <stdio.h>
#include <unistd.h>
#include <vigra/utilities.hxx>

using std::list;

namespace vigra {

template <class Iterator>
class CachedFileSequentialAccessIteratorPolicy
{
public:
    typedef Iterator BaseType;
    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::difference_type::MoveX difference_type;
    typedef typename Iterator::reference reference;
    typedef typename Iterator::index_reference index_reference;
    typedef typename Iterator::pointer pointer;
    typedef typename Iterator::iterator_category iterator_category;

    static void initialize(BaseType & d) {
        width_ = d.i->width();
        height_ = d.i->height();
    }

    static reference dereference(BaseType const & d) {
        return *d;
    }

    static index_reference dereference(BaseType const & d, difference_type n) {
        int dy = n / width_;
        int dx = n % width_;
        if (d.x + dx >= width_) {dy++; dx -= width_;};
        else if (d.x + dx < 0) {dy--; dx += width_;};
        return d(dx, dy);
    }

    static bool equal(BaseType const & d1, BaseType const & d2) {
        return d1 == d2;
    }

    static bool less(BaseType const & d1, BaseType const & d2) {
        return (d1.y*width_ + d1.x) < (d2.y*width_ + d2.x);
    }

    static difference_type difference(BaseType const & d1, BaseType const & d2) {
        return (d1.y*width_ + d1.x) - (d2.y*width_ + d2.x);
    }

    static void increment(BaseType & d) {
        ++d.x;
        if (d.x == width_) {
            d.x = 0;
            ++d.y;
        }
    }

    static void decrement(BaseType & d) {
        --d.x;
        if (d.x < 0) {
            d.x = width_ - 1;
            --d.y;
        }
    }

    static void advance(BaseType & d, difference_type n) {
        int dy = n / width_;
        int dx = n % width_;
        d.x += dx;
        d.y += dy;
        if (d.x >= width_) {++d.y; d.x -= width_;}
        if (d.x < 0) {--dy.y; d.x += width_;}
    }

private:
    int width_;
    int height_;

};

template <class IMAGEITERATOR, class PIXELTYPE, class REFERENCE, class POINTER>
class CachedFileImageIteratorBase
{
public:
    typedef CachedFileImageIteratorBase<IMAGEITERATOR,
            PIXELTYPE, REFERENCE, POINTER> self_type;
    typedef CachedFileImage<PIXELTYPE> image_type;
    typedef PIXELTYPE value_type;
    typedef PIXELTYPE PixelType;
    typedef REFERENCE reference;
    typedef REFERENCE index_reference;
    typedef POINTER pointer;
    typedef Diff2D difference_type;
    typedef image_traverser_tag iterator_category;
    typedef RowIterator<IMAGEITERATOR> row_iterator;
    typedef ColumnIterator<IMAGEITERATOR> column_iterator;
    typedef int MoveX;
    typedef int MoveY;

    MoveX x;
    MoveY y;
    image_type *i;

    IMAGEITERATOR & operator+=(difference_type const & s) {
        x += s.x;
        y += s.y;
        return static_cast<IMAGEITERATOR &>(*this);
    }

    IMAGEITERATOR & operator-=(difference_type const & s) {
        x -= s.x;
        y -= s.y;
        return static_cast<IMAGEITERATOR &>(*this);
    }

    IMAGEITERATOR operator+(difference_type const & s) const {
        IMAGEITERATOR ret(static_cast<IMAGEITERATOR const &>(*this));
        ret += s;
        return ret;
    }

    IMAGEITERATOR operator-(difference_type const & s) const {
        IMAGEITERATOR ret(static_cast<IMAGEITERATOR const &>(*this));
        ret -= s;
        return ret;
    }

    difference_type operator-(CachedFileImageIteratorBase const & rhs) const {
        return difference_type(x-rhs.x, y-rhs.y);
    }

    bool operator==(CachedFileImageIteratorBase const & rhs) const {
        return (x == rhs.x) && (y == rhs.y);
    }

    bool operator!=(CachedFileImageIteratorBase const & rhs) const {
        return (x != rhs.x) || (y != rhs.y);
    }

    reference operator*() const {
        return (*i)(x, y);
    }

    // pointer is supposed to be a weak_ptr
    pointer operator->() const {
        return (*i)[y] + x;
    }

    index_reference operator[](difference_type const & d) const {
        return (*i)(x+d.x, y+d.y);
    }

    index_reference operator()(int dx, int dy) const {
        return (*i)(x+dx, y+dy);
    }

    // pointer is supposed to be a weak_ptr
    pointer operator[](int dy) const {
        return (*i)[y + dy] + x;
    }

    row_iterator rowIterator() const {
        return row_iterator(*this);
    }

    column_iterator columnIterator() const {
        return column_iterator(*this);
    }

protected:

    CachedFileImageIteratorBase(int X, int Y, image_type *I) : x(X), y(Y), i(I) { }

    CachedFileImageIteratorBase() : x(0), y(0), i(NULL) { }

};
    
template <class PIXELTYPE>
class CachedFileImageIterator
: public CachedFileImageIteratorBase<CachedFileImageIterator<PIXELTYPE>,
                PIXELTYPE, PIXELTYPE &, PIXELTYPE *>
// FIXME this needs to be a weak_ptr    ^^^^^^^^^^^
// in case someone uses the iterator to get a pointer to cached data.
{
public:

    typedef CachedFileImageIteratorBase<CachedFileImageIterator,
            PIXELTYPE, PIXELTYPE &, PIXELTYPE *> Base;

    CachedFileImageIterator(int x, int y, CachedFileImage<PIXELTYPE> *i)
    : Base(x, y, i)
    {}

    CachedFileImageIterator()
    : Base(0, 0, NULL)
    {}

};

template <class PIXELTYPE>
class ConstCachedFileImageIterator
: public CachedFileImageIteratorBase<ConstCachedFileImageIterator<PIXELTYPE>,
                PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *>
// FIXME this needs to be a weak_ptr          ^^^^^^^^^^^^^^^^^
// in case someone uses the iterator to get a pointer to cached data.
{
public:

    typedef CachedFileImageIteratorBase<ConstCachedFileImageIterator,
            PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *> Base;
    // FIXME this needs to be a weak_ptr  ^^^^^^^^^^^^^^^^^

    ConstCachedFileImageIterator(int x, int y, CachedFileImage<PIXELTYPE> *i)
    : Base(x, y, i)
    {}

    ConstCachedFileImageIterator(CachedFileImageIterator<PIXELTYPE> const & rhs)
    : Base(rhs.x, rhs.y, rhs.i)
    {}

    ConstCachedFileImageIterator()
    : Base(0, 0, NULL)
    {}

    ConstCachedFileImageIterator &
    operator=(CachedFileImageIterator<PIXELTYPE> const & rhs)
    {
        x = rhs.x;
        y = rhs.y;
        i = rhs.i;
        return *this;
    }

};

template <class T> struct IteratorTraits;

template <class PIXELTYPE>
class CachedFileImage {
public:

    typedef PIXELTYPE value_type;
    typedef PIXELTYPE PixelType;
    typedef PIXELTYPE & reference;
    typedef PIXELTYPE const & const_reference;
    // FIXME these need to be weak_ptrs
    typedef PIXELTYPE * pointer;
    typedef PIXELTYPE const * const_pointer;
    typedef CachedFileImageIterator<PIXELTYPE> traverser;
    typedef ConstCachedFileImageIterator<PIXELTYPE> const_traverser;
    typedef IteratorAdaptor<CachedFileSequentialAccessIteratorPolicy<traverser> > iterator;
    typedef IteratorAdaptor<CachedFileSequentialAccessIteratorPolicy<const_traverser> > const_iterator;
    typedef Diff2D difference_type;
    typedef Size2D size_type;
    typedef typename IteratorTraits<traverser>::DefaultAccessor Accessor;
    typedef typename IteratorTraits<const_traverser>::DefaultAccessor ConstAccessor;

    struct Allocator {
        static value_type * allocate(int n) {
            return (value_type *)::operator new(n*sizeof(value_type));
        }
        static void deallocate(value_type * p) {
            ::operator delete(p);
        }
    };

    CachedFileImage() {
        initMembers();
    }

    CachedFileImage(int width, int height) {
        vigra_precondition((width >= 0) && (height >= 0),
                "CachedFileImage::CachedFileImage(int width, int height): "
                "width and height must be >= 0.\n");
        initMembers();
        resize(width, height, value_type());
    }

    explicit CachedFileImage(difference_type const & size) {
        vigra_precondition((size.x >= 0) && (size.y >= 0),
                "CachedFileImage::CachedIfelImage(Diff2D size): "
                "size.x and size.y must be >= 0.\n");
        initMembers();
        resize(size.x, size.y, value_type());
    }

    CachedFileImage(int width, int height, value_type const & d) {
        vigra_precondition((width >= 0) && (height >= 0),
                "CachedFileImage::CachedFileImage(int width, int height, value_type const & ): "
                "width and height must be >= 0.\n");
        initMembers();
        resize(width, height, d);
    }

    CachedFileImage(const CachedFileImage & rhs) {
        initMembers();
        resizeCopy(rhs);
    }

    ~CachedFileImage() {
        deallocate();
    }

    CachedFileImage & operator=(const CachedFileImage &rhs);
    CachedFileImage & init(value_type const & pixel);

    void resize(int width, int height) {
        resize(width, height, value_type());
    }

    void resize(difference_type const & size) {
        resize(size.x, size.y, value_type());
    }

    void resize(int width, int height, value_type const & d);
    void resizeCopy(const CachedFileImage & rhs);
    void swap( CachedFileImage<PIXELTYPE>& rhs );

    int width() const {
        return width_;
    }

    int height() const {
        return height_;
    }

    size_type size() const {
        return size_type(width(), height());
    }

    bool isInside(difference_type const & d) const {
        return d.x >= 0 && d.y >= 0 &&
               d.x < width() && d.y < height();
    }

    reference operator[](difference_type const & d) {
        return (getLinePointer(d.y))[d.x];
    }

    const_reference operator[](difference_type const & d) const {
        return (getLinePointer(d.y))[d.x];
    }

    reference operator()(int dx, int dy) {
        return (getLinePointer(dy))[dx];
    }

    const_reference operator()(int dx, int dy) const {
        return (getLinePointer(dy))[dx];
    }

    // dangerous - needs to return a weak_ptr
    pointer operator[](int dy) {
        return getLinePointer(dy);
    }

    // dangerous - needs to return a weak_ptr
    const_pointer operator[](int dy) const {
        return getLinePointer(dy);
    }

    traverser upperLeft() {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::upperLeft(): image must have non-zero size.");
        return traverser(0, 0, this);
    }

    traverser lowerRight() {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::lowerRight(): image must have non-zero size.");
        return traverser(width_, height_, this);
    }

    const_traverser upperLeft() const {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::upperLeft(): image must have non-zero size.");
        return const_traverser(0, 0, this);
    }

    const_traverser lowerRight() const {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::lowerRight(): image must have non-zero size.");
        return const_traverser(width_, height_, this);
    }

    iterator begin() {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::begin(): image must have non-zero size.");
        return iterator(traverser(0, 0, this));
    }

    iterator end() {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::end(): image must have non-zero size.");
        return iterator(traverser(0, height_, this));
    }

    const_iterator begin() const {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::begin(): image must have non-zero size.");
        return const_iterator(const_traverser(0, 0, this));
    }

    const_iterator end() const {
        vigra_precondition(width_ > 0 && height_ > 0,
                "CachedFileImage::end(): image must have non-zero size.");
        return const_iterator(const_traverser(0, height_, this));
    }

    Accessor accessor() {
        return Accessor();
    }

    ConstAccessor accessor() const {
        return ConstAccessor();
    }

private:

    void deallocate();

    void initLineStartArray();
    
    // obtain a pointer to the beginning of a line.
    PIXELTYPE * getLinePointer(int dy);

    // Free space, if necessary, by swapping out a block of lines to the file.
    void swapLeastRecentlyUsedBlock();

    // Lazy creation of tmp file for swapping image data to disk
    void initTmpfile();

    inline int lineToBlockNumber(int line) const {
        return line / linesPerBlocksize_;
    }

    inline int blockToFirstLineNumber(int block) const {
        return block * linesPerBlocksize_;
    }

    void initMembers() {
        linesPerBlocksize_ = 0;
        blocksAllowed_ = 0;
        blockLRU_ = NULL:
        lines_ = NULL:
        width_ = 0;
        height_ = 0;
        tmpFile_ = NULL:
        tmpFilename_ = NULL:
    }

    // how many image lines are loaded in one fread
    int linesPerBlocksize_;
    // how many blocks may exist in memory at one time
    int blocksAllowed_;

    // lru replacement policy
    // most recently used block is at start of list.
    // least recently used block is at end of list.
    list<int> *blockLRU_;

    PIXELTYPE ** lines_;
    int width_, height_;

    FILE *tmpFile_;
    char *tmpFilename_;
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::deallocate() {
    delete blockLRU_;
    if (lines_ != NULL) {
        // Go through lines and delete any allocated memory there.
        for (int line = 0; line < height_; line++) {
            PIXELTYPE *p = lines_[line];
            if (p != NULL) {
                for (int column = 0; column < width_; column++) {
                    (p[column]).~PIXELTYPE();
                }
                Allocator::deallocate(p);
            }
        }
        delete[] lines_;
    }
    if (tmpFile_ != NULL) {
        fclose(tmpFile_);
        unlink(tmpFilename_);
    }
    delete[] tmpFilename_;
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::initLineStartArray() {

    // Number of lines to load in one block.
    linesPerBlocksize_ = ceil(CACHED_FILE_IMAGE_BLOCKSIZE / (width_ * sizeof(PIXELTYPE)));
    // a block must be at least one line.
    linesPerBlocksize_ = max(1, linesPerBlocksize_);

    // FIXME some mechanism for determining how much of the image is going to be in memory.
    // implement a CachedFileImageDirector
    blocksAllowed_ = CACHED_FILE_IMAGE_CACHELINES;

    // Cap blocksAllowed_ at the actual number of blocks needed.
    int blocksNeeded = ceil(height_ / linesPerBlocksize_);
    blocksAllowed_ = min(blocksAllowed_, blocksNeeded);

    // Create the blockLRU list.
    blockLRU_ = new list<int>();

    lines_ = new PIXELTYPE*[height_];

    // Allocate mem for the first linesPerBlocksize_*blocksAllowed_ lines.
    int line = 0;
    for (int block = 0; block < blocksAllowed_; block++) {
        blockLRU_->push_front(block);
        for (int subblock = 0; subblock < linesPerBlocksize_; subblock++, line++) {
            if (line < height_) {
                lines_[line] = Allocator::allocate(width_);
            }
        }
    }

    // All remaining lines (if any) are null (swapped out)
    for (; line < height_; line++) {
        lines_[line] = NULL:
    }

    return;
};

template <class PIXELTYPE>
PIXELTYPE * CachedFileImage<PIXELTYPE>::getLinePointer(int dy) {
    int blockNumber = lineToBlockNumber(dy);
    PIXELTYPE *line = lines_[dy];

    // Check if line dy is swapped out.
    if (line == NULL) {
        int firstLineInBlock = blockToFirstLineNumber(blockNumber);

        // Make space for new block.
        swapLeastRecentlyUsedBlock();

        // Find the right spot in the file.
        off_t offset = width_ * firstLineInBlock * sizeof(PIXELTYPE);
        if (fseeko(tmpFile_, offset, SEEK_SET) != 0) {
            vigra_fail(strerror(errno));
        }

        // Allocate lines for new block.
        for (int l = 0; l < linesPerBlocksize_; l++) {
            int absoluteLineNumber = l + firstLineInBlock;
            if (absoluteLineNumber >= height_) break;
            lines_[absoluteLineNumber] = Allocator::allocate(width_);

            // fill the line with data from the file.
            if (fread(lines_[absoluteLineNumber], sizeof(PIXELTYPE), width_, tmpFile_)
                    < width_) {
                vigra_fail("CachedFileImage: error reading from image backing file.\n");
            }
        }

        // Mark this block as most recently used.
        blockLRU_->push_front(blockNumber);

        return lines_[dy];
    }
    else {
        // make blockNumber the least recently used block.
        blockLRU_->remove(blockNumber);
        blockLRU_->push_front(blockNumber);
        return line;
    }
};


template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::swapLeastRecentlyUsedBlock() {
    int block = blockLRU_->back();
    blockLRU_->pop_back();

    int firstLineInBlock = blockToFirstLineNumber(blockNumber);

    // Lazy init the temp file.
    if (tmpFile_ == NULL) initTmpfile();

    // Find the right spot in the file.
    off_t offset = width_ * firstLineInBlock * sizeof(PIXELTYPE);
    if (fseeko(tmpFile_, offset, SEEK_SET) != 0) {
        vigra_fail(strerror(errno));
    }

    for (int l = 0; l < linesPerBlocksize_; l++) {
        int absoluteLineNumber = l + firstLineInBlock;
        if (absoluteLineNumber >= height_) break;
        if (fwrite(lines_[absoluteLineNumber], sizeof(PIXELTYPE), width_, tmpFile_)
                != width_) {
            vigra_fail("CachedFileImage: error writing to image backing file.\n");
        }
        // Deallocate line
        PIXELTYPE *p = lines_[absoluteLineNumber];
        for (int column = 0; column <= width_; column++) {
            (p[column]).~PIXELTYPE();
        }
        Allocator::deallocate(p);
        lines_[absoluteLineNumber] = NULL;
    }
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::initTmpfile() {
    char filenameTemplate[] = ".enblend_tmpXXXXXX";

    int tmpFD = mkstemp(filenameTemplate);
    if (tmpFD < 0) {
        vigra_fail("CachedFileImage: unable to create image backing file.\n");
    }

    int filenameTemplateLength = strlen(filenameTemplate) + 1;
    tmpFilename_ = new char[filenameTemplateLength];
    strncpy(tmpFilename_, filenameTemplate, filenameTemplateLength);

    tmpFile_ = fdopen(tmpFD, "wb+");
    if (tmpFile_ == NULL) {
        vigra_fail(strerror(errno));
    }
};

template <class PIXELTYPE>
CachedFileImage<PIXELTYPE> &
CachedFileImage<PIXELTYPE>::operator=(const CachedFileImage<PIXELTYPE> & rhs) {
    if (this != &rhs) {
        if ((width() != rhs.width()) ||
                (height() != rhs.height())) {
            resizeCopy(rhs);
        }
        else {
            const_iterator is = rhs.begin();
            const_iterator iend = rhs.end();
            iterator id = begin();
            for(; is != iend; ++is, ++id) *id = *is;
        }
    }
    return *this;
};

template <class PIXELTYPE>
CachedFileImage<PIXELTYPE> &
CachedFileImage<PIXELTYPE>::init(value_type const & pixel) {
    iterator i = begin();
    iterator iend = end();

    for(; i != iend; ++i) *i = pixel;

    return *this;
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::resize(int width, int height, value_type const & d) {
    vigra_precondition((width >= 0) && (height >= 0),
            "CachedFileImage::resize(int width, int height, value_type const &): "
            "width and height must be >= 0.\n");
    deallocate();
    initMembers();
    width_ = width;
    height_ = height;
    initLineStartArray();
    init(d);
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::resizeCopy(const CachedFileImage & rhs) {
    deallocate();
    initMembers();
    if (rhs.width() * rhs.height() > 0) {
        width_ = rhs.width();
        height_ = rhs.height();
        initLineStartArray();

        const_iterator is = rhs.begin();
        const_iterator iend = rhs.end();
        iterator id = begin();
        for(; is != iend; ++is, ++id) *id = *is;
    }
};

template <class PIXELTYPE>
void swap( CachedFileImage<PIXELTYPE>& rhs ) {
    if (&rhs != this) {
        std::swap(linesPerBlocksize_, rhs.linesPerBlocksize_);
        std::swap(blocksAllowed_, rhs.blocksAllowed_);
        std::swap(blockLRU_, rhs.blockLRU_);
        std::swap(lines_, rhs.lines_);
        std::swap(width_, rhs.width_);
        std::swap(height_, rhs.height_);
        std::swap(tmpFile_, rhs.tmpFile_);
        std::swap(tmpFilename_, rhs.tmpFilename_);
    }
};

/********************************************************/
/*                                                      */
/*              argument object factories               */
/*                                                      */
/********************************************************/

template <class PixelType, class Accessor>
inline triple<typename CachedFileImage<PixelType>::const_traverser, 
              typename CachedFileImage<PixelType>::const_traverser, Accessor>
srcImageRange(CachedFileImage<PixelType> const & img, Accessor a)
{
    return triple<typename CachedFileImage<PixelType>::const_traverser, 
                  typename CachedFileImage<PixelType>::const_traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline pair<typename CachedFileImage<PixelType>::const_traverser, Accessor>
srcImage(CachedFileImage<PixelType> const & img, Accessor a)
{
    return pair<typename CachedFileImage<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline triple<typename CachedFileImage<PixelType>::traverser, 
              typename CachedFileImage<PixelType>::traverser, Accessor>
destImageRange(CachedFileImage<PixelType> & img, Accessor a)
{
    return triple<typename CachedFileImage<PixelType>::traverser, 
                  typename CachedFileImage<PixelType>::traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline pair<typename CachedFileImage<PixelType>::traverser, Accessor>
destImage(CachedFileImage<PixelType> & img, Accessor a)
{
    return pair<typename CachedFileImage<PixelType>::traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename CachedFileImage<PixelType>::const_traverser, Accessor>
maskImage(CachedFileImage<PixelType> const & img, Accessor a)
{
    return pair<typename CachedFileImage<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

/****************************************************************/

template <class PixelType>
inline triple<typename CachedFileImage<PixelType>::const_traverser, 
              typename CachedFileImage<PixelType>::const_traverser, 
              typename CachedFileImage<PixelType>::ConstAccessor>
srcImageRange(CachedFileImage<PixelType> const & img)
{
    return triple<typename CachedFileImage<PixelType>::const_traverser, 
                  typename CachedFileImage<PixelType>::const_traverser, 
                  typename CachedFileImage<PixelType>::ConstAccessor>(img.upperLeft(),
                                                                 img.lowerRight(),
                                                                 img.accessor());
}

template <class PixelType>
inline pair< typename CachedFileImage<PixelType>::const_traverser, 
             typename CachedFileImage<PixelType>::ConstAccessor>
srcImage(CachedFileImage<PixelType> const & img)
{
    return pair<typename CachedFileImage<PixelType>::const_traverser, 
                typename CachedFileImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                                               img.accessor());
}

template <class PixelType>
inline triple< typename CachedFileImage<PixelType>::traverser, 
               typename CachedFileImage<PixelType>::traverser, 
               typename CachedFileImage<PixelType>::Accessor>
destImageRange(CachedFileImage<PixelType> & img)
{
    return triple<typename CachedFileImage<PixelType>::traverser, 
                  typename CachedFileImage<PixelType>::traverser, 
                  typename CachedFileImage<PixelType>::Accessor>(img.upperLeft(),
                                                            img.lowerRight(),
                                                            img.accessor());
}

template <class PixelType>
inline pair< typename CachedFileImage<PixelType>::traverser, 
             typename CachedFileImage<PixelType>::Accessor>
destImage(CachedFileImage<PixelType> & img)
{
    return pair<typename CachedFileImage<PixelType>::traverser, 
                typename CachedFileImage<PixelType>::Accessor>(img.upperLeft(), 
                                                          img.accessor());
}

template <class PixelType>
inline pair< typename CachedFileImage<PixelType>::const_traverser, 
             typename CachedFileImage<PixelType>::ConstAccessor>
maskImage(CachedFileImage<PixelType> const & img)
{
    return pair<typename CachedFileImage<PixelType>::const_traverser, 
                typename CachedFileImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                                               img.accessor());
}

} // namespace vigra

#endif /* VIGRA_EXT_CACHEDFILEIMAGE_HXX */
