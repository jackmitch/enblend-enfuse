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

#include <assert.h>
#include <errno.h>
#include <map>
#include <iostream>
#include <list>
#include <stdio.h>
#include <unistd.h>
#include <boost/pool/pool.hpp>
#include <vigra/utilities.hxx>

using std::cout;
using std::endl;
using std::list;
using std::map;

namespace vigra {

class CachedFileImageBase {
public:
    virtual int numBlocksAllocated() const = 0;
    virtual int numBlocksNeeded() const = 0;
    virtual void swapOutBlock() const = 0;
};

template <class PIXELTYPE> class CachedFileImage;

/** A singleton that manages memory for several CachedFileImages.
 */
class CachedFileImageDirector {
public:

    ~CachedFileImageDirector() {
        if (!imageList.empty()) {
            cout << "Cleaning up temporary files." << endl;
            list<CachedFileImageBase const *>::iterator i;
            for (i = imageList.begin(); i != imageList.end(); i++) {
                delete *i;
            }
        }
        delete pool;
    }

    // Obtain a reference to the singleton.
    static CachedFileImageDirector &v() {
        static CachedFileImageDirector instance;
        return instance;
    }

    // Obtain a pointer to a blocksize chunk of memory.
    void* allocateBlock() {
        void *block = pool->malloc();
        pool->set_next_size(1);
        if (block == NULL) throw std::bad_alloc();
        return block;
    }

    void deallocateBlock(void* block) {
        if (block != NULL) pool->free(block);
    }

    // Set the number of bytes this director manages.
    void setAllocation(long long bytes) {
        // This may not be changed after images have been created.
        assert(imageList.empty());
        managedBytes = bytes;
        // Recalculate the number of blocks available.
        managedBlocks = (int)ceil(managedBytes / (double)blocksize);
        blocksAvailable = managedBlocks;
    }

    // Set the cache block size. This is the minimum amount that is
    // moved between the caches and the backing files.
    void setBlockSize(int bytes) {
        // This may not be changed after images have been created.
        assert(imageList.empty());
        blocksize = bytes;
        delete pool;
        pool = new boost::pool<>(blocksize);
        // Recalculate the number of blocks available.
        managedBlocks = (int)ceil(managedBytes / (double)blocksize);
        blocksAvailable = managedBlocks;
    }

    // Get the cache block size.
    int getBlockSize() {
        return blocksize;
    }

    int getManagedBlocks() {
        return managedBlocks;
    }

    int getBlocksAvailable() {
        return blocksAvailable;
    }

    // Request a certain number of blocks for a new image.
    // Returns the number of blocks the new image may use.
    int requestBlocksForNewImage(int blocks, CachedFileImageBase const * image) {
        int blocksAllocated = 0;
        if (blocksAvailable >= blocks) {
            // Plenty of blocks available.
            // Give the image all the blocks it wants.
            blocksAllocated = blocks;
            blocksAvailable -= blocks;
        } else if (blocksAvailable > 0) {
            // Not enough blocks available.
            // Give the image as many blocks as are available.
            blocksAllocated = blocksAvailable;
            blocksAvailable = 0;
        } else {
            // Zero blocks available.
            // Try to free a block by forcing another image to swap.
            blocksAvailable += freeBlock();
            if (blocksAvailable == 0) {
                // Attempt to free a block was a failure.
                vigra_fail("CachedFileImageDirector::requestBlocksForNewImage(): "
                        "no blocks available and attempt to free blocks failed.");
            }
            blocksAllocated = blocksAvailable;
            blocksAvailable = 0;
        }

        // Register the new image.
        // By placing it at the back, it is marked as the
        // most-recently-swapped image.
        imageList.push_back(image);

        imageToMissMap[image] = 0LL;

        return blocksAllocated;
    }

    // Unregister a CachedFileImage with the director and return its blocks
    // to the pool.
    void returnBlocksUnregisterImage(int blocks, CachedFileImageBase const * image) {
        //cout << "returning " << blocks << " blocks" << endl;
        blocksAvailable += blocks;
        imageList.remove(image);
    }

    // Tell the director that a cache miss has occured for this image.
    // The director may decide to allocate more blocks to the image.
    // If so it returns a nonzero number.
    int registerCacheMiss(CachedFileImageBase const * image) {
        cacheMisses++;
        imageToMissMap[image]++;

        // Remove image from list.
        imageList.remove(image);

        if (blocksAvailable == 0) {
            // Try to free a block from an image that has not
            // missed recently. If this fails then we just won't
            // give the calling image more blocks.
            blocksAvailable += freeBlock();
        }

        // Add image to back of list.
        // This marks it most-recently-swapped
        imageList.push_back(image);

        if (blocksAvailable > 0) {
            // There are more blocks available to give out.
            // Give the image one more block.
            blocksAvailable--;
            return 1;
        } else {
            return 0;
        }
    }

    // How many cache misses have occured in total.
    long long getCacheMisses() {
        return cacheMisses;
    }

    long long getCacheMisses(CachedFileImageBase const * image) {
        return imageToMissMap[image];
    }

    void resetCacheMisses() {
        cacheMisses = 0LL;
        map<CachedFileImageBase const *, long long>::iterator i;
        for (i = imageToMissMap.begin(); i != imageToMissMap.end(); i++) {
            (*i).second = 0LL;
        }
    }

    void printStats() {
        cout << "Summary: cache misses="
             << cacheMisses
             << " blocks managed="
             << managedBlocks
             << " allocated="
             << (managedBlocks - blocksAvailable)
             << " free="
             << blocksAvailable
             << endl;
    }

    void printStats(const char * imageName, const CachedFileImageBase * image) {
        cout << imageName << ":"
             << " cache misses=" << imageToMissMap[image]
             << " blocks allocated=" << image->numBlocksAllocated()
             << " / " << image->numBlocksNeeded()
             << " required"
             << endl;
    }

    void printStats(const char * imageName, const int imageNumber,
            const CachedFileImageBase * image) {
        cout << imageName << imageNumber << ":"
             << " cache misses=" << imageToMissMap[image]
             << "  blocks allocated=" << image->numBlocksAllocated()
             << " / " << image->numBlocksNeeded()
             << " required"
             << endl;
    }

protected:
    CachedFileImageDirector()
    : blocksize(2<<20),
      managedBytes(1LL<<30),
      cacheMisses(0),
      imageList(),
      imageToMissMap()
    {
        // Recalculate the number of blocks available.
        managedBlocks = (int)ceil(managedBytes / (double)blocksize);
        blocksAvailable = managedBlocks;
        pool = new boost::pool<>(blocksize);
    }

    int freeBlock() {
        // Try to free a block from an image.
        // Check least-recently-missed images first.
        list<CachedFileImageBase const *>::iterator i;
        for (i = imageList.begin(); i != imageList.end(); i++) {
            CachedFileImageBase const * image = *i;
            if (image->numBlocksAllocated() > 0) {
                // Image is currently using blocks.
                image->swapOutBlock();
                // Mark image as most-recently-swapped.
                //imageList.erase(i);
                //imageList.push_back(image);
                return 1;
            }
        }

        // No blocks could be freed from other images.
        return 0;
    }

    int blocksize;
    long long managedBytes;
    int managedBlocks;
    int blocksAvailable;
    long long cacheMisses;

    // Pool for blocks
    boost::pool<> *pool;

    // List of images.
    // Front is least-recently-missed image
    // Back is most-recently-missed image
    list<CachedFileImageBase const *> imageList;

    // Map of images to cache misses
    map<CachedFileImageBase const *, long long> imageToMissMap;

};

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

    static void initialize(BaseType & d) { }

    static reference dereference(BaseType const & d) {
        return *d;
    }

    static index_reference dereference(BaseType const & d, difference_type n) {
        int width = d.i->width();
        int dy = n / width;
        int dx = n % width;
        if (d.x + dx >= width) {dy++; dx -= width;}
        else if (d.x + dx < 0) {dy--; dx += width;}
        return d(dx, dy);
    }

    static bool equal(BaseType const & d1, BaseType const & d2) {
        int width1 = d1.i->width();
        int width2 = d2.i->width();
        return (d1.y*width1 + d1.x) == (d2.y*width2 + d2.x);
    }

    static bool less(BaseType const & d1, BaseType const & d2) {
        int width1 = d1.i->width();
        int width2 = d2.i->width();
        return (d1.y*width1 + d1.x) < (d2.y*width2 + d2.x);
    }

    static difference_type difference(BaseType const & d1, BaseType const & d2) {
        int width1 = d1.i->width();
        int width2 = d2.i->width();
        return (d1.y*width1 + d1.x) - (d2.y*width2 + d2.x);
    }

    static void increment(BaseType & d) {
        ++d.x;
        if (d.x == d.i->width()) {
            d.x = 0;
            ++d.y;
        }
    }

    static void decrement(BaseType & d) {
        --d.x;
        if (d.x < 0) {
            d.x = d.i->width() - 1;
            --d.y;
        }
    }

    static void advance(BaseType & d, difference_type n) {
        int width = d.i->width();
        int dy = n / width;
        int dx = n % width;
        d.x += dx;
        d.y += dy;
        if (d.x >= width) {++d.y; d.x -= width;}
        if (d.x < 0) {--dy.y; d.x += width;}
    }

};

template <class IMAGEITERATOR, class IMAGETYPE, class PIXELTYPE, class REFERENCE, class POINTER>
class CachedFileImageIteratorBase
{
public:
    typedef CachedFileImageIteratorBase<IMAGEITERATOR,
            IMAGETYPE, PIXELTYPE, REFERENCE, POINTER> self_type;
    typedef IMAGETYPE image_type;
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
        return row_iterator(static_cast<IMAGEITERATOR const &>(*this));
    }

    column_iterator columnIterator() const {
        return column_iterator(static_cast<IMAGEITERATOR const &>(*this));
    }

protected:

    CachedFileImageIteratorBase(const int X, const int Y, image_type * const I) : x(X), y(Y), i(I) { }

    CachedFileImageIteratorBase() : x(0), y(0), i(NULL) { }

};
    
template <class PIXELTYPE>
class CachedFileImageIterator
: public CachedFileImageIteratorBase<CachedFileImageIterator<PIXELTYPE>,
                CachedFileImage<PIXELTYPE>,
                PIXELTYPE, PIXELTYPE &, PIXELTYPE *>
// FIXME this needs to be a weak_ptr    ^^^^^^^^^^^
// in case someone uses the iterator to get a pointer to cached data.
{
public:

    typedef CachedFileImageIteratorBase<CachedFileImageIterator,
            CachedFileImage<PIXELTYPE>,
            PIXELTYPE, PIXELTYPE &, PIXELTYPE *> Base;

    CachedFileImageIterator(const int x, const int y, CachedFileImage<PIXELTYPE> * const i)
    : Base(x, y, i)
    {}

    CachedFileImageIterator()
    : Base(0, 0, NULL)
    {}

};

template <class PIXELTYPE>
class ConstCachedFileImageIterator
: public CachedFileImageIteratorBase<ConstCachedFileImageIterator<PIXELTYPE>,
                const CachedFileImage<PIXELTYPE>,
                PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *>
// FIXME this needs to be a weak_ptr          ^^^^^^^^^^^^^^^^^
// in case someone uses the iterator to get a pointer to cached data.
{
public:

    typedef CachedFileImageIteratorBase<ConstCachedFileImageIterator,
            const CachedFileImage<PIXELTYPE>,
            PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *> Base;
    // FIXME this needs to be a weak_ptr  ^^^^^^^^^^^^^^^^^

    ConstCachedFileImageIterator(const int x, const int y, const CachedFileImage<PIXELTYPE> * const i)
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

//class CachedFileImageBase {
//public:
//    virtual int numBlocksAllocated() const = 0;
//    virtual void swapOutBlock() const = 0;
//};

template <class PIXELTYPE>
class CachedFileImage : public CachedFileImageBase {
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

    //struct Allocator {
    //    static value_type * allocate(int n) {
    //        return (value_type *)::operator new(n*sizeof(value_type));
    //    }
    //    static void deallocate(value_type * p) {
    //        ::operator delete(p);
    //    }
    //};

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

    virtual ~CachedFileImage() {
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
        //cout << "dirty" << endl;
        return (getLinePointerDirty(d.y))[d.x];
    }

    const_reference operator[](difference_type const & d) const {
        //cout << "clean" << endl;
        return (getLinePointer(d.y))[d.x];
    }

    reference operator()(int dx, int dy) {
        //cout << "dirty" << endl;
        return (getLinePointerDirty(dy))[dx];
    }

    const_reference operator()(int dx, int dy) const {
        //cout << "clean" << endl;
        return (getLinePointer(dy))[dx];
    }

    // dangerous - needs to return a weak_ptr
    pointer operator[](int dy) {
        //cout << "dirty" << endl;
        return getLinePointerDirty(dy);
    }

    // dangerous - needs to return a weak_ptr
    const_pointer operator[](int dy) const {
        //cout << "clean" << endl;
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

    int numBlocksAllocated() const {
        return blocksAllocated_;
    }

    int numBlocksNeeded() const {
        return blocksNeeded_;
    }

    void swapOutBlock() const {
        swapLeastRecentlyUsedBlock();
    }

private:

    PIXELTYPE initPixel;

    void deallocate();

    void initLineStartArray();
    
    // obtain a pointer to the beginning of a line.
    // split into two functions for efficiency.
    // getLinePointer can then be inlined, and we only incur the function call overhead
    // on cache misses.
    PIXELTYPE * getLinePointer(int dy) const;
    PIXELTYPE * getLinePointerDirty(int dy);
    PIXELTYPE * getLinePointerCacheMiss(int dy) const;

    // Free space, if necessary, by swapping out a block of lines to the file.
    void swapLeastRecentlyUsedBlock() const;

    // Lazy creation of tmp file for swapping image data to disk
    void initTmpfile() const;

    //inline int lineToBlockNumber(int line) const {
    //    return line / linesPerBlocksize_;
    //}

    //inline int blockToFirstLineNumber(int block) const {
    //    return block * linesPerBlocksize_;
    //}

    void initMembers() {
        initPixel = value_type();
        linesPerBlocksize_ = 0;
        blocksAllocated_ = 0;
        blocksNeeded_ = 0;
        blockLRU_ = NULL;
        lines_ = NULL;
        blockIsClean_ = NULL;
        blockInFile_ = NULL;
        width_ = 0;
        height_ = 0;
        tmpFile_ = NULL;
        tmpFilename_ = NULL;
    }

    // how many image lines are loaded in one fread
    int linesPerBlocksize_;
    // How many blocks are currently cached in memory.
    mutable int blocksAllocated_;
    int blocksNeeded_;

    // lru replacement policy
    // most recently used block is at start of list.
    // least recently used block is at end of list.
    mutable list<int> *blockLRU_;

    mutable PIXELTYPE ** lines_;
    mutable bool * blockIsClean_;
    mutable bool * blockInFile_;

    int width_, height_;

    mutable FILE *tmpFile_;
    mutable char *tmpFilename_;
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::deallocate() {
    CachedFileImageDirector::v().returnBlocksUnregisterImage(blocksAllocated_, this);
    delete blockLRU_;
    if (lines_ != NULL) {
        // Go through lines and delete any allocated memory there.
        //for (int line = 0; line < height_; line++) {
        //    PIXELTYPE *p = lines_[line];
        //    if (p != NULL) {
        //        for (int column = 0; column < width_; column++) {
        //            (p[column]).~PIXELTYPE();
        //        }
        //        Allocator::deallocate(p);
        //    }
        //}
        int line = 0;
        for (int block = 0; block < blocksNeeded_; block++) {
            int firstLineInBlock = line;
            for (int subblock = 0; subblock < linesPerBlocksize_; subblock++, line++) {
                if (line >= height_) break;
                PIXELTYPE *p = lines_[line];
                if (p != NULL) {
                    for (int column = 0; column < width_; column++) {
                        (p[column]).~PIXELTYPE();
                    }
                }
            }
            CachedFileImageDirector::v().deallocateBlock(lines_[firstLineInBlock]);
            if (line >= height_) break;
        }
        delete[] lines_;
    }
    delete[] blockIsClean_;
    delete[] blockInFile_;
    if (tmpFile_ != NULL) {
        fclose(tmpFile_);
        unlink(tmpFilename_);
    }
    delete[] tmpFilename_;
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::initLineStartArray() {

    // Number of lines to load in one block.
    linesPerBlocksize_ = (int)floor(
            ((double)CachedFileImageDirector::v().getBlockSize())
            / (width_ * sizeof(PIXELTYPE)));

    if (linesPerBlocksize_ <= 0) {
        throw vigra::InvariantViolation("Image cache block size is too small.");
    }

    blocksNeeded_ = (int)ceil(((double)height_) / linesPerBlocksize_);

    int blocksAllowed = CachedFileImageDirector::v().requestBlocksForNewImage(
            blocksNeeded_, this);

    // Create the blockLRU list.
    blockLRU_ = new list<int>();

    lines_ = new PIXELTYPE*[height_];
    blockIsClean_ = new bool[blocksNeeded_];
    blockInFile_ = new bool[blocksNeeded_];

    // Initialize blockIsClean / blockInFile vectors.
    for (int block = 0; block < blocksNeeded_; block++) {
        blockIsClean_[block] = true;
        blockInFile_[block] = false;
    }

    // Allocate mem for the first linesPerBlocksize_*blocksAllowed lines.
    int line = 0;
    for (int block = 0; block < blocksAllowed; block++) {
        blockLRU_->push_front(block);
        blocksAllocated_++;

        // Get a block from the director.
        PIXELTYPE* blockStart = (PIXELTYPE*)CachedFileImageDirector::v().allocateBlock();

        // Divide the block up amongst the lines in the block.
        for (int subblock = 0; subblock < linesPerBlocksize_; subblock++, line++, blockStart+=width_) {
            if (line >= height_) break;
            //lines_[line] = Allocator::allocate(width_);
            lines_[line] = blockStart;
            std::uninitialized_fill_n(lines_[line], width_, initPixel);
        }
        if (line >= height_) break;
    }

    // All remaining lines (if any) are null (swapped out)
    for (; line < height_; line++) {
        lines_[line] = NULL;
    }

    return;
};

template <class PIXELTYPE>
inline PIXELTYPE * CachedFileImage<PIXELTYPE>::getLinePointerDirty(int dy) {
    PIXELTYPE *line = lines_[dy];

    // Check if line dy is swapped out.
    if (line == NULL) line = getLinePointerCacheMiss(dy);

    // Mark this block as dirty.
    blockIsClean_[dy / linesPerBlocksize_] = false;

    return line;
};

template <class PIXELTYPE>
inline PIXELTYPE * CachedFileImage<PIXELTYPE>::getLinePointer(int dy) const {
    PIXELTYPE *line = lines_[dy];

    // Check if line dy is swapped out.
    if (line == NULL) {
        return getLinePointerCacheMiss(dy);
    }
    else {
        // make blockNumber the least recently used block.
        // this remove function call really sucks
        //blockLRU_->remove(blockNumber);
        //blockLRU_->push_front(blockNumber);
        return line;
    }
};

template <class PIXELTYPE>
PIXELTYPE * CachedFileImage<PIXELTYPE>::getLinePointerCacheMiss(int dy) const {
    int blockNumber = dy / linesPerBlocksize_; // lineToBlockNumber(dy);
    int firstLineInBlock = blockNumber * linesPerBlocksize_; // blockToFirstLineNumber(blockNumber);

    int moreBlocks = CachedFileImageDirector::v().registerCacheMiss(this);
    if (moreBlocks == 0 && blocksAllocated_ == 0) {
        vigra_fail("CachedFileImage::getLinePointerCacheMiss(): no blocks available "
                " and attempt to free blocks failed.");
    } else if (moreBlocks == 0) {
        // Make space for new block.
        swapLeastRecentlyUsedBlock();
    }
    //cout << "swapping in block " << blockNumber << endl;

    // Allocate a block.
    PIXELTYPE* blockStart = (PIXELTYPE*)CachedFileImageDirector::v().allocateBlock();

    int numLinesInBlock = min(height_, firstLineInBlock + linesPerBlocksize_) - firstLineInBlock;
    int pixelsToRead = numLinesInBlock * width_;

    if (blockInFile_[blockNumber]) {
        // Find the right spot in the file.
        off_t offset = width_ * firstLineInBlock * sizeof(PIXELTYPE);
        if (fseeko(tmpFile_, offset, SEEK_SET) != 0) {
            vigra_fail(strerror(errno));
        }
        // Fill the block with data from the file.
        int itemsRead = fread(blockStart, sizeof(PIXELTYPE), pixelsToRead, tmpFile_);
        if (itemsRead < pixelsToRead) {
            vigra_fail("CachedFileImage: error reading from image backing file.\n");
        }
    }
    else {
        // File does not have data for this block.
        // Fill lines with initPixel.
        std::uninitialized_fill_n(blockStart, pixelsToRead, initPixel);
    }

    // Divide the block up amongst the lines in the block.
    for (int l = 0; l < linesPerBlocksize_; l++, blockStart+=width_) {
        int absoluteLineNumber = l + firstLineInBlock;
        if (absoluteLineNumber >= height_) break;
        lines_[absoluteLineNumber] = blockStart;
    }
    //    //// Allocate lines for new block.
    //    //for (int l = 0; l < linesPerBlocksize_; l++) {
    //    //    int absoluteLineNumber = l + firstLineInBlock;
    //    //    if (absoluteLineNumber >= height_) break;
    //    //    lines_[absoluteLineNumber] = Allocator::allocate(width_);

    //    //    // fill the line with data from the file.
    //    //    int itemsRead = fread(lines_[absoluteLineNumber], sizeof(PIXELTYPE), width_, tmpFile_);

    //    //    if (itemsRead < width_) {
    //    //        vigra_fail("CachedFileImage: error reading from image backing file.\n");
    //    //    }
    //    //}
    //}
    //else {
    //    // File does not have data for this block. Create new lines and fill them with initPixel.
    //    // Allocate lines for new block.

    //    // Allocate a block.
    //    PIXELTYPE* blockStart = (PIXELTYPE*)CachedFileImageDirector::v().allocateBlock();

    //    // Divide the block up amongst the lines in the block.
    //    for (int l = 0; l < linesPerBlocksize_; l++, blockStart+=width_) {
    //        int absoluteLineNumber = l + firstLineInBlock;
    //        if (absoluteLineNumber >= height_) break;
    //        lines_[absoluteLineNumber] = blockStart;
    //        std::uninitialized_fill_n(lines_[absoluteLineNumber], width_, initPixel);
    //    }

    //    //for (int l = 0; l < linesPerBlocksize_; l++) {
    //    //    int absoluteLineNumber = l + firstLineInBlock;
    //    //    if (absoluteLineNumber >= height_) break;
    //    //    lines_[absoluteLineNumber] = Allocator::allocate(width_);

    //    //    // fill with initPixel.
    //    //    std::uninitialized_fill_n(lines_[absoluteLineNumber], width_, initPixel);
    //    //}
    //}

    // Mark this block as clean.
    blockIsClean_[blockNumber] = true;

    // Mark this block as most recently used.
    blockLRU_->push_front(blockNumber);
    blocksAllocated_++;

    return lines_[dy];
};


template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::swapLeastRecentlyUsedBlock() const {
    if (blocksAllocated_ == 0) return;

    int blockNumber = blockLRU_->back();

    //list<int>::iterator listIterator = blockLRU_->begin();
    //for(; listIterator != blockLRU_->end(); listIterator++) {
    //    cout << *listIterator << " ";
    //}
    //cout << endl << "swapping out block " << blockNumber << endl;
    //cout << "linesPerBlocksize=" << linesPerBlocksize_ << endl;

    blockLRU_->pop_back();

    int firstLineInBlock = blockNumber * linesPerBlocksize_;
    PIXELTYPE *blockStart = lines_[firstLineInBlock];

    // If block is dirty, swap it to the file.
    if (!blockIsClean_[blockNumber]) {
        blockInFile_[blockNumber] = true;

        // Lazy init the temp file.
        if (tmpFile_ == NULL) initTmpfile();

        // Find the right spot in the file.
        off_t offset = width_ * firstLineInBlock * sizeof(PIXELTYPE);
        if (fseeko(tmpFile_, offset, SEEK_SET) != 0) {
            vigra_fail(strerror(errno));
        }

        int numLinesInBlock = min(height_, firstLineInBlock + linesPerBlocksize_) - firstLineInBlock;
        int pixelsToWrite = numLinesInBlock * width_;
        int itemsWritten = fwrite(blockStart, sizeof(PIXELTYPE), pixelsToWrite, tmpFile_);
        if (itemsWritten < pixelsToWrite) {
            vigra_fail("CachedFileImage: error writing to image backing file.\n");
        }
    }

    for (int l = 0; l < linesPerBlocksize_; l++) {
        int absoluteLineNumber = l + firstLineInBlock;
        if (absoluteLineNumber >= height_) break;
        //cout << "swapping line " << absoluteLineNumber << endl;
        //if (fwrite(lines_[absoluteLineNumber], sizeof(PIXELTYPE), width_, tmpFile_)
        //        != (unsigned int)width_) {
        //    vigra_fail("CachedFileImage: error writing to image backing file.\n");
        //}
        // Deallocate line
        PIXELTYPE *p = lines_[absoluteLineNumber];
        for (int column = 0; column <= width_; column++) {
            //FIXME if pixel type is not a simple data type and this destructor actually does
            // something, then we are in big trouble.
            (p[column]).~PIXELTYPE();
        }
        //Allocator::deallocate(p);
        lines_[absoluteLineNumber] = NULL;
    }

    CachedFileImageDirector::v().deallocateBlock(blockStart);

    //}
    //else {
    //    // Block is clean - just deallocate it.
    //    PIXELTYPE *blockStart = lines
    //    for (int l = 0; l < linesPerBlocksize_; l++) {
    //        int absoluteLineNumber = l + firstLineInBlock;
    //        if (absoluteLineNumber >= height_) break;
    //        // Deallocate line
    //        PIXELTYPE *p = lines_[absoluteLineNumber];
    //        for (int column = 0; column <= width_; column++) {
    //            //FIXME if pixel type is not a simple data type and this destructor actually does
    //            // something, then we are in big trouble.
    //            (p[column]).~PIXELTYPE();
    //        }
    //        Allocator::deallocate(p);
    //        lines_[absoluteLineNumber] = NULL;
    //    }
    //}

    blocksAllocated_--;
};

template <class PIXELTYPE>
void CachedFileImage<PIXELTYPE>::initTmpfile() const {
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
    //if (setvbuf(tmpFile_, NULL, _IONBF, 0) != 0) {
    //    vigra_fail(strerror(errno));
    //}
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
    initPixel = pixel;

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
    initPixel = d;
    width_ = width;
    height_ = height;
    initLineStartArray();
    //this should do a fast init.
    //init(d);
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
void CachedFileImage<PIXELTYPE>::swap( CachedFileImage<PIXELTYPE>& rhs ) {
    if (&rhs != this) {
        std::swap(initPixel, rhs.initPixel);
        std::swap(linesPerBlocksize_, rhs.linesPerBlocksize_);
        std::swap(blocksAllocated_, rhs.blocksAllocated_);
        std::swap(blocksNeeded_, rhs.blocksNeeded_);
        std::swap(blockLRU_, rhs.blockLRU_);
        std::swap(lines_, rhs.lines_);
        std::swap(blockIsClean_, rhs.blockIsClean_);
        std::swap(blockInFile_, rhs.blockInFile_);
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
