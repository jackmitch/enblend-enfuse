/*
 * Copyright (C) 2013-2014 Christoph L. Spiel
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
#ifndef MMAP_IMAGE_H_INCLUDED
#define MMAP_IMAGE_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <array>
#include <map>
#include <set>
#include <stdexcept>            // std::exception
#include <string>
#include <vector>

#include <vigra/basicimageview.hxx>
#include <vigra/diff2d.hxx>

#include "openmp_lock.h"


namespace mmap_image
{
    namespace sys
    {
        class error : public std::exception
        {
        public:
            error() = delete;
            explicit error(int an_error_number) throw();
            error(const error& an_error) throw();
            error& operator=(const error& an_error) throw();
            ~error() throw();

            const char* what() const throw() {return message_;}
            int number() const throw() {return error_number_;}

        private:
            int error_number_;
            char* message_;
        }; // class error


        std::string get_tmpdir();
    } // namespace sys


    // Map memory backed by a file.
    class MemoryMapFile
    {
    public:
        MemoryMapFile() = delete;
        MemoryMapFile(const MemoryMapFile&) = delete;
        MemoryMapFile& operator=(const MemoryMapFile&) = delete;

        MemoryMapFile(size_t a_size);
        MemoryMapFile(size_t a_size, const std::string& a_filename);
        ~MemoryMapFile();

        bool resize(size_t a_size, bool may_move = false);

        size_t size() const {return size_;}
        void* data() const {return map_;}

        static std::vector<std::string> open_map_files();

    protected:
        void initialize();
        void finalize();

    private:
        size_t size_;
        std::string filename_;
        const int desc_;
        void* map_;

        static omp::lock open_files_lock_;
        static std::set<std::string> open_files_;
    }; // end class MemoryMapFile


    class BuddyHistogram
    {
    public:
        explicit BuddyHistogram(unsigned a_number_of_bins = 48U);

        unsigned number_of_bins() const;
        unsigned total_count() const;

        void insert(size_t a_value);
        unsigned count(unsigned a_bin_index) const;
        unsigned cumulative_count(unsigned a_bin_index) const;
        bool has_overflown() const;

    private:
        unsigned total_count_;
        std::vector<unsigned> histogram_;
        bool overflow_;
    }; // end class BuddyHistogram


    // DESIGN DECISIONS
    //
    // 1. We want to derive from vigra::BasicImageView to get the
    //    almost the power of vigra::BasicImage without having to
    //    override lots of methods.
    //
    // 2. vigra::BasicImageView needs the start address of the image
    //    in the call to the constructor.  The address cannot be
    //    changed later; vigra::BasicImageView<>::data() is a reader.
    //
    // 3. To enforce correct evaluation order we first derive
    //    (private) from MemoryDirector which gives us the address.
    //
    // Side Note: We do not loose too much typing precision by
    //            returning the addresses as void*.  That way, we can
    //            keep MemoryDirector non-templated.

    class MemoryDirector
    {
    public:
        MemoryDirector() = delete;
        MemoryDirector(const MemoryDirector&) = delete;
        MemoryDirector& operator=(const MemoryDirector&) = delete;

        explicit MemoryDirector(size_t a_size);

        virtual ~MemoryDirector();

        virtual bool hold_in_core(size_t a_size);

        static bool force_in_core(bool do_enforce); // short-ciruit any algorithm

        bool is_in_core() const {return is_in_core_;}
        void* base_address() const {return is_in_core_ ? image_data_ : map_file_->data();}

    private:
        enum hcr                      // Hold-in-Core Rule ids ("HCR")
        {
            HCR_ON_DISK,              // Unconditionally prefer out-of-core storage
            HCR_IN_CORE,              // Unconditionally use in-core memory
            HCR_SINGLE_THRESHOLD      // Decide on a single size threshold
        };

        static bool do_enforce_in_core_;
        enum hcr rule_id_;

        const bool is_in_core_;

        MemoryMapFile* map_file_;
        void* image_data_;
    }; // end class MemoryDirector


    class ScopedInCore
    {
    public:
        ScopedInCore() = delete;
        ScopedInCore(MemoryDirector* a_memory_manager);
        virtual ~ScopedInCore();

    private:
        MemoryDirector* memory_manager_;
        bool previous_force_state_;
    }; // end class ScopedInCore


    struct MemoryDirectorStatistics
    {
        size_t number;                // number of allocations
        size_t number_out_of_core;    // number of allocations outside of core

        size_t size;                  // size of the allocations in bytes
        size_t size_out_of_core;      // size of the out-of-core allocations
        BuddyHistogram histogram;     // histogram of sizes

        void show_allocation(const std::string& a_label) const;
        void show_histogram(const std::string& a_label) const;
    }; // end struct MemoryDirectorStatistics


    class MemoryDirectorWithStatistics : public MemoryDirector
    {
        typedef MemoryDirector super;
        typedef std::map<void*, size_t> address_map_t;

    public:
        MemoryDirectorWithStatistics() = delete;
        MemoryDirectorWithStatistics(const MemoryDirectorWithStatistics&) = delete;
        MemoryDirectorWithStatistics& operator=(const MemoryDirectorWithStatistics&) = delete;

        explicit MemoryDirectorWithStatistics(size_t a_size);

        ~MemoryDirectorWithStatistics();

        bool hold_in_core(size_t a_size) override;

        enum // three MemoryDirector Statistics ("MDS")
        {
            MDS_TOTAL,                  // grand totals
            MDS_CURRENT,                // usage right now
            MDS_MAX,                    // maxima aka high-water marks
            MDS_Number_Of_Statistics
        };

        static std::array<MemoryDirectorStatistics, MDS_Number_Of_Statistics> statistics;
        static void dump_allocations();
        static void dump_statistics();

    private:
        enum hcr                               // Hold-in-Core Rule ids
        {
            HCR_ON_DISK,                       // Unconditionally prefer out-of-core storage
            HCR_IN_CORE,                       // Unconditionally use in-core memory
            HCR_SINGLE_THRESHOLD,              // Decide on a single size threshold
            HCR_DOUBLE_THRESHOLD_AND_POOL_SIZE // Decide on two thresholds and the current
                                               // size of the in-core memory
        };

        enum hcr rule_id_;
        address_map_t sizes_;
    }; // end class MemoryDirectorWithStatistics


    template <class pixel_type>
    class MemoryMappedImage :
        private MemoryDirectorWithStatistics,
        public vigra::BasicImageView<pixel_type>
    {
    public:
        typedef vigra::BasicImageView<pixel_type> super;

        MemoryMappedImage() = delete;
        MemoryMappedImage(const MemoryMappedImage&) = delete;
        MemoryMappedImage& operator=(const MemoryMappedImage&) = delete;

        MemoryMappedImage(const vigra::Diff2D& a_size);
        MemoryMappedImage(int a_width, int a_height);
    }; // end class MemoryMappedImage
} // end namespace mmap_image


#endif // MMAP_IMAGE_H_INCLUDED

// Local Variables:
// mode: c++
// End:
