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


#include <cassert>
#include <cstring>              // strcpy, strlen
#include <iomanip>
#include <iostream>
#include <memory>               // std::unique_ptr
#include <numeric>              // std::accumulate

#include <fcntl.h>              // O_CREAT, O_RDWR, open(), ...
#include <stdlib.h>             // mkostemp()
#include <sys/stat.h>           // fstat()


#ifdef _WIN32
#include <mman.h>
#include <io.h>

#define _USE_MATH_DEFINES
#define NOMINMAX
#define VC_EXTRALEAN
#include <windows.h>            // SYSTEM_INFO, GetNativeSystemInfo
#undef DIFFERENCE

#if defined _MSC_VER
#define __attribute__(x)
#endif

#else

#include <sys/mman.h>           // mmap(), munmap(), ...
#include <unistd.h>             // _SC_PAGESIZE, lseek(), sysconf()
#endif


#include <vigra/rgbvalue.hxx>
#include <vigra/basicimage.hxx>

#include "parameter.h"

#include "mmap_image.h"


namespace mmap_image
{
#ifdef _WIN32
    enum
    {
        MAP_NORESERVE = 0x4000,
        MADV_RANDOM = 1
    };
#endif


    namespace sys
    {
        error::error(int an_error_number) throw() :
            error_number_(an_error_number), message_(nullptr)
        {
            char* m = strerror(an_error_number);
            message_ = new char[strlen(m) + 1U];
            strcpy(message_, m);
        }


        error::error(const error& an_error) throw() :
            error_number_(an_error.error_number_),
            message_(new char[strlen(an_error.message_) + 1U])
        {
            strcpy(message_, an_error.message_);
        }


        error&
        error::operator=(const error& an_error) throw()
        {
            if (this != &an_error)
            {
                error_number_ = an_error.error_number_;
                delete [] message_;
                message_ = new char[strlen(an_error.message_) + 1U];
                strcpy(message_, an_error.message_);
            }

            return *this;
        }


        error::~error() throw()
        {
            delete [] message_;
        }


#ifdef _WIN32
#ifndef off64_t
#define off64_t off_t
#endif

        inline static size_t
        pagesize()
        {
            SYSTEM_INFO info;
            GetNativeSystemInfo(&info);

            return info.dwPageSize;
        }


        // Forced write to disk
        static int
        fsync(int fd)
        {
            HANDLE h = (HANDLE)_get_osfhandle(fd);
            DWORD err;

            if (h == INVALID_HANDLE_VALUE)
            {
                errno = EBADF;
                return -1;
            }
            if (!FlushFileBuffers(h))
            {
                /* Windows error -> Unix */
                err = GetLastError();
                switch (err)
                {
                case ERROR_INVALID_HANDLE:
                    errno = EINVAL;
                    break;
                default:
                    errno = EIO;
                }
                return -1;
            }
            return 0;
        }


        inline static int
        fdatasync(int fd)
        {
            return fsync(fd);
        }


        // Set size with 64-bit support
        static int
        ftruncate(int fd, off64_t length)
        {
            HANDLE h = (HANDLE)_get_osfhandle(fd);
            LARGE_INTEGER l;
            LARGE_INTEGER o;

            if (h == INVALID_HANDLE_VALUE)
            {
                errno = EBADF;
                return -1;
            }
            l.QuadPart = length;
            if (!SetFilePointerEx(h, l, &o, FILE_BEGIN))
            {
                return -1;
            }
            if (!SetEndOfFile(h))
            {
                return -1;
            }

            return 0;
        }

#else

        inline static size_t
        pagesize()
        {
            return sysconf(_SC_PAGESIZE);
        }
#endif


        enum
        {
            OPEN_FLAGS = O_CREAT | O_TRUNC | O_RDWR,
            FILE_PERMISSIONS = 0600,
            PROTECTION = PROT_READ | PROT_WRITE,
            MAPPING_FLAGS = MAP_SHARED | MAP_NORESERVE
        };


        inline static size_t
        block_size()
        {
            return std::max(static_cast<size_t>(16 * 4096), pagesize());
        }


        static off_t
        get_filesize(int a_descriptor)
        {
            std::unique_ptr<struct stat> buffer(new struct stat);

            if (fstat(a_descriptor, buffer.get()) == -1)
            {
                throw error(errno);
            }

            return buffer->st_size;
        }


        static void
        initialize_file(int a_descriptor, size_t a_size)
        {
            if (lseek(a_descriptor, off_t(), SEEK_SET) == -1)
            {
                throw error(errno);
            }

            const size_t buffer_size = block_size();
            std::unique_ptr<char> buffer(new char[buffer_size]);
            memset(buffer.get(), int(), buffer_size);

            // IMPLEMENTATION NOTE: We must add 2 buffer_sizes,
            // one is at the beginning to give us an address at a
            // page boundary; the second absorbs the truncation of
            // the integral division.
            const size_t n = a_size / buffer_size + 2U;

#ifdef HAVE_ASYNC_IO
            struct aiocb aio;

            aio.aio_fildes = a_descriptor;
            aio.aio_buf = buffer.get();
            aio.aio_nbytes = buffer_size;
            aio.aio_sigevent = SIGEV_NONE; // do not signal any aio_write() call

            for (size_t i = 0U; i != n; ++i)
            {
                aio.aio_offset = i * buffer_size;
                if (aio_write(&aio) == -1)
                {
                    throw error(errno);
                }
            }

            std::array<struct aiocb*, 1U> aios{&aio};
            if (aio_suspend(aios.data(), aios.size(), nullptr) == -1)
            {
                throw error(errno);
            }
#else
            for (size_t i = 0U; i != n; ++i)
            {
                if (write(a_descriptor, buffer.get(), buffer_size) == -1)
                {
                    throw error(errno);
                }
            }

            if (fdatasync(a_descriptor) == -1)
            {
                throw error(errno);
            }
#endif
        }


        static void
        resize_file(int a_descriptor, off_t a_size)
        {
            const size_t current_size = get_filesize(a_descriptor) / block_size();
            // See for comment in initialize_file() why we have to add 2.
            const size_t new_size = a_size /  block_size() + 2U;

            if (new_size < current_size)
            {
                if (ftruncate(a_descriptor, new_size * block_size()) == -1)
                {
                    throw error(errno);
                }
                if (fsync(a_descriptor) == -1) // prefer fsync() as it flushes meta data, too
                {
                    throw error(errno);
                }
            }
            else if (new_size > current_size)
            {
                initialize_file(a_descriptor, a_size);
            }
            else
            {
                ;                   // file size is unchanged -- do nothing
            }
        }


        static int
        open_file(const std::string& a_filename, int a_flag_set, int a_permission_set, size_t a_size)
        {
            const int descriptor = open(a_filename.c_str(), a_flag_set, a_permission_set);

            if (descriptor == -1)
            {
                throw error(errno);
            }
            else
            {
                initialize_file(descriptor, a_size);

                return descriptor;
            }
        }


        static void*
        map_memory(size_t a_size, int a_protection_request, int a_flag_set, int a_descriptor)
        {
            void* map = mmap(nullptr, a_size, a_protection_request, a_flag_set, a_descriptor, off_t());

            if (reinterpret_cast<long int>(map) == -1L)
            {
                throw error(errno);
            }

            return map;
        }


        static void*
        remap_memory(void* an_address, size_t an_old_size, size_t a_new_size,
                     __attribute__((unused)) int a_protection_request,
                     __attribute__((unused)) int a_flag_set,
                     __attribute__((unused)) int a_descriptor,
                     bool may_move)
        {
#ifdef _GNU_SOURCE
            void* map = mremap(an_address, an_old_size, a_new_size, may_move ? MREMAP_MAYMOVE : 0);
#else
            if (munmap(an_address, an_old_size))
            {
                throw error(errno);
            }

            void* map = mmap(may_move ? nullptr : an_address, a_new_size,
                             a_protection_request, a_flag_set, a_descriptor, off_t());
#endif

            if (reinterpret_cast<long int>(map) == -1L)
            {
                throw error(errno);
            }

            return map;
        }


#ifdef _WIN32
        static std::string
        map_filename_template()
        {
            // On Windows we are using the function GetTempFileName
            // this function does only support a prefix with three
            // characters.
            return std::string("en");
        }


        /* extern */ std::string
        get_tmpdir()
        {
            char lpPathBuffer[MAX_PATH];
            const DWORD dwRetVal = GetTempPath(MAX_PATH, lpPathBuffer);

            if (dwRetVal >= MAX_PATH || dwRetVal == 0)
            {
                throw std::runtime_error("could not get path to temporary directory");
            }

            return std::string(lpPathBuffer);
        }


        inline static void
        advise_memory(void* an_address, size_t a_size, int an_advice)
        {
            // empty
        }


        static int
        open_temporary_file(/* output */ std::string& a_filename, const std::string& a_filename_template,
                            __attribute__((unused)) int a_flag_set, size_t a_size)
        {
            char temporary_filename[MAX_PATH];
            if (GetTempFileName(get_tmpdir().c_str(), a_filename_template.c_str(), 0, temporary_filename) == 0)
            {
                throw std::runtime_error("Could not create tempfile name.");
            };

            std::unique_ptr<char> writable_template(new char[strlen(temporary_filename) + 1U]);
            strcpy(writable_template.get(), temporary_filename);

            const int descriptor = open(temporary_filename, _O_RDWR | _O_BINARY | _O_TEMPORARY);

            if (descriptor == -1)
            {
                throw error(errno);
            }

            a_filename = writable_template.get();
            initialize_file(descriptor, a_size);

            return descriptor;
        }

#else

        static std::string
        map_filename_template()
        {
            enum {HOSTNAME_SIZE = 256U}; // HOSTNAME_SIZE >= HOST_NAME_MAX
            std::ostringstream oss;
            std::array<char, HOSTNAME_SIZE> hostname;

            oss << "en_mmap_image_";
            if (gethostname(hostname.data(), HOSTNAME_SIZE) == 0)
            {
                oss << hostname.data() << '_';
            }
            oss  << getpid() << "_XXXXXX";

            return oss.str();
        };


        /* extern */ std::string
        get_tmpdir()
        {
            const char* tmpdir = getenv("TMPDIR");

            if (tmpdir && *tmpdir != 0)
            {
                return std::string(tmpdir) + "/";
            }
            else
            {
                return std::string("/tmp/");
            }
        }


        inline static void
        advise_memory(void* an_address, size_t a_size, int an_advice)
        {
            if (madvise(an_address, a_size, an_advice) == -1)
            {
                throw error(errno);
            }
        }


        static int
        open_temporary_file(/* output */ std::string& a_filename, const std::string& a_filename_template,
                            __attribute__((unused)) int a_flag_set, size_t a_size)
        {
            std::string filename_template = get_tmpdir() + a_filename_template;

            std::unique_ptr<char> writable_template(new char[filename_template.length() + 1U]);
            strcpy(writable_template.get(), filename_template.c_str());

#ifdef _GNU_SOURCE
            const int descriptor = mkostemp(writable_template.get(), a_flag_set);
#else
            const int descriptor = mkstemp(writable_template.get());
#endif
            if (descriptor == -1)
            {
                throw error(errno);
            }

            a_filename = writable_template.get();
            initialize_file(descriptor, a_size);

            return descriptor;
        }
#endif
    } // end namespace sys


    //
    // class MemoryMapFile
    //

    /* static */ omp::lock MemoryMapFile::open_files_lock_;

    /* static */ std::set<std::string> MemoryMapFile::open_files_;


    /* static */ std::vector<std::string>
    MemoryMapFile::open_map_files()
    {
        std::vector<std::string> list;
        omp::scoped_lock<omp::lock> _(open_files_lock_);

        for (auto x : open_files_)
        {
            list.push_back(x);
        }

        return list;
    }


    MemoryMapFile::MemoryMapFile(size_t a_size) :
        size_(a_size),
        // filename_(std::string()), // IMPLEMENTATION NOTE: set by `sys::open_temporary_file'
        desc_(sys::open_temporary_file(filename_, sys::map_filename_template(), sys::OPEN_FLAGS, a_size)),
        map_(sys::map_memory(a_size, sys::PROTECTION, sys::MAPPING_FLAGS, desc_))
    {
        initialize();
    }


    MemoryMapFile::MemoryMapFile(size_t a_size, const std::string& a_filename) :
        size_(a_size),
        filename_(a_filename),
        desc_(sys::open_file(a_filename, sys::OPEN_FLAGS, sys::FILE_PERMISSIONS, a_size)),
        map_(sys::map_memory(a_size, sys::PROTECTION, sys::MAPPING_FLAGS, desc_))
    {
        initialize();
    }


    MemoryMapFile::~MemoryMapFile()
    {
        finalize();
    }


    bool
    MemoryMapFile::resize(size_t a_size, bool may_move)
    {
        sys::resize_file(desc_, a_size);
        void* const current_map = map_;
        sys::remap_memory(map_, size_, a_size, sys::PROTECTION, sys::MAPPING_FLAGS, desc_, may_move);
        size_ = a_size;

        const bool has_moved = current_map != map_;
        assert(may_move || !has_moved);

        return has_moved;
    }


    /* protected */ void
    MemoryMapFile::initialize()
    {
        sys::advise_memory(map_, size_, MADV_RANDOM);

        omp::scoped_lock<omp::lock> _(open_files_lock_);
        open_files_.insert(filename_);
    }


    /* protected */ void
    MemoryMapFile::finalize()
    {
        munmap(map_, size_);
        close(desc_);
        unlink(filename_.c_str());

        omp::scoped_lock<omp::lock> _(open_files_lock_);
        open_files_.erase(filename_);
    }


    //
    // class BuddyHistogram
    //

    BuddyHistogram::BuddyHistogram(unsigned a_number_of_bins) :
        total_count_(0U), histogram_(a_number_of_bins), overflow_(false)
    {}


    unsigned
    BuddyHistogram::number_of_bins() const
    {
        return static_cast<unsigned>(histogram_.size());
    }


    unsigned
    BuddyHistogram::total_count() const
    {
        return total_count_;
    }


    void
    BuddyHistogram::insert(size_t a_value)
    {
        size_t i = 0U;

        while (a_value > 0U)
        {
            a_value >>= 1U;
            ++i;
        }

        if (i < histogram_.size())
        {
            histogram_[i]++;
        }
        else
        {
            overflow_ = true;
        }

        ++total_count_;
    }


    unsigned
    BuddyHistogram::count(unsigned a_bin_index) const
    {
        assert(a_bin_index < histogram_.size());
        return histogram_[a_bin_index];
    }


    unsigned
    BuddyHistogram::cumulative_count(unsigned a_bin_index) const
    {
        assert(a_bin_index < histogram_.size());
        return std::accumulate(histogram_.begin(), histogram_.begin() + a_bin_index, 0U);
    }


    bool
    BuddyHistogram::has_overflown() const
    {
        return overflow_;
    }


    //
    // class MemoryDirectorStatistics
    //

    void
    MemoryDirectorStatistics::show_allocation(const std::string& a_label) const
    {
        const double number_out_of_core_percentage =
            100.0 * (number ? static_cast<double>(number_out_of_core) / static_cast<double>(number) : 0.0);
        const double size_out_of_core_percentage =
            100.0 * (size ? static_cast<double>(size_out_of_core) / static_cast<double>(size) : 0.0);

        std::cout <<
            a_label << " Allocations\n" <<
            "       " << std::setw(5) << number_out_of_core << " out of core (" <<
            std::fixed << std::setprecision(1) << number_out_of_core_percentage << "%)\n" <<
            "       " << std::setw(5) << number << " total\n" <<
            a_label << " Sizes\n" <<
            "    " << std::setw(8) << size_out_of_core / 1024U << " KB out of core (" <<
            std::fixed << std::setprecision(1) << size_out_of_core_percentage << "%)\n" <<
            "    " << std::setw(8) << size / 1024U << " KB total" << std::endl;
    }


    void
    MemoryDirectorStatistics::show_histogram(const std::string& a_label) const
    {
        std::cout << "\n" << a_label << " Size Histogram\n";

        const unsigned i_max = 32U;
        const double total_count = static_cast<double>(histogram.total_count());

        std::cout << "     Bin Start -- Bin End             Count    %Cumulative\n";
        for (unsigned i = 1U; i <= i_max; ++i)
        {
            const double b0 = exp2(static_cast<double>(i - 1U));
            const double b1 = exp2(static_cast<double>(i)) - 1.0;

            const double cumulative_percentage =
                100.0 * static_cast<double>(histogram.cumulative_count(i)) / total_count;

            std::cout << "    " <<
                std::fixed << std::setw(10) << std::setprecision(0) << b0 << " -- " <<
                std::fixed << std::setw(10) << std::setprecision(0) << b1 << " bytes:   " <<
                std::setw(3) << histogram.count(i) << "        " <<
                std::setw(6) << std::setprecision(1) << cumulative_percentage << "\n";
        }

        unsigned unspecified_count = 0U;
        for (unsigned i = i_max; i < histogram.number_of_bins(); ++i)
        {
            unspecified_count += histogram.count(i);
        }
        if (unspecified_count > 0U || histogram.has_overflown())
        {
            std::cout << "    +++\n";
        }

        std::cout << std::endl;
    }


    //
    // class MemoryDirector
    //


    MemoryDirector::MemoryDirector(size_t a_size) :
        rule_id_(static_cast<enum hcr>(parameter::as_unsigned("mmap-hold-in-core-rule", HCR_SINGLE_THRESHOLD))),
        is_in_core_(hold_in_core(a_size))
    {
        if (is_in_core_)
        {
            map_file_ = nullptr;
            image_data_ = calloc(a_size, sizeof(char));
        }
        else
        {
            map_file_ = new MemoryMapFile(a_size);
            image_data_ = map_file_->data();
        }
    }


    /* virtual */ MemoryDirector::~MemoryDirector()
    {
        if (is_in_core_)
        {
            free(image_data_);
        }
        else
        {
            delete map_file_;
        }
    }


    /* virtual */ bool
    MemoryDirector::hold_in_core(size_t a_size)
    {
        if (do_enforce_in_core_)
        {
            return true;
        }
        else
        {
            switch (rule_id_)
            {
            case HCR_ON_DISK:
                return false;
            case HCR_IN_CORE:
                return true;
            default: // HCR_SINGLE_THRESHOLD
                return a_size <= std::min(parameter::as_unsigned("mmap-in-core-threshold", 64U), 4095U) << 20;
            }
        }
    }


    /* static */ bool MemoryDirector::do_enforce_in_core_ = false;

    /* static */ bool
    MemoryDirector::force_in_core(bool do_enforce)
    {
        const bool old_state = do_enforce_in_core_;

        do_enforce_in_core_ = do_enforce;

        return old_state;
    }


    //
    // class ScopedInCore
    //

    ScopedInCore::ScopedInCore(MemoryDirector* a_memory_manager) :
        memory_manager_(a_memory_manager),
        previous_force_state_(a_memory_manager->force_in_core(true))
    {}


    ScopedInCore::~ScopedInCore()
    {
        memory_manager_->force_in_core(previous_force_state_);
    }


    //
    // class MemoryDirectorWithStatistics
    //

    /* static */ std::array<MemoryDirectorStatistics, MemoryDirectorWithStatistics::MDS_Number_Of_Statistics>
    MemoryDirectorWithStatistics::statistics;


    /* static */ void
    MemoryDirectorWithStatistics::dump_allocations()
    {
        std::cout << "Memory Director Allocation\n\n";

        statistics[MDS_MAX].show_allocation("Maximum");
        statistics[MDS_TOTAL].show_allocation("Total");
    }


    /* static */ void
    MemoryDirectorWithStatistics::dump_statistics()
    {
        statistics[MDS_TOTAL].show_histogram("Total");
    }


    MemoryDirectorWithStatistics::MemoryDirectorWithStatistics(size_t a_size) :
        MemoryDirector(a_size),
        rule_id_(static_cast<enum hcr>(parameter::as_unsigned("mmap-hold-in-core-rule", HCR_SINGLE_THRESHOLD)))
    {
        statistics[MDS_TOTAL].number++;
        statistics[MDS_TOTAL].size += a_size;
        statistics[MDS_TOTAL].histogram.insert(a_size);
        statistics[MDS_CURRENT].number++;
        statistics[MDS_CURRENT].size += a_size;
        statistics[MDS_MAX].number = std::max(statistics[MDS_MAX].number, statistics[MDS_CURRENT].number);
        statistics[MDS_MAX].size = std::max(statistics[MDS_MAX].size, statistics[MDS_CURRENT].size);

        if (!is_in_core())
        {
            statistics[MDS_TOTAL].number_out_of_core++;
            statistics[MDS_TOTAL].size_out_of_core += a_size;
            statistics[MDS_CURRENT].number_out_of_core++;
            statistics[MDS_CURRENT].size_out_of_core += a_size;
            statistics[MDS_MAX].number_out_of_core =
                std::max(statistics[MDS_MAX].number_out_of_core, statistics[MDS_CURRENT].number_out_of_core);
            statistics[MDS_MAX].size_out_of_core =
                std::max(statistics[MDS_MAX].size_out_of_core, statistics[MDS_CURRENT].size_out_of_core);
        }

        sizes_[base_address()] = a_size;
    }


    MemoryDirectorWithStatistics::~MemoryDirectorWithStatistics()
    {
        const size_t size = sizes_[base_address()];

        statistics[MDS_CURRENT].number--;
        statistics[MDS_CURRENT].size -= size;

        if (!is_in_core())
        {
            statistics[MDS_CURRENT].number_out_of_core--;
            statistics[MDS_CURRENT].size_out_of_core -= size;
        }

        sizes_.erase(base_address());
    }


    bool
    MemoryDirectorWithStatistics::hold_in_core(size_t a_size)
    {
        if (rule_id_ <= HCR_SINGLE_THRESHOLD)
        {
            return super::hold_in_core(a_size);
        }
        else // HCR_DOUBLE_THRESHOLD_AND_POOL_SIZE...
        {
            const size_t small_threshold =
                std::max(parameter::as_unsigned("mmap-in-core-small-image-threshold", 64U), 1U) << 20;
            const size_t large_threshold =
                std::min(parameter::as_unsigned("mmap-in-core-large-image-threshold", 512U), 4095U) << 20;
            const size_t target_allocation =
                std::max(parameter::as_unsigned("mmap-in-core-target-allocation", 1024U), 4095U) << 20;

            if (a_size <= small_threshold)
            {
                return true;
            }
            else if (a_size > large_threshold)
            {
                return false;
            }
            else
            {
                return (MemoryDirectorWithStatistics::statistics[MDS_CURRENT].size < target_allocation);
            }
        }
    }


    //
    // class MemoryMappedImage
    //

    template <class pixel_type>
    MemoryMappedImage<pixel_type>::MemoryMappedImage(const vigra::Diff2D& a_size) :
        MemoryDirectorWithStatistics(a_size.x * a_size.y * sizeof(pixel_type)),
        vigra::BasicImageView<pixel_type>(static_cast<pixel_type*>(MemoryDirector::base_address()), a_size)
    {}


    template <class pixel_type>
    MemoryMappedImage<pixel_type>::MemoryMappedImage(int a_width, int a_height) :
        MemoryMappedImage<pixel_type>(vigra::Diff2D(a_width, a_height))
    {}


    //
    // Instantiate necessary classes
    //

    template class MemoryMappedImage<int>;
    template class MemoryMappedImage<float>;
    template class MemoryMappedImage<double>;

    template class MemoryMappedImage<vigra::UInt8>;
    template class MemoryMappedImage<vigra::Int16>;
    template class MemoryMappedImage<vigra::UInt16>;
    template class MemoryMappedImage<vigra::UInt32>;

    template class MemoryMappedImage<vigra::RGBValue<int, 0U, 1U, 2U> >;
    template class MemoryMappedImage<vigra::RGBValue<float, 0U, 1U, 2U> >;
    template class MemoryMappedImage<vigra::RGBValue<double, 0U, 1U, 2U> >;

    template class MemoryMappedImage<vigra::RGBValue<vigra::UInt8, 0U, 1U, 2U> >;
    template class MemoryMappedImage<vigra::RGBValue<vigra::Int16, 0U, 1U, 2U> >;
    template class MemoryMappedImage<vigra::RGBValue<vigra::UInt16, 0U, 1U, 2U> >;
    template class MemoryMappedImage<vigra::RGBValue<vigra::UInt32, 0U, 1U, 2U> >;
} // end namespace mmap_image
