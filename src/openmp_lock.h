/*
 * Copyright (C) 2013-2017 Christoph L. Spiel
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
#ifndef OPENMP_LOCK_H_INCLUDED_
#define OPENMP_LOCK_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "openmp_def.h"


namespace omp
{
#ifdef OPENMP

    class lock
    {
    public:
        lock() {omp_init_lock(&lock_);}
        ~lock() {omp_destroy_lock(&lock_);}

        lock(const lock&) = delete;
        lock& operator=(const lock&) = delete;

        void set() {omp_set_lock(&lock_);}
        void unset() {omp_unset_lock(&lock_);}
        bool test() {return omp_test_lock(&lock_) != 0;}

    private:
        omp_lock_t lock_;
    };


    class nestable_lock
    {
    public:
        nestable_lock() {omp_init_nest_lock(&lock_);}
        ~nestable_lock() {omp_destroy_nest_lock(&lock_);}

        nestable_lock(const nestable_lock&) = delete;
        nestable_lock& operator=(const nestable_lock&) = delete;

        void set() {omp_set_nest_lock(&lock_);}
        void unset() {omp_unset_nest_lock(&lock_);}
        bool test() {return omp_test_nest_lock(&lock_) != 0;}

    private:
        omp_nest_lock_t lock_;
    };

#else

    class lock
    {
    public:
        lock() {}
        ~lock() {}

        lock(const lock&) = delete;
        lock& operator=(const lock&) = delete;

        void set() {}
        void unset() {}
        bool test() {return false;}
    };


    class nestable_lock
    {
    public:
        nestable_lock() {}
        ~nestable_lock() {}

        nestable_lock(const nestable_lock&) = delete;
        nestable_lock& operator=(const nestable_lock&) = delete;

        void set() {}
        void unset() {}
        bool test() {return false;}
    };

#endif // OPENMP


    template <class basic_lock>
    class scoped_lock
    {
    public:
        scoped_lock() = delete;

        explicit scoped_lock(basic_lock& a_lock) : lock_(&a_lock) {lock_->set();}
        ~scoped_lock() {lock_->unset();}

        scoped_lock(const scoped_lock&) = delete;
        scoped_lock& operator=(const scoped_lock&) = delete;

        void set() {lock_->set();}
        void unset() {lock_->unset();}
        bool test() {return lock_->test();}

    private:
        basic_lock* const lock_;
    }; // class scoped_lock
} // namespace omp


#endif // OPENMP_LOCK_H_INCLUDED_

// Local Variables:
// mode: c++
// End:
