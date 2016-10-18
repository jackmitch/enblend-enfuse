/*
 * Copyright (C) 2016 Dr. Christoph L. Spiel
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

#ifndef ALLOCATE_H_INLCUDED
#define ALLOCATE_H_INLCUDED


#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace allocate
{
    namespace detail
    {
        class positional_hint
        {
        public:
            positional_hint() = delete;
            explicit positional_hint(std::size_t a_position) : position_(a_position) {}

            std::size_t position() const {return position_;}

            const char* append_location_to(const char* a_message) const
            {
                std::stringstream message;
                message << a_message << " at position " << position();
                return message.str().c_str();
            }

            const char* append_location_to(const std::string& a_message) const
            {
                return append_location_to(a_message.c_str());
            }

        private:
            const std::size_t position_;
        }; // class positional_hint
    } // namespace detail


    class not_initialized : public std::exception, private detail::positional_hint
    {
    public:
        not_initialized() = delete;
        explicit not_initialized(std::size_t an_uninitialized_position) :
            detail::positional_hint(an_uninitialized_position) {}

        virtual const char* what() const noexcept
        {
            return detail::positional_hint::append_location_to("uninitialized array element");
        }

        std::size_t position() const {return detail::positional_hint::position();}
    }; // class not_initialized


    class already_initialized : public std::exception, private detail::positional_hint
    {
    public:
        already_initialized() = delete;
        explicit already_initialized(std::size_t an_uninitialized_position) :
            detail::positional_hint(an_uninitialized_position) {}

        virtual const char* what() const noexcept
        {
            return detail::positional_hint::append_location_to("array element already initialized");
        }

        std::size_t position() const {return detail::positional_hint::position();}
    }; // class already_initialized


    class out_of_range : public std::exception, private detail::positional_hint
    {
    public:
        out_of_range() = delete;
        out_of_range(std::size_t an_out_of_range_position, std::size_t a_size) :
            detail::positional_hint(an_out_of_range_position), array_size_(a_size) {}

        virtual const char* what() const noexcept
        {
            std::stringstream message;
            message << "out of range (0.." << (array_size() - 1) << ") access";
            return detail::positional_hint::append_location_to(message.str());
        }

        std::size_t position() const {return detail::positional_hint::position();}
        std::size_t array_size() const {return array_size_;}

    private:
        const std::size_t array_size_;
    }; // class out_of_range


    template <typename t>
    struct default_deallocator
    {
        typedef t value_type;

        default_deallocator() {}
        void operator()(value_type*) const {}
    }; // template struct default_deallocator


    template <typename t, typename deallocator = default_deallocator<t>>
    class array
    {
    public:
        typedef t value_type;
        typedef deallocator deallocator_type;

        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

        typedef value_type& reference;
        typedef const value_type& const_reference;
        typedef value_type* pointer;
        typedef const value_type* const_pointer;

        typedef t* iterator;
        typedef const t* const_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        array() = delete;
        array(const array&) = delete;
        array& operator=(const array&) = delete;

        explicit array(size_t an_object_count) :
            size_(an_object_count),
            base_(::operator new[](an_object_count * sizeof(t))),
            initialized_(an_object_count)
        {
            initialization_order_.reserve(an_object_count);
        }

        ~array()
        {
            finalize();
            ::operator delete[](base_);
        }

        // Not implemented yet:
        //     reference operator[](size_type a_position) {...}

        const_reference operator[](size_type a_position) const
        {
            if (is_initialized(a_position))
            {
                return *(begin() + a_position);
            }
            else
            {
                throw not_initialized(a_position);
            }
        }

        // Not implemented yet:
        //     reference at(size_type a_position) {...}

        const_reference at(size_type a_position) const
        {
            if (a_position >= size_)
            {
                throw out_of_range(a_position, size());
            }
            else
            {
                return operator[](a_position);
            }
        }

        bool empty() const {return size_ == size_type();}
        size_type size() const {return size_;}

        iterator begin() {return iterator(base_);}
        iterator end() {return begin() + size_;}
        const_iterator begin() const {return iterator(base_);}
        const_iterator end() const {return begin() + size_;}

        reverse_iterator rbegin() {return reverse_iterator(end());}
        reverse_iterator rend() {return reverse_iterator(begin());}
        const_reverse_iterator rbegin() const {return const_reverse_iterator(end());}
        const_reverse_iterator rend() const {return const_reverse_iterator(begin());}

        void mark_as_initialized(size_type a_position)
        {
            if (is_initialized(a_position))
            {
                throw already_initialized(a_position);
            }
            else
            {
                initialized_.at(a_position) = true;
                initialization_order_.push_back(a_position);
            }
        }

        void mark_as_initialized(iterator an_iterator) {mark_as_initialized(an_iterator - begin());}

        bool is_initialized(size_type a_position) const {return initialized_.at(a_position);}
        bool is_initialized(iterator an_iterator) const {return is_initialized(an_iterator - begin());}

    private:
        void finalize()
        {
            deallocator_type deallocate;

            // Implementation Note: We apply the destructors in the
            // reverse order of construction.  This is the sole reason
            // why we have `initialization_order_'.
            for (auto position = initialization_order_.rbegin();
                 position != initialization_order_.rend();
                 ++position)
            {
                const iterator value(begin() + *position);
                value->~t();
                deallocate(value);
            }
        }

        const size_type size_;
        void* const base_;
        std::vector<bool> initialized_;
        std::vector<size_type> initialization_order_;
    }; // template class array
} // namespace allocate


#endif // ALLOCATE_H_INLCUDED

// Local Variables:
// mode: c++
// End:
