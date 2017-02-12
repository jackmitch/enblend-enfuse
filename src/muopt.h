/*
 * Copyright (C) 2012-2017 Dr. Christoph L. Spiel
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
#ifndef MUOPT_H_INCLUDED_
#define MUOPT_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#define PRAGMA(m_token_sequence) _Pragma(#m_token_sequence)


#ifdef __GNUC__

#define PREFETCH(m_addr) __builtin_prefetch(m_addr)
#define HINTED_PREFETCH(m_addr, m_rw_hint, m_temporal_locality_hint) \
    __builtin_prefetch((m_addr), (m_rw_hint), (m_temporal_locality_hint))

#define EXPECT_RESULT(m_condition, m_expected_result) \
    __builtin_expect((m_condition), static_cast<int>(m_expected_result))


#ifdef __ICC
// See e.g.
//     https://software.intel.com/en-us/articles/data-alignment-to-assist-vectorization
#define ASSUME_ALIGNED_F_(m_pointer_expression, m_alignment) \
    (__assume_aligned(m_pointer_expression, m_alignment), m_pointer_expression)
#else
#define ASSUME_ALIGNED_F_(m_pointer_expression, m_alignment) \
    __builtin_assume_aligned(m_pointer_expression, m_alignment)
#endif // __ICC


#ifdef DEBUG
#define ASSUME_ALIGNED(m_pointer_expression, m_alignment) \
    (assert((reinterpret_cast<size_t>(m_pointer_expression) & ((m_alignment) - 1U)) == size_t()), \
     static_cast<decltype (m_pointer_expression)> \
     (ASSUME_ALIGNED_F_(m_pointer_expression, m_alignment)))
#else
#define ASSUME_ALIGNED(m_pointer_expression, m_alignment) \
    static_cast<decltype (m_pointer_expression)> \
    (ASSUME_ALIGNED_F_(m_pointer_expression, m_alignment))
#endif // DEBUG


#ifdef __ICC
#define ASSUME_NO_VECTOR_DEPENDENCY PRAGMA(ivdep)
#else
#define ASSUME_NO_VECTOR_DEPENDENCY PRAGMA(GCC ivdep)
#endif // __ICC


#ifdef __clang__
#define VECTORIZE_LOOP PRAGMA(clang loop vectorize(enable) interleave(enable))
#define VECTORIZE_LOOP_USING(...) PRAGMA(clang loop __VA_ARGS__)
#define NO_VECTORIZE_LOOP PRAGMA(clang loop vectorize(disable) interleave(disable))
#else
#define VECTORIZE_LOOP PRAGMA(vector)
#define VECTORIZE_LOOP_USING(...) PRAGMA(vector __VA_ARGS__)
#define NO_VECTORIZE_LOOP PRAGMA(novector)
#endif // __clang__


#if defined(__ICC)
#define ASSUME_F_(m_condition) __assume(m_condition)
#elif defined(__clang__)
// See e.g.
//     http://clang.llvm.org/docs/LanguageExtensions.html
#if __has_builtin(__builtin_assume)
#define ASSUME_F_(m_condition) __builtin_assume(m_condition)
#else
#define ASSUME_F_(m_condition)
#endif // __has_builtin
#else
#define ASSUME_F_(m_condition)
#endif

#ifdef DEBUG
#define ASSUME(m_condition) (assert(m_condition), ASSUME_F_(m_condition))
#else
#define ASSUME(m_condition) ASSUME_F_(m_condition)
#endif // DEBUG


#else // neither GCC family, CLang family, nor ICC


#define PREFETCH(m_addr)
#define HINTED_PREFETCH(m_addr, m_rw_hint, m_temporal_locality_hint)

#define EXPECT_RESULT(m_condition, m_expected_result) (m_condition)

#ifdef DEBUG
#define ASSUME(m_condition) assert(m_condition)
#define ASSUME_ALIGNED(m_pointer_expression, m_alignment) \
    (assert((reinterpret_cast<size_t>(m_pointer_expression) & ((m_alignment) - 1U)) == size_t()), \
     (m_pointer_expression))
#else
#define ASSUME(m_condition)
#define ASSUME_ALIGNED(m_pointer_expression, m_alignment) (m_pointer_expression)
#endif // DEBUG

#define ASSUME_NO_VECTOR_DEPENDENCY

#define VECTORIZE_LOOP
#define VECTORIZE_LOOP_USING(...)
#define NO_VECTORIZE_LOOP


#endif // __GNUC__


typedef enum {
    PREPARE_FOR_READ,
    PREPARE_FOR_WRITE
} rw_hint;


// If data is only touched once, or if the dataset is smaller than the
// cache, prefer the non-temporal version; otherwise use one of the
// temporal versions.
typedef enum {
    // Fetch data into the first way of the L1/L2 cache, minimizing cache pollution.
    NO_TEMPORAL_LOCALITY,

    // Fetch data into the least-recently-used way of the ...
    LOW_TEMPORAL_LOCALITY,      // ... L3 cache?
    MEDIUM_TEMPORAL_LOCALITY,   // ... L2/L3 cache?
    HIGH_TEMPORAL_LOCALITY      // ... L1/L2/L3 cache just as a normal load would do.
} temporal_locality_hint;


#endif // MUOPT_H_INCLUDED_

// Local Variables:
// mode: c++
// End:
