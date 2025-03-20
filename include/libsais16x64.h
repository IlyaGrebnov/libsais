/*--

This file is a part of libsais, a library for linear time suffix array,
longest common prefix array and burrows wheeler transform construction.

   Copyright (c) 2021-2025 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information.

--*/

#ifndef LIBSAIS16X64_H
#define LIBSAIS16X64_H 1

#define LIBSAIS16X64_VERSION_MAJOR   2
#define LIBSAIS16X64_VERSION_MINOR   9
#define LIBSAIS16X64_VERSION_PATCH   1
#define LIBSAIS16X64_VERSION_STRING  "2.9.1"

#ifdef _WIN32
    #ifdef LIBSAIS_SHARED
        #ifdef LIBSAIS_EXPORTS
            #define LIBSAIS16X64_API __declspec(dllexport)
        #else
            #define LIBSAIS16X64_API __declspec(dllimport)
        #endif
    #else
        #define LIBSAIS16X64_API
    #endif
#else
    #define LIBSAIS16X64_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

    #include <stdint.h>

    /**
    * Constructs the suffix array of a given 16-bit string.
    * @param T [0..n-1] The input 16-bit string.
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the given 16-bit string.
    * @param fs The extra space available at the end of SA array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64(const uint16_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq);

    /**
    * Constructs the generalized suffix array (GSA) of a given 16-bit string set.
    * @param T [0..n-1] The input 16-bit string set using 0 as separators (T[n-1] must be 0).
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the given 16-bit string set.
    * @param fs The extra space available at the end of SA array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_gsa(const uint16_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq);

    /**
    * Constructs the suffix array of a given integer array.
    * Note, during construction input array will be modified, but restored at the end if no errors occurred.
    * @param T [0..n-1] The input integer array.
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the integer array.
    * @param k The alphabet size of the input integer array.
    * @param fs Extra space available at the end of SA array (can be 0, but 4k or better 6k is recommended for optimal performance).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_long(int64_t * T, int64_t * SA, int64_t n, int64_t k, int64_t fs);

#if defined(LIBSAIS_OPENMP)
    /**
    * Constructs the suffix array of a given 16-bit string in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string.
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the given 16-bit string.
    * @param fs The extra space available at the end of SA array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_omp(const uint16_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq, int64_t threads);

    /**
    * Constructs the generalized suffix array (GSA) of a given 16-bit string set in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string set using 0 as separators (T[n-1] must be 0).
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the given 16-bit string set.
    * @param fs The extra space available at the end of SA array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_gsa_omp(const uint16_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq, int64_t threads);

    /**
    * Constructs the suffix array of a given integer array in parallel using OpenMP.
    * Note, during construction input array will be modified, but restored at the end if no errors occurred.
    * @param T [0..n-1] The input integer array.
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the integer array.
    * @param k The alphabet size of the input integer array.
    * @param fs Extra space available at the end of SA array (can be 0, but 4k or better 6k is recommended for optimal performance).
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_long_omp(int64_t * T, int64_t * SA, int64_t n, int64_t k, int64_t fs, int64_t threads);
#endif

    /**
    * Constructs the burrows-wheeler transformed 16-bit string (BWT) of a given 16-bit string.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n-1+fs] The temporary array.
    * @param n The length of the given 16-bit string.
    * @param fs The extra space available at the end of A array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @return The primary index if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_bwt(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, int64_t fs, int64_t * freq);

    /**
    * Constructs the burrows-wheeler transformed 16-bit string (BWT) of a given 16-bit string with auxiliary indexes.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n-1+fs] The temporary array.
    * @param n The length of the given 16-bit string.
    * @param fs The extra space available at the end of A array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @param r The sampling rate for auxiliary indexes (must be power of 2).
    * @param I [0..(n-1)/r] The output auxiliary indexes.
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_bwt_aux(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, int64_t fs, int64_t * freq, int64_t r, int64_t * I);

#if defined(LIBSAIS_OPENMP)
    /**
    * Constructs the burrows-wheeler transformed 16-bit string (BWT) of a given 16-bit string in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n-1+fs] The temporary array.
    * @param n The length of the given 16-bit string.
    * @param fs The extra space available at the end of A array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return The primary index if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_bwt_omp(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, int64_t fs, int64_t * freq, int64_t threads);

    /**
    * Constructs the burrows-wheeler transformed 16-bit string (BWT) of a given 16-bit string with auxiliary indexes in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n-1+fs] The temporary array.
    * @param n The length of the given 16-bit string.
    * @param fs The extra space available at the end of A array (0 should be enough for most cases).
    * @param freq [0..65535] The output 16-bit symbol frequency table (can be NULL).
    * @param r The sampling rate for auxiliary indexes (must be power of 2).
    * @param I [0..(n-1)/r] The output auxiliary indexes.
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_bwt_aux_omp(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, int64_t fs, int64_t * freq, int64_t r, int64_t * I, int64_t threads);
#endif

    /**
    * Constructs the original 16-bit string from a given burrows-wheeler transformed 16-bit string (BWT) with primary index.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n] The temporary array (NOTE, temporary array must be n + 1 size).
    * @param n The length of the given 16-bit string.
    * @param freq [0..65535] The input 16-bit symbol frequency table (can be NULL).
    * @param i The primary index.
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_unbwt(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, const int64_t * freq, int64_t i);

    /**
    * Constructs the original 16-bit string from a given burrows-wheeler transformed 16-bit string (BWT) with auxiliary indexes.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n] The temporary array (NOTE, temporary array must be n + 1 size).
    * @param n The length of the given 16-bit string.
    * @param freq [0..65535] The input 16-bit symbol frequency table (can be NULL).
    * @param r The sampling rate for auxiliary indexes (must be power of 2).
    * @param I [0..(n-1)/r] The input auxiliary indexes.
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_unbwt_aux(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, const int64_t * freq, int64_t r, const int64_t * I);

#if defined(LIBSAIS_OPENMP)
    /**
    * Constructs the original 16-bit string from a given burrows-wheeler transformed 16-bit string (BWT) with primary index in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n] The temporary array (NOTE, temporary array must be n + 1 size).
    * @param n The length of the given 16-bit string.
    * @param freq [0..65535] The input 16-bit symbol frequency table (can be NULL).
    * @param i The primary index.
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_unbwt_omp(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, const int64_t * freq, int64_t i, int64_t threads);

    /**
    * Constructs the original 16-bit string from a given burrows-wheeler transformed 16-bit string (BWT) with auxiliary indexes in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string.
    * @param U [0..n-1] The output 16-bit string (can be T).
    * @param A [0..n] The temporary array (NOTE, temporary array must be n + 1 size).
    * @param n The length of the given 16-bit string.
    * @param freq [0..65535] The input 16-bit symbol frequency table (can be NULL).
    * @param r The sampling rate for auxiliary indexes (must be power of 2).
    * @param I [0..(n-1)/r] The input auxiliary indexes.
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_unbwt_aux_omp(const uint16_t * T, uint16_t * U, int64_t * A, int64_t n, const int64_t * freq, int64_t r, const int64_t * I, int64_t threads);
#endif

    /**
    * Constructs the permuted longest common prefix array (PLCP) of a given 16-bit string and a suffix array.
    * @param T [0..n-1] The input 16-bit string.
    * @param SA [0..n-1] The input suffix array.
    * @param PLCP [0..n-1] The output permuted longest common prefix array.
    * @param n The length of the 16-bit string and the suffix array.
    * @return 0 if no error occurred, -1 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_plcp(const uint16_t * T, const int64_t * SA, int64_t * PLCP, int64_t n);

    /**
    * Constructs the permuted longest common prefix array (PLCP) of a given 16-bit string set and a generalized suffix array (GSA).
    * @param T [0..n-1] The input 16-bit string set using 0 as separators (T[n-1] must be 0).
    * @param SA [0..n-1] The input generalized suffix array.
    * @param PLCP [0..n-1] The output permuted longest common prefix array.
    * @param n The length of the string set and the generalized suffix array.
    * @return 0 if no error occurred, -1 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_plcp_gsa(const uint16_t * T, const int64_t * SA, int64_t * PLCP, int64_t n);

    /**
    * Constructs the longest common prefix array (LCP) of a given permuted longest common prefix array (PLCP) and a suffix array.
    * @param PLCP [0..n-1] The input permuted longest common prefix array.
    * @param SA [0..n-1] The input suffix array or generalized suffix array (GSA).
    * @param LCP [0..n-1] The output longest common prefix array (can be SA).
    * @param n The length of the permuted longest common prefix array and the suffix array.
    * @return 0 if no error occurred, -1 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_lcp(const int64_t * PLCP, const int64_t * SA, int64_t * LCP, int64_t n);

#if defined(LIBSAIS_OPENMP)
    /**
    * Constructs the permuted longest common prefix array (PLCP) of a given 16-bit string and a suffix array in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string.
    * @param SA [0..n-1] The input suffix array.
    * @param PLCP [0..n-1] The output permuted longest common prefix array.
    * @param n The length of the 16-bit string and the suffix array.
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_plcp_omp(const uint16_t * T, const int64_t * SA, int64_t * PLCP, int64_t n, int64_t threads);

    /**
    * Constructs the permuted longest common prefix array (PLCP) of a given 16-bit string set and a generalized suffix array (GSA) in parallel using OpenMP.
    * @param T [0..n-1] The input 16-bit string set using 0 as separators (T[n-1] must be 0).
    * @param SA [0..n-1] The input generalized suffix array.
    * @param PLCP [0..n-1] The output permuted longest common prefix array.
    * @param n The length of the string set and the generalized suffix array.
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_plcp_gsa_omp(const uint16_t * T, const int64_t * SA, int64_t * PLCP, int64_t n, int64_t threads);

    /**
    * Constructs the longest common prefix array (LCP) of a given permuted longest common prefix array (PLCP) and a suffix array in parallel using OpenMP.
    * @param PLCP [0..n-1] The input permuted longest common prefix array.
    * @param SA [0..n-1] The input suffix array or generalized suffix array (GSA).
    * @param LCP [0..n-1] The output longest common prefix array (can be SA).
    * @param n The length of the permuted longest common prefix array and the suffix array.
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 otherwise.
    */
    LIBSAIS16X64_API int64_t libsais16x64_lcp_omp(const int64_t * PLCP, const int64_t * SA, int64_t * LCP, int64_t n, int64_t threads);
#endif

#ifdef __cplusplus
}
#endif

#endif
