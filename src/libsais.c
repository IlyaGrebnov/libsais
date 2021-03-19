/*--

This file is a part of libsais, a library for linear time
suffix array and burrows wheeler transform construction.

   Copyright (c) 2021 Ilya Grebnov <ilya.grebnov@gmail.com>

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

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "libsais.h"

#define INT_BIT                         (32)
#define ALPHABET_SIZE                   (1 << CHAR_BIT)
#define SUFFIX_GROUP_BIT                (INT_BIT - 1)
#define SUFFIX_GROUP_MARKER             (1 << (SUFFIX_GROUP_BIT - 1))

#define BUCKETS_INDEX2(_c, _s)          (((_c) << 1) + (_s))
#define BUCKETS_INDEX4(_c, _s)          (((_c) << 2) + (_s))

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
    #define FORCEINLINE inline __attribute__((__always_inline__))
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
    #define RESTRICT __restrict
    #define FORCEINLINE inline __forceinline
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if defined(__has_builtin)
    #if __has_builtin(__builtin_prefetch)
        #define HAS_BUILTIN_PREFECTCH
    #endif
#elif defined(__GNUC__) && __GNUC__ > 3
    #define HAS_BUILTIN_PREFECTCH
#endif 

#if defined(HAS_BUILTIN_PREFECTCH)
    #define libsais_prefetch(address) __builtin_prefetch((const void *)(address), 0, 0)
    #define libsais_prefetchw(address) __builtin_prefetch((const void *)(address), 1, 0)
#elif defined (_M_IX86) || defined (_M_AMD64)
    #include <intrin.h>
    #define libsais_prefetch(address) _mm_prefetch((const void *)(address), _MM_HINT_NTA)
    #define libsais_prefetchw(address) _m_prefetchw((const void *)(address))
#elif defined (_M_ARM)
    #include <intrin.h>
    #define libsais_prefetch(address) __prefetch((const void *)(address))
    #define libsais_prefetchw(address) __prefetchw((const void *)(address))
#elif defined (_M_ARM64)
    #include <intrin.h>
    #define libsais_prefetch(address) __prefetch2((const void *)(address), 1)
    #define libsais_prefetchw(address) __prefetch2((const void *)(address), 17)
#else
    #error Your compiler, configuration or platform is not supported.
#endif

static FORCEINLINE void * libsais_align_up(const void * address, size_t alignment)
{
    return (void *)((((intptr_t)address) + ((intptr_t)alignment) - 1) & (-((intptr_t)alignment)));
}

static FORCEINLINE void * libsais_aligned_malloc(size_t size, size_t alignment)
{
    void * address = malloc(size + sizeof(short) + alignment - 1);
    if (address != NULL)
    {
        void * aligned_address = libsais_align_up((void *)((intptr_t)address + (intptr_t)(sizeof(short))), alignment);
        ((short *)aligned_address)[-1] = (short)((intptr_t)aligned_address - (intptr_t)address);

        return aligned_address;
    }

    return NULL;
}

static FORCEINLINE void libsais_aligned_free(void * aligned_address)
{
    if (aligned_address != NULL)
    {
        free((void *)((intptr_t)aligned_address - ((short *)aligned_address)[-1]));
    }
}

static int libsais_gather_lms_suffixes_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n)
{
    const ptrdiff_t prefetch_distance = 128;

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= 3; i -= 4)
    {
        libsais_prefetch(&T[i - prefetch_distance]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((s & 3) == 1);
        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((s & 3) == 1);
        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((s & 3) == 1);
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
    }

    return n - 1 - m;
}

static int libsais_gather_lms_suffixes_32s(const int * RESTRICT T, int * RESTRICT SA, int n)
{
    const ptrdiff_t prefetch_distance = 32;

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= 3; i -= 4)
    {
        libsais_prefetch(&T[i - prefetch_distance]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((s & 3) == 1);
        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((s & 3) == 1);
        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((s & 3) == 1);
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
    }

    return n - 1 - m;
}

static int libsais_gather_compacted_lms_suffixes_32s(const int * RESTRICT T, int * RESTRICT SA, int n)
{
    const ptrdiff_t prefetch_distance = 32;

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= 3; i -= 4)
    {
        libsais_prefetch(&T[i - prefetch_distance]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((ptrdiff_t)(s & 3) == (c0 >= 0));
        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((ptrdiff_t)(s & 3) == (c1 >= 0));
        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((ptrdiff_t)(s & 3) == (c0 >= 0));
        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((ptrdiff_t)(s & 3) == (c1 >= 0));
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((ptrdiff_t)(s & 3) == (c1 >= 0));
    }

    return n - 1 - m;
}

static void libsais_count_lms_suffixes_32s_2k(const int * RESTRICT T, int n, int k, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    memset(buckets, 0, 2 * (size_t)k * sizeof(int));

    int             i   = n - 2;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetch(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3], 0)]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1)));
        buckets[BUCKETS_INDEX2((size_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1)));
        buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;

        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1)));
        buckets[BUCKETS_INDEX2((size_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1)));
        buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1)));
        buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;
    }

    buckets[BUCKETS_INDEX2((size_t)c0, 0)]++;
}

static int libsais_count_and_gather_lms_suffixes_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 128;

    memset(buckets, 0, 4 * ALPHABET_SIZE * sizeof(int));

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= 3; i -= 4)
    {
        libsais_prefetch(&T[i - prefetch_distance]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c0, s & 3)]++;

        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]++;

        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c0, s & 3)]++;

        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]++;
    }

    buckets[BUCKETS_INDEX4((size_t)c0, (s << 1) & 3)]++;

    return n - 1 - m;
}

static int libsais_count_and_gather_lms_suffixes_32s_4k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    memset(buckets, 0, 4 * (size_t)k * sizeof(int));

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetch(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 0], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 1], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 2], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 3], 0)]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c0, s & 3)]++;

        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]++;

        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c0, s & 3)]++;

        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]++;
    }

    buckets[BUCKETS_INDEX4((size_t)c0, (s << 1) & 3)]++;

    return n - 1 - m;
}

static int libsais_count_and_gather_lms_suffixes_32s_2k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    memset(buckets, 0, 2 * (size_t)k * sizeof(int));

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetch(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3], 0)]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX2((size_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;

        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX2((size_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((s & 3) == 1);
        buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;
    }

    buckets[BUCKETS_INDEX2((size_t)c0, 0)]++;

    return n - 1 - m;
}

static int libsais_count_and_gather_compacted_lms_suffixes_32s_2k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    memset(buckets, 0, 2 * (size_t)k * sizeof(int));

    int             i   = n - 2;
    int             m   = n - 1;
    size_t          s   = 1;
    ptrdiff_t       c0  = T[n - 1];
    ptrdiff_t       c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetch(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0] & INT_MAX, 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1] & INT_MAX, 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2] & INT_MAX, 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3] & INT_MAX, 0)]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((ptrdiff_t)(s & 3) == (c0 >=0));
        c0 &= INT_MAX; buckets[BUCKETS_INDEX2((size_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 0; m -= ((ptrdiff_t)(s & 3) == (c1 >= 0));
        c1 &= INT_MAX; buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;

        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); SA[m] = i - 1; m -= ((ptrdiff_t)(s & 3) == (c0 >= 0));
        c0 &= INT_MAX; buckets[BUCKETS_INDEX2((size_t)c0, (s & 3) == 1)]++;

        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i - 2; m -= ((ptrdiff_t)(s & 3) == (c1 >= 0));
        c1 &= INT_MAX; buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); SA[m] = i + 1; m -= ((ptrdiff_t)(s & 3) == (c1 >= 0));
        c1 &= INT_MAX; buckets[BUCKETS_INDEX2((size_t)c1, (s & 3) == 1)]++;
    }

    c0 &= INT_MAX; buckets[BUCKETS_INDEX2((size_t)c0, 0)]++;

    return n - 1 - m;
}

static void libsais_count_suffixes_32s(const int * RESTRICT T, int n, int k, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    memset(buckets, 0, (size_t)k * sizeof(int));

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - 7; i < j; i += 8)
    {
        libsais_prefetch(&T[i + prefetch_distance]);

        buckets[T[i + 0]]++;
        buckets[T[i + 1]]++;
        buckets[T[i + 2]]++;
        buckets[T[i + 3]]++;
        buckets[T[i + 4]]++;
        buckets[T[i + 5]]++;
        buckets[T[i + 6]]++;
        buckets[T[i + 7]]++;
    }

    for (j += 7; i < j; i += 1)
    {
        buckets[T[i]]++;
    }
}

static void libsais_initialize_buckets_start_and_end_8u(int * RESTRICT buckets)
{
    int * RESTRICT bucket_start = &buckets[6 * ALPHABET_SIZE];
    int * RESTRICT bucket_end   = &buckets[7 * ALPHABET_SIZE];

    ptrdiff_t i, j; int sum = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4(UCHAR_MAX, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
    {
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];
        bucket_end[j] = sum;
    }
}

static void libsais_initialize_buckets_start_and_end_32s_6k(int k, int * RESTRICT buckets)
{
    int * RESTRICT bucket_start = &buckets[4 * k];
    int * RESTRICT bucket_end   = &buckets[5 * k];

    ptrdiff_t i, j; int sum = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
    {
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];
        bucket_end[j] = sum;
    }
}

static void libsais_initialize_buckets_start_and_end_32s_4k(int k, int * RESTRICT buckets)
{
    int * RESTRICT bucket_start = &buckets[2 * k];
    int * RESTRICT bucket_end   = &buckets[3 * k];

    ptrdiff_t i, j; int sum = 0;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    { 
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        bucket_end[j] = sum;
    }
}

static void libsais_initialize_buckets_end_32s_2k(int k, int * RESTRICT buckets)
{
    ptrdiff_t i; int sum0 = 0;
    for (i = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX2((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0))
    { 
        sum0 += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)]; buckets[i + BUCKETS_INDEX2(0, 0)] = sum0;
    }
}

static void libsais_initialize_buckets_start_32s_1k(int k, int * RESTRICT buckets)
{
    ptrdiff_t i; int sum = 0;
    for (i = 0; i <= (ptrdiff_t)k - 1; i += 1) { int tmp = buckets[i]; buckets[i] = sum; sum += tmp; }
}

static void libsais_initialize_buckets_end_32s_1k(int k, int * RESTRICT buckets)
{
    ptrdiff_t i; int sum = 0;
    for (i = 0; i <= (ptrdiff_t)k - 1; i += 1) { sum += buckets[i]; buckets[i] = sum; }
}

static void libsais_initialize_buckets_start_and_end_32s_2k(int k, int * RESTRICT buckets)
{
    ptrdiff_t i, j;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    {
        buckets[j] = buckets[i];
    }

    buckets[k] = 0; memcpy(&buckets[k + 1], buckets, ((size_t)k - 1) * sizeof(int));
}

static int libsais_initialize_buckets_for_lms_suffixes_radix_sort_8u(const unsigned char * RESTRICT T, int * RESTRICT buckets, int first_lms_suffix)
{
    {
        size_t      s = 0;
        ptrdiff_t   c0 = T[first_lms_suffix];
        ptrdiff_t   c1 = 0;

        for (; --first_lms_suffix >= 0; )
        {
            c1 = c0; c0 = T[first_lms_suffix]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1)));
            buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]--;
        }

        buckets[BUCKETS_INDEX4((size_t)c0, (s << 1) & 3)]--;
    }

    {
        int * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

        ptrdiff_t i, j; int sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4(UCHAR_MAX, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
        { 
            temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum; sum += buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 3)]; temp_bucket[j] = sum;
        }

        return sum;
    }
}

static void libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_2k(const int * RESTRICT T, int k, int * RESTRICT buckets, int first_lms_suffix)
{
    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 0)]++;
    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 1)]--;

    ptrdiff_t i; int sum0 = 0, sum1 = 0;
    for (i = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX2((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0))
    { 
        sum0 += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        sum1 += buckets[i + BUCKETS_INDEX2(0, 1)];
        
        buckets[i + BUCKETS_INDEX2(0, 0)] = sum0;
        buckets[i + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static int libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(const int * RESTRICT T, int k, int * RESTRICT buckets, int first_lms_suffix)
{
    {
        size_t      s = 0;
        ptrdiff_t   c0 = T[first_lms_suffix];
        ptrdiff_t   c1 = 0;

        for (; --first_lms_suffix >= 0; )
        {
            c1 = c0; c0 = T[first_lms_suffix]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1)));
            buckets[BUCKETS_INDEX4((size_t)c1, s & 3)]--;
        }

        buckets[BUCKETS_INDEX4((size_t)c0, (s << 1) & 3)]--;
    }

    {
        int * RESTRICT temp_bucket = &buckets[4 * k];

        ptrdiff_t i, j; int sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
        { 
            temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum; sum += buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 3)]; temp_bucket[j] = sum;
        }

        return sum;
    }
}

static void libsais_initialize_buckets_for_radix_and_partial_sorting_32s_4k(const int * RESTRICT T, int k, int * RESTRICT buckets, int first_lms_suffix)
{
    int * RESTRICT bucket_start = &buckets[2 * k];
    int * RESTRICT bucket_end   = &buckets[3 * k];

    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 0)]++;
    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 1)]--;

    ptrdiff_t i, j; int sum0 = 0, sum1 = 0;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    { 
        bucket_start[j] = sum1;

        sum0 += buckets[i + BUCKETS_INDEX2(0, 1)];
        sum1 += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        buckets[i + BUCKETS_INDEX2(0, 1)] = sum0;

        bucket_end[j] = sum1;
    }
}

static void libsais_radix_sort_lms_suffixes_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int m, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];

    ptrdiff_t i, j;
    for (i = (ptrdiff_t)n - 1, j = (ptrdiff_t)n - (ptrdiff_t)m + prefetch_distance + 3; i > j; i -= 4)
    {
        libsais_prefetch(&SA[i - 2 * prefetch_distance]);

        libsais_prefetch(&T[SA[i - prefetch_distance - 0]]);
        libsais_prefetch(&T[SA[i - prefetch_distance - 1]]);
        libsais_prefetch(&T[SA[i - prefetch_distance - 2]]);
        libsais_prefetch(&T[SA[i - prefetch_distance - 3]]);

        int p0 = SA[i - 0]; SA[--induction_bucket[BUCKETS_INDEX2(T[p0], 0)]] = p0;
        int p1 = SA[i - 1]; SA[--induction_bucket[BUCKETS_INDEX2(T[p1], 0)]] = p1;
        int p2 = SA[i - 2]; SA[--induction_bucket[BUCKETS_INDEX2(T[p2], 0)]] = p2;
        int p3 = SA[i - 3]; SA[--induction_bucket[BUCKETS_INDEX2(T[p3], 0)]] = p3;
    }

    for (j -= prefetch_distance + 3; i > j; i -= 1)
    {
        int p = SA[i]; SA[--induction_bucket[BUCKETS_INDEX2(T[p], 0)]] = p;
    }
}

static void libsais_radix_sort_lms_suffixes_32s_2k(const int * RESTRICT T, int * RESTRICT SA, int n, int m, int * RESTRICT induction_bucket)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i, j;
    for (i = (ptrdiff_t)n - 1, j = (ptrdiff_t)n - (ptrdiff_t)m + 2 * prefetch_distance + 3; i > j; i -= 4)
    {
        libsais_prefetch(&SA[i - 3 * prefetch_distance]);
        
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 0]]);
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 1]]);
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 2]]);
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 3]]);

        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 0]], 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 1]], 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 2]], 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 3]], 0)]);

        int p0 = SA[i - 0]; SA[--induction_bucket[BUCKETS_INDEX2(T[p0], 0)]] = p0;
        int p1 = SA[i - 1]; SA[--induction_bucket[BUCKETS_INDEX2(T[p1], 0)]] = p1;
        int p2 = SA[i - 2]; SA[--induction_bucket[BUCKETS_INDEX2(T[p2], 0)]] = p2;
        int p3 = SA[i - 3]; SA[--induction_bucket[BUCKETS_INDEX2(T[p3], 0)]] = p3;
    }

    for (j -= 2 * prefetch_distance + 3; i > j; i -= 1)
    {
        int p = SA[i]; SA[--induction_bucket[BUCKETS_INDEX2(T[p], 0)]] = p;
    }
}

static int libsais_radix_sort_lms_suffixes_32s_1k(const int * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int             i = n - 2;
    int             m = 0;
    size_t          s = 1;
    ptrdiff_t       c0 = T[n - 1];
    ptrdiff_t       c1 = 0;
    ptrdiff_t       c2 = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetch(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[T[i - prefetch_distance - 0]]);
        libsais_prefetchw(&buckets[T[i - prefetch_distance - 1]]);
        libsais_prefetchw(&buckets[T[i - prefetch_distance - 2]]);
        libsais_prefetchw(&buckets[T[i - prefetch_distance - 3]]);

        c1 = T[i - 0]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c0]] = i + 1; m++; }
        
        c0 = T[i - 1]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c1]] = i - 0; m++; }

        c1 = T[i - 2]; s = (s << 1) + (size_t)(c1 > (c0 - (ptrdiff_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c0]] = i - 1; m++; }

        c0 = T[i - 3]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c1]] = i - 2; m++; }
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; s = (s << 1) + (size_t)(c0 > (c1 - (ptrdiff_t)(s & 1))); 
        if ((s & 3) == 1) { SA[--buckets[c2 = c1]] = i + 1; m++; }
    }

    if (m > 1)
    {
        SA[buckets[c2]] = 0;
    }

    return m;
}

static void libsais_radix_sort_set_markers_32s(int * RESTRICT SA, int k, int * RESTRICT induction_bucket, int marker)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)k - 1 - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetch(&induction_bucket[BUCKETS_INDEX2(i + 2 * prefetch_distance, 0)]);

        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 0, 0)]]);
        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 1, 0)]]);
        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 2, 0)]]);
        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 3, 0)]]);

        SA[induction_bucket[BUCKETS_INDEX2(i + 0, 0)]] |= marker;
        SA[induction_bucket[BUCKETS_INDEX2(i + 1, 0)]] |= marker;
        SA[induction_bucket[BUCKETS_INDEX2(i + 2, 0)]] |= marker;
        SA[induction_bucket[BUCKETS_INDEX2(i + 3, 0)]] |= marker;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[induction_bucket[BUCKETS_INDEX2(i, 0)]] |= marker;
    }
}

static void libsais_initialize_buckets_for_partial_sorting_8u(const unsigned char * RESTRICT T, int * RESTRICT buckets, int first_lms_suffix, int left_suffixes_count)
{
    int * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

    buckets[BUCKETS_INDEX4((size_t)T[first_lms_suffix], 1)]++;

    ptrdiff_t i, j; int sum0 = left_suffixes_count + 1, sum1 = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4(UCHAR_MAX, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    { 
        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;

        sum0 += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 2)];
        sum1 += buckets[i + BUCKETS_INDEX4(0, 1)];

        buckets[j + BUCKETS_INDEX2(0, 0)] = sum0;
        buckets[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static void libsais_initialize_buckets_for_partial_sorting_32s_6k(const int * RESTRICT T, int k, int * RESTRICT buckets, int first_lms_suffix, int left_suffixes_count)
{
    int * RESTRICT temp_bucket = &buckets[4 * k];

    buckets[BUCKETS_INDEX4((size_t)T[first_lms_suffix], 1)]++;

    ptrdiff_t i, j; int sum0 = left_suffixes_count + 1, sum1 = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4((ptrdiff_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    { 
        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;

        sum0 += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 2)];
        sum1 += buckets[i + BUCKETS_INDEX4(0, 1)];

        buckets[j + BUCKETS_INDEX2(0, 0)] = sum0;
        buckets[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static int libsais_partial_sorting_scan_left_to_right_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets, int left_suffixes_count, int d)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];
    int * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    SA[induction_bucket[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])]++] = (n - 1) | INT_MIN;
    distinct_names[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])] = ++d;

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)left_suffixes_count - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetch(&SA[i + 2 * prefetch_distance]);

        libsais_prefetch(&T[SA[i + prefetch_distance + 0] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i + prefetch_distance + 0] & INT_MAX] - 2);
        libsais_prefetch(&T[SA[i + prefetch_distance + 1] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i + prefetch_distance + 1] & INT_MAX] - 2);

        int p0 = SA[i + 0]; d += (p0 < 0); p0 &= INT_MAX; int v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] >= T[p0 - 1]);
        SA[induction_bucket[v0]++] = (p0 - 1) | ((distinct_names[v0] != d) << (INT_BIT - 1)); distinct_names[v0] = d;

        int p1 = SA[i + 1]; d += (p1 < 0); p1 &= INT_MAX; int v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] >= T[p1 - 1]);
        SA[induction_bucket[v1]++] = (p1 - 1) | ((distinct_names[v1] != d) << (INT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; d += (p < 0); p &= INT_MAX; int v = BUCKETS_INDEX2(T[p - 1], T[p - 2] >= T[p - 1]);
        SA[induction_bucket[v]++] = (p - 1) | ((distinct_names[v] != d) << (INT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static int libsais_partial_sorting_scan_left_to_right_32s_6k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets, int left_suffixes_count, int d)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[4 * k];
    int * RESTRICT distinct_names   = &buckets[2 * k];

    SA[induction_bucket[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])]++] = (n - 1) | INT_MIN;
    distinct_names[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])] = ++d;

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)left_suffixes_count - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetch(&SA[i + 3 * prefetch_distance]);

        libsais_prefetch(&T[SA[i + 2 * prefetch_distance + 0] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i + 2 * prefetch_distance + 0] & INT_MAX] - 2);
        libsais_prefetch(&T[SA[i + 2 * prefetch_distance + 1] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i + 2 * prefetch_distance + 1] & INT_MAX] - 2);

        int p0 = SA[i + prefetch_distance + 0] & INT_MAX; int v0 = BUCKETS_INDEX2(T[p0 - (p0 > 0)], 0);
        libsais_prefetchw(&induction_bucket[v0]); libsais_prefetchw(&distinct_names[v0]);

        int p1 = SA[i + prefetch_distance + 1] & INT_MAX; int v1 = BUCKETS_INDEX2(T[p1 - (p1 > 0)], 0);
        libsais_prefetchw(&induction_bucket[v1]); libsais_prefetchw(&distinct_names[v1]);

        int p2 = SA[i + 0]; d += (p2 < 0); p2 &= INT_MAX; int v2 = BUCKETS_INDEX2(T[p2 - 1], T[p2 - 2] >= T[p2 - 1]);
        SA[induction_bucket[v2]++] = (p2 - 1) | ((distinct_names[v2] != d) << (INT_BIT - 1)); distinct_names[v2] = d;

        int p3 = SA[i + 1]; d += (p3 < 0); p3 &= INT_MAX; int v3 = BUCKETS_INDEX2(T[p3 - 1], T[p3 - 2] >= T[p3 - 1]);
        SA[induction_bucket[v3]++] = (p3 - 1) | ((distinct_names[v3] != d) << (INT_BIT - 1)); distinct_names[v3] = d;
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; d += (p < 0); p &= INT_MAX; int v = BUCKETS_INDEX2(T[p - 1], T[p - 2] >= T[p - 1]);
        SA[induction_bucket[v]++] = (p - 1) | ((distinct_names[v] != d) << (INT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static int libsais_partial_sorting_scan_left_to_right_32s_4k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets, int d)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[2 * k];
    int * RESTRICT distinct_names   = &buckets[0 * k];

    SA[induction_bucket[T[n - 1]]++] = (n - 1) | ((T[n - 2] < T[n - 1]) << (INT_BIT - 1)) | SUFFIX_GROUP_MARKER;
    distinct_names[BUCKETS_INDEX2(T[n - 1], T[n - 2] < T[n - 1])] = ++d;

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 3 * prefetch_distance]);

        int s0 = SA[i + 2 * prefetch_distance + 0]; const int * Ts0 = &T[s0 & ~SUFFIX_GROUP_MARKER] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL); Ts0--; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i + 2 * prefetch_distance + 1]; const int * Ts1 = &T[s1 & ~SUFFIX_GROUP_MARKER] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL); Ts1--; libsais_prefetch(s1 > 0 ? Ts1 : NULL);
        int s2 = SA[i + 1 * prefetch_distance + 0]; if (s2 > 0) { const ptrdiff_t Ts2 = T[(s2 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts2]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts2, 0)]); }
        int s3 = SA[i + 1 * prefetch_distance + 1]; if (s3 > 0) { const ptrdiff_t Ts3 = T[(s3 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts3]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts3, 0)]); }

        int p0 = SA[i + 0]; SA[i + 0] = p0 & INT_MAX;
        if (p0 > 0)
        {
            SA[i + 0] = 0; d += (p0 >> (SUFFIX_GROUP_BIT - 1)); p0 &= ~SUFFIX_GROUP_MARKER; int v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] < T[p0 - 1]);
            SA[induction_bucket[T[p0 - 1]]++] = (p0 - 1) | ((T[p0 - 2] < T[p0 - 1]) << (INT_BIT - 1)) | ((distinct_names[v0] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v0] = d;
        }

        int p1 = SA[i + 1]; SA[i + 1] = p1 & INT_MAX;
        if (p1 > 0)
        {
            SA[i + 1] = 0; d += (p1 >> (SUFFIX_GROUP_BIT - 1)); p1 &= ~SUFFIX_GROUP_MARKER; int v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] < T[p1 - 1]);
            SA[induction_bucket[T[p1 - 1]]++] = (p1 - 1) | ((T[p1 - 2] < T[p1 - 1]) << (INT_BIT - 1)) | ((distinct_names[v1] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v1] = d;
        }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; SA[i] = p & INT_MAX;
        if (p > 0)
        {
            SA[i] = 0; d += (p >> (SUFFIX_GROUP_BIT - 1)); p &= ~SUFFIX_GROUP_MARKER; int v = BUCKETS_INDEX2(T[p - 1], T[p - 2] < T[p - 1]);
            SA[induction_bucket[T[p - 1]]++] = (p - 1) | ((T[p - 2] < T[p - 1]) << (INT_BIT - 1)) | ((distinct_names[v] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v] = d;
        }
    }

    return d;
}

static void libsais_partial_sorting_scan_left_to_right_32s_1k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT induction_bucket)
{
    const ptrdiff_t prefetch_distance = 32;

    SA[induction_bucket[T[n - 1]]++] = (n - 1) | ((T[n - 2] < T[n - 1]) << (INT_BIT - 1));

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 3 * prefetch_distance]);

        int s0 = SA[i + 2 * prefetch_distance + 0]; const int * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i + 2 * prefetch_distance + 1]; const int * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL);
        int s2 = SA[i + 1 * prefetch_distance + 0]; if (s2 > 0) { libsais_prefetchw(&induction_bucket[T[s2 - 1]]); libsais_prefetch(&T[s2] - 2); }
        int s3 = SA[i + 1 * prefetch_distance + 1]; if (s3 > 0) { libsais_prefetchw(&induction_bucket[T[s3 - 1]]); libsais_prefetch(&T[s3] - 2); }

        int p0 = SA[i + 0]; SA[i + 0] = p0 & INT_MAX; if (p0 > 0) { SA[i + 0] = 0; SA[induction_bucket[T[p0 - 1]]++] = (p0 - 1) | ((T[p0 - 2] < T[p0 - 1]) << (INT_BIT - 1)); }
        int p1 = SA[i + 1]; SA[i + 1] = p1 & INT_MAX; if (p1 > 0) { SA[i + 1] = 0; SA[induction_bucket[T[p1 - 1]]++] = (p1 - 1) | ((T[p1 - 2] < T[p1 - 1]) << (INT_BIT - 1)); }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; SA[i] = p & INT_MAX; if (p > 0) { SA[i] = 0; SA[induction_bucket[T[p - 1]]++] = (p - 1) | ((T[p - 2] < T[p - 1]) << (INT_BIT - 1)); }
    }
}

static void libsais_partial_sorting_shift_markers_8u(int * RESTRICT SA, const int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    const int * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];
    
    ptrdiff_t c;
    for (c = BUCKETS_INDEX2(UCHAR_MAX, 0); c >= BUCKETS_INDEX2(1, 0); c -= BUCKETS_INDEX2(1, 0))
    {
        ptrdiff_t i, j; int s = INT_MIN;
        for (i = (ptrdiff_t)temp_bucket[c] - 1, j = (ptrdiff_t)buckets[c - BUCKETS_INDEX2(1, 0)] + 3; i >= j; i -= 4)
        {
            libsais_prefetchw(&SA[i - prefetch_distance]);

            int p0 = SA[i - 0], q0 = (p0 & INT_MIN) ^ s; s = s ^ q0; SA[i - 0] = p0 ^ q0;
            int p1 = SA[i - 1], q1 = (p1 & INT_MIN) ^ s; s = s ^ q1; SA[i - 1] = p1 ^ q1;
            int p2 = SA[i - 2], q2 = (p2 & INT_MIN) ^ s; s = s ^ q2; SA[i - 2] = p2 ^ q2;
            int p3 = SA[i - 3], q3 = (p3 & INT_MIN) ^ s; s = s ^ q3; SA[i - 3] = p3 ^ q3;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            int p = SA[i], q = (p & INT_MIN) ^ s; s = s ^ q; SA[i] = p ^ q;
        }
    }
}

static void libsais_partial_sorting_shift_markers_32s_6k(int * RESTRICT SA, int k, const int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    const int * RESTRICT temp_bucket = &buckets[4 * k];
    
    ptrdiff_t c;
    for (c = BUCKETS_INDEX2((ptrdiff_t)k - 1, 0); c >= BUCKETS_INDEX2(1, 0); c -= BUCKETS_INDEX2(1, 0))
    {
        ptrdiff_t i, j; int s = INT_MIN;
        for (i = (ptrdiff_t)temp_bucket[c] - 1, j = (ptrdiff_t)buckets[c - BUCKETS_INDEX2(1, 0)] + 3; i >= j; i -= 4)
        {
            libsais_prefetchw(&SA[i - prefetch_distance]);

            int p0 = SA[i - 0], q0 = (p0 & INT_MIN) ^ s; s = s ^ q0; SA[i - 0] = p0 ^ q0;
            int p1 = SA[i - 1], q1 = (p1 & INT_MIN) ^ s; s = s ^ q1; SA[i - 1] = p1 ^ q1;
            int p2 = SA[i - 2], q2 = (p2 & INT_MIN) ^ s; s = s ^ q2; SA[i - 2] = p2 ^ q2;
            int p3 = SA[i - 3], q3 = (p3 & INT_MIN) ^ s; s = s ^ q3; SA[i - 3] = p3 ^ q3;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            int p = SA[i], q = (p & INT_MIN) ^ s; s = s ^ q; SA[i] = p ^ q;
        }
    }
}

static void libsais_partial_sorting_shift_markers_32s_4k(int * RESTRICT SA, int n)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i; int s = SUFFIX_GROUP_MARKER;
    for (i = (ptrdiff_t)n - 1; i >= 3; i -= 4)
    {
        libsais_prefetchw(&SA[i - prefetch_distance]);

        int p0 = SA[i - 0], q0 = ((p0 & SUFFIX_GROUP_MARKER) ^ s) & ((p0 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q0; SA[i - 0] = p0 ^ q0;
        int p1 = SA[i - 1], q1 = ((p1 & SUFFIX_GROUP_MARKER) ^ s) & ((p1 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q1; SA[i - 1] = p1 ^ q1;
        int p2 = SA[i - 2], q2 = ((p2 & SUFFIX_GROUP_MARKER) ^ s) & ((p2 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q2; SA[i - 2] = p2 ^ q2;
        int p3 = SA[i - 3], q3 = ((p3 & SUFFIX_GROUP_MARKER) ^ s) & ((p3 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q3; SA[i - 3] = p3 ^ q3;
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i], q = ((p & SUFFIX_GROUP_MARKER) ^ s) & ((p > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q; SA[i] = p ^ q;
    }
}

static int libsais_partial_sorting_scan_right_to_left_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets, int first_lms_suffix, int left_suffixes_count, int d)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    int * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    ptrdiff_t i, j;
    for (i = (ptrdiff_t)n - (ptrdiff_t)first_lms_suffix - 1, j = (ptrdiff_t)left_suffixes_count + 1 + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetch(&SA[i - 2 * prefetch_distance]);

        libsais_prefetch(&T[SA[i - prefetch_distance - 0] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i - prefetch_distance - 0] & INT_MAX] - 2);
        libsais_prefetch(&T[SA[i - prefetch_distance - 1] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i - prefetch_distance - 1] & INT_MAX] - 2);

        int p0 = SA[i - 0]; d += (p0 < 0); p0 &= INT_MAX; int v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]);
        SA[--induction_bucket[v0]] = (p0 - 1) | ((distinct_names[v0] != d) << (INT_BIT - 1)); distinct_names[v0] = d;

        int p1 = SA[i - 1]; d += (p1 < 0); p1 &= INT_MAX; int v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]);
        SA[--induction_bucket[v1]] = (p1 - 1) | ((distinct_names[v1] != d) << (INT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        int p = SA[i]; d += (p < 0); p &= INT_MAX; int v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
        SA[--induction_bucket[v]] = (p - 1) | ((distinct_names[v] != d) << (INT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static int libsais_partial_sorting_scan_right_to_left_32s_6k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets, int first_lms_suffix, int left_suffixes_count, int d)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[0 * k];
    int * RESTRICT distinct_names   = &buckets[2 * k];

    ptrdiff_t i, j;
    for (i = (ptrdiff_t)n - (ptrdiff_t)first_lms_suffix - 1, j = (ptrdiff_t)left_suffixes_count + 1 + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetch(&SA[i - 3 * prefetch_distance]);

        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 0] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 0] & INT_MAX] - 2);
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 1] & INT_MAX] - 1);
        libsais_prefetch(&T[SA[i - 2 * prefetch_distance - 1] & INT_MAX] - 2);

        int p0 = SA[i - prefetch_distance - 0] & INT_MAX; int v0 = BUCKETS_INDEX2(T[p0 - (p0 > 0)], 0);
        libsais_prefetchw(&induction_bucket[v0]); libsais_prefetchw(&distinct_names[v0]);

        int p1 = SA[i - prefetch_distance - 1] & INT_MAX; int v1 = BUCKETS_INDEX2(T[p1 - (p1 > 0)], 0);
        libsais_prefetchw(&induction_bucket[v1]); libsais_prefetchw(&distinct_names[v1]);

        int p2 = SA[i - 0]; d += (p2 < 0); p2 &= INT_MAX; int v2 = BUCKETS_INDEX2(T[p2 - 1], T[p2 - 2] > T[p2 - 1]);
        SA[--induction_bucket[v2]] = (p2 - 1) | ((distinct_names[v2] != d) << (INT_BIT - 1)); distinct_names[v2] = d;

        int p3 = SA[i - 1]; d += (p3 < 0); p3 &= INT_MAX; int v3 = BUCKETS_INDEX2(T[p3 - 1], T[p3 - 2] > T[p3 - 1]);
        SA[--induction_bucket[v3]] = (p3 - 1) | ((distinct_names[v3] != d) << (INT_BIT - 1)); distinct_names[v3] = d;
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        int p = SA[i]; d += (p < 0); p &= INT_MAX; int v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
        SA[--induction_bucket[v]] = (p - 1) | ((distinct_names[v] != d) << (INT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static int libsais_partial_sorting_scan_right_to_left_32s_4k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets, int d)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[3 * k];
    int * RESTRICT distinct_names   = &buckets[0 * k];

    ptrdiff_t i;
    for (i = (ptrdiff_t)n - 1; i >= 2 * prefetch_distance + 1; i -= 2)
    {
        libsais_prefetchw(&SA[i - 3 * prefetch_distance]);

        int s0 = SA[i - 2 * prefetch_distance - 0]; const int * Ts0 = &T[s0 & ~SUFFIX_GROUP_MARKER] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL); Ts0--; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i - 2 * prefetch_distance - 1]; const int * Ts1 = &T[s1 & ~SUFFIX_GROUP_MARKER] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL); Ts1--; libsais_prefetch(s1 > 0 ? Ts1 : NULL);
        int s2 = SA[i - 1 * prefetch_distance - 0]; if (s2 > 0) { const ptrdiff_t Ts2 = T[(s2 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts2]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts2, 0)]); }
        int s3 = SA[i - 1 * prefetch_distance - 1]; if (s3 > 0) { const ptrdiff_t Ts3 = T[(s3 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts3]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts3, 0)]); }

        int p0 = SA[i - 0];
        if (p0 > 0)
        {
            SA[i - 0] = 0; d += (p0 >> (SUFFIX_GROUP_BIT - 1)); p0 &= ~SUFFIX_GROUP_MARKER; int v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]);
            SA[--induction_bucket[T[p0 - 1]]] = (p0 - 1) | ((T[p0 - 2] > T[p0 - 1]) << (INT_BIT - 1)) | ((distinct_names[v0] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v0] = d;
        }

        int p1 = SA[i - 1];
        if (p1 > 0)
        {
            SA[i - 1] = 0; d += (p1 >> (SUFFIX_GROUP_BIT - 1)); p1 &= ~SUFFIX_GROUP_MARKER; int v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]);
            SA[--induction_bucket[T[p1 - 1]]] = (p1 - 1) | ((T[p1 - 2] > T[p1 - 1]) << (INT_BIT - 1)) | ((distinct_names[v1] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v1] = d;
        }
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i];
        if (p > 0)
        {
            SA[i] = 0; d += (p >> (SUFFIX_GROUP_BIT - 1)); p &= ~SUFFIX_GROUP_MARKER; int v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
            SA[--induction_bucket[T[p - 1]]] = (p - 1) | ((T[p - 2] > T[p - 1]) << (INT_BIT - 1)) | ((distinct_names[v] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v] = d;
        }
    }

    return d;
}

static void libsais_partial_sorting_scan_right_to_left_32s_1k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT induction_bucket)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i;
    for (i = (ptrdiff_t)n - 1; i >= 2 * prefetch_distance + 1; i -= 2)
    {
        libsais_prefetchw(&SA[i - 3 * prefetch_distance]);

        int s0 = SA[i - 2 * prefetch_distance - 0]; const int * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i - 2 * prefetch_distance - 1]; const int * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL);
        int s2 = SA[i - 1 * prefetch_distance - 0]; if (s2 > 0) { libsais_prefetchw(&induction_bucket[T[s2 - 1]]); libsais_prefetch(&T[s2] - 2); }
        int s3 = SA[i - 1 * prefetch_distance - 1]; if (s3 > 0) { libsais_prefetchw(&induction_bucket[T[s3 - 1]]); libsais_prefetch(&T[s3] - 2); }

        int p0 = SA[i - 0]; if (p0 > 0) { SA[i - 0] = 0; SA[--induction_bucket[T[p0 - 1]]] = (p0 - 1) | ((T[p0 - 2] > T[p0 - 1]) << (INT_BIT - 1)); }
        int p1 = SA[i - 1]; if (p1 > 0) { SA[i - 1] = 0; SA[--induction_bucket[T[p1 - 1]]] = (p1 - 1) | ((T[p1 - 2] > T[p1 - 1]) << (INT_BIT - 1)); }
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i]; if (p > 0) { SA[i] = 0; SA[--induction_bucket[T[p - 1]]] = (p - 1) | ((T[p - 2] > T[p - 1]) << (INT_BIT - 1)); }
    }
}

static void libsais_partial_sorting_gather_lms_suffixes_32s_4k(int * RESTRICT SA, int n)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i, j, l;
    for (i = 0, j = (ptrdiff_t)n - 3, l = 0; i < j; i += 4)
    {
        libsais_prefetch(&SA[i + prefetch_distance]);

        int s0 = SA[i + 0]; SA[l] = (s0 - SUFFIX_GROUP_MARKER) & (~SUFFIX_GROUP_MARKER); l += (s0 < 0);
        int s1 = SA[i + 1]; SA[l] = (s1 - SUFFIX_GROUP_MARKER) & (~SUFFIX_GROUP_MARKER); l += (s1 < 0);
        int s2 = SA[i + 2]; SA[l] = (s2 - SUFFIX_GROUP_MARKER) & (~SUFFIX_GROUP_MARKER); l += (s2 < 0);
        int s3 = SA[i + 3]; SA[l] = (s3 - SUFFIX_GROUP_MARKER) & (~SUFFIX_GROUP_MARKER); l += (s3 < 0);
    }

    for (j += 3; i < j; i += 1)
    {
        int s = SA[i]; SA[l] = (s - SUFFIX_GROUP_MARKER) & (~SUFFIX_GROUP_MARKER); l += (s < 0);
    }
}

static void libsais_partial_sorting_gather_lms_suffixes_32s_1k(int * RESTRICT SA, int n)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i, j, l;
    for (i = 0, j = (ptrdiff_t)n - 3, l = 0; i < j; i += 4)
    {
        libsais_prefetch(&SA[i + prefetch_distance]);

        int s0 = SA[i + 0]; SA[l] = s0 & INT_MAX; l += (s0 < 0);
        int s1 = SA[i + 1]; SA[l] = s1 & INT_MAX; l += (s1 < 0);
        int s2 = SA[i + 2]; SA[l] = s2 & INT_MAX; l += (s2 < 0);
        int s3 = SA[i + 3]; SA[l] = s3 & INT_MAX; l += (s3 < 0);
    }

    for (j += 3; i < j; i += 1)
    {
        int s = SA[i]; SA[l] = s & INT_MAX; l += (s < 0);
    }
}

static void libsais_induce_partial_order_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets, int first_lms_suffix, int left_suffixes_count)
{
    memset(&buckets[2 * ALPHABET_SIZE], 0, 2 * ALPHABET_SIZE * sizeof(int));

    int d = libsais_partial_sorting_scan_left_to_right_8u(T, SA, n, buckets, left_suffixes_count, 0);
    libsais_partial_sorting_shift_markers_8u(SA, buckets);
    libsais_partial_sorting_scan_right_to_left_8u(T, SA, n, buckets, first_lms_suffix, left_suffixes_count, d);
}

static void libsais_induce_partial_order_32s_6k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets, int first_lms_suffix, int left_suffixes_count)
{
    memset(&buckets[2 * k], 0, 2 * (size_t)k * sizeof(int));

    int d = libsais_partial_sorting_scan_left_to_right_32s_6k(T, SA, n, k, buckets, left_suffixes_count, 0);
    libsais_partial_sorting_shift_markers_32s_6k(SA, k, buckets);
    libsais_partial_sorting_scan_right_to_left_32s_6k(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count, d);
}

static void libsais_induce_partial_order_32s_4k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(int));

    int d = libsais_partial_sorting_scan_left_to_right_32s_4k(T, SA, n, k, buckets, 0);
    libsais_partial_sorting_shift_markers_32s_4k(SA, n);
    libsais_partial_sorting_scan_right_to_left_32s_4k(T, SA, n, k, buckets, d);
    libsais_partial_sorting_gather_lms_suffixes_32s_4k(SA, n);
}

static void libsais_induce_partial_order_32s_2k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    libsais_partial_sorting_scan_left_to_right_32s_1k(T, SA, n, k, &buckets[1 * k]);
    libsais_partial_sorting_scan_right_to_left_32s_1k(T, SA, n, k, &buckets[0 * k]);
    libsais_partial_sorting_gather_lms_suffixes_32s_1k(SA, n);
}

static void libsais_induce_partial_order_32s_1k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_start_32s_1k(k, buckets);
    libsais_partial_sorting_scan_left_to_right_32s_1k(T, SA, n, k, buckets);

    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_end_32s_1k(k, buckets);
    libsais_partial_sorting_scan_right_to_left_32s_1k(T, SA, n, k, buckets);

    libsais_partial_sorting_gather_lms_suffixes_32s_1k(SA, n);
}

static int libsais_renumber_and_gather_lms_suffixes_8u(int * RESTRICT SA, int n, int m)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT SAm = &SA[m];

    memset(SAm, 0, ((size_t)n >> 1) * sizeof(int));

    ptrdiff_t i, j; int name = 0;
    for (i = 0, j = (ptrdiff_t)m - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetch(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 0] & INT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 1] & INT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 2] & INT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 3] & INT_MAX) >> 1]);

        int p0 = SA[i + 0]; SAm[(p0 & INT_MAX) >> 1] = name | INT_MIN; name += p0 < 0;
        int p1 = SA[i + 1]; SAm[(p1 & INT_MAX) >> 1] = name | INT_MIN; name += p1 < 0;
        int p2 = SA[i + 2]; SAm[(p2 & INT_MAX) >> 1] = name | INT_MIN; name += p2 < 0;
        int p3 = SA[i + 3]; SAm[(p3 & INT_MAX) >> 1] = name | INT_MIN; name += p3 < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        int p = SA[i]; SAm[(p & INT_MAX) >> 1] = name | INT_MIN; name += p < 0;
    }

    if (name < m)
    {
        ptrdiff_t l;
        for (i = (ptrdiff_t)m + ((ptrdiff_t)n >> 1) - 1, j = (ptrdiff_t)m + 3, l = (ptrdiff_t)n - 1; i >= j; i -= 4)
        {
            libsais_prefetch(&SA[i - prefetch_distance]);

            int s0 = SA[i - 0]; SA[l] = s0 & INT_MAX; l -= s0 < 0;
            int s1 = SA[i - 1]; SA[l] = s1 & INT_MAX; l -= s1 < 0;
            int s2 = SA[i - 2]; SA[l] = s2 & INT_MAX; l -= s2 < 0;
            int s3 = SA[i - 3]; SA[l] = s3 & INT_MAX; l -= s3 < 0;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            int s = SA[i]; SA[l] = s & INT_MAX; l -= s < 0;
        }
    }
    else
    {
        for (i = 0; i < m; i += 1) { SA[i] &= INT_MAX; }
    }

    return name;
}

static int libsais_renumber_and_mark_distinct_lms_suffixes_32s_4k(int * RESTRICT SA, int n, int m)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT SAm = &SA[m];

    memset(SAm, 0, ((size_t)n >> 1) * sizeof(int));

    ptrdiff_t i, j; int p0, p1, p2, p3 = -1, name = 1;
    for (i = 0, j = (ptrdiff_t)m - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 0] & INT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 1] & INT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 2] & INT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 3] & INT_MAX) >> 1]);

        p0 = SA[i + 0]; SAm[(SA[i + 0] = p0 & INT_MAX) >> 1] = name | (p0 & p3 & INT_MIN); name += p0 < 0;
        p1 = SA[i + 1]; SAm[(SA[i + 1] = p1 & INT_MAX) >> 1] = name | (p1 & p0 & INT_MIN); name += p1 < 0;
        p2 = SA[i + 2]; SAm[(SA[i + 2] = p2 & INT_MAX) >> 1] = name | (p2 & p1 & INT_MIN); name += p2 < 0;
        p3 = SA[i + 3]; SAm[(SA[i + 3] = p3 & INT_MAX) >> 1] = name | (p3 & p2 & INT_MIN); name += p3 < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        p2 = p3; p3 = SA[i]; SAm[(SA[i] = p3 & INT_MAX) >> 1] = name | (p3 & p2 & INT_MIN); name += p3 < 0;
    }
    
    if (name <= m)
    {
        p3 = -1;
        for (i = m, j = (ptrdiff_t)m + ((ptrdiff_t)n >> 1) - 3; i < j; i += 4)
        {
            libsais_prefetchw(&SA[i + prefetch_distance]);

            p0 = SA[i + 0]; SA[i + 0] = p0 & (p3 | INT_MAX); p0 = (p0 == 0) ? p3 : p0;
            p1 = SA[i + 1]; SA[i + 1] = p1 & (p0 | INT_MAX); p1 = (p1 == 0) ? p0 : p1;
            p2 = SA[i + 2]; SA[i + 2] = p2 & (p1 | INT_MAX); p2 = (p2 == 0) ? p1 : p2;
            p3 = SA[i + 3]; SA[i + 3] = p3 & (p2 | INT_MAX); p3 = (p3 == 0) ? p2 : p3;
        }

        for (j += 3; i < j; i += 1)
        {
            p2 = p3; p3 = SA[i]; SA[i] = p3 & (p2 | INT_MAX); p3 = (p3 == 0) ? p2 : p3;
        }
    }

    return name - 1;
}

static int libsais_renumber_and_mark_distinct_lms_suffixes_32s_1k(int * RESTRICT T, int * RESTRICT SA, int n, int m)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT SAm = &SA[m];

    {
        libsais_gather_lms_suffixes_32s(T, SA, n);

        memset(&SA[m], 0, ((size_t)n - (size_t)m - (size_t)m) * sizeof(int));

        ptrdiff_t i, j;
        for (i = (ptrdiff_t)n - (ptrdiff_t)m, j = (ptrdiff_t)n - 1 - prefetch_distance - 3; i < j; i += 4)
        {
            libsais_prefetch(&SA[i + 2 * prefetch_distance]);

            libsais_prefetchw(&SAm[((unsigned int)SA[i + prefetch_distance + 0]) >> 1]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + prefetch_distance + 1]) >> 1]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + prefetch_distance + 2]) >> 1]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + prefetch_distance + 3]) >> 1]);

            SAm[((unsigned int)SA[i + 0]) >> 1] = SA[i + 1] - SA[i + 0] + 1 + INT_MIN;
            SAm[((unsigned int)SA[i + 1]) >> 1] = SA[i + 2] - SA[i + 1] + 1 + INT_MIN;
            SAm[((unsigned int)SA[i + 2]) >> 1] = SA[i + 3] - SA[i + 2] + 1 + INT_MIN;
            SAm[((unsigned int)SA[i + 3]) >> 1] = SA[i + 4] - SA[i + 3] + 1 + INT_MIN;
        }

        for (j += prefetch_distance + 3; i < j; i += 1)
        {
            SAm[((unsigned int)SA[i]) >> 1] = SA[i + 1] - SA[i] + 1 + INT_MIN;
        }

        SAm[((unsigned int)SA[n - 1]) >> 1] = 1 + INT_MIN;
    }

    {
        ptrdiff_t i, j;
        for (i = 0, j = (ptrdiff_t)(n >> 1) - 3; i < j; i += 4)
        {
            libsais_prefetchw(&SAm[i + prefetch_distance]);

            SAm[i + 0] = (SAm[i + 0] < 0 ? SAm[i + 0] : 0) & INT_MAX;
            SAm[i + 1] = (SAm[i + 1] < 0 ? SAm[i + 1] : 0) & INT_MAX;
            SAm[i + 2] = (SAm[i + 2] < 0 ? SAm[i + 2] : 0) & INT_MAX;
            SAm[i + 3] = (SAm[i + 3] < 0 ? SAm[i + 3] : 0) & INT_MAX;
        }

        for (j += 3; i < j; i += 1)
        {
            SAm[i] = (SAm[i] < 0 ? SAm[i] : 0) & INT_MAX;
        }
    }

    int name = 1;

    {
        ptrdiff_t i, j, p = SA[0], plen = SAm[p >> 1]; int pdiff = INT_MIN;
        for (i = 1, j = m - prefetch_distance - 1; i < j; i += 2)
        {
            libsais_prefetch(&SA[i + 2 * prefetch_distance]);
            
            libsais_prefetchw(&SAm[((unsigned int)SA[i + prefetch_distance + 0]) >> 1]); libsais_prefetch(&T[((unsigned int)SA[i + prefetch_distance + 0])]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + prefetch_distance + 1]) >> 1]); libsais_prefetch(&T[((unsigned int)SA[i + prefetch_distance + 1])]);

            ptrdiff_t q = SA[i + 0], qlen = SAm[q >> 1]; int qdiff = INT_MIN;
            if (plen == qlen) { ptrdiff_t l = 0; do { if (T[p + l] != T[q + l]) { break; } } while (++l < qlen); qdiff = (l - qlen) & INT_MIN; }
            SAm[p >> 1] = name | (pdiff & qdiff); name += (qdiff < 0);

            p = SA[i + 1]; plen = SAm[p >> 1]; pdiff = INT_MIN;
            if (qlen == plen) { ptrdiff_t l = 0; do { if (T[q + l] != T[p + l]) { break; } } while (++l < plen); pdiff = (l - plen) & INT_MIN; }
            SAm[q >> 1] = name | (qdiff & pdiff); name += (pdiff < 0);
        }

        for (j += prefetch_distance + 1; i < j; i += 1)
        {
            ptrdiff_t q = SA[i], qlen = SAm[q >> 1]; int qdiff = INT_MIN;
            if (plen == qlen) { ptrdiff_t l = 0; do { if (T[p + l] != T[q + l]) { break; } } while (++l < plen); qdiff = (l - plen) & INT_MIN; }
            SAm[p >> 1] = name | (pdiff & qdiff); name += (qdiff < 0);

            p = q; plen = qlen; pdiff = qdiff;
        }

        SAm[p >> 1] = name | pdiff; name++;
    }

    if (name <= m)
    {
        ptrdiff_t i, j; int p0, p1, p2, p3 = -1;
        for (i = m, j = (ptrdiff_t)m + ((ptrdiff_t)n >> 1) - 3; i < j; i += 4)
        {
            libsais_prefetchw(&SA[i + prefetch_distance]);

            p0 = SA[i + 0]; SA[i + 0] = p0 & (p3 | INT_MAX); p0 = (p0 == 0) ? p3 : p0;
            p1 = SA[i + 1]; SA[i + 1] = p1 & (p0 | INT_MAX); p1 = (p1 == 0) ? p0 : p1;
            p2 = SA[i + 2]; SA[i + 2] = p2 & (p1 | INT_MAX); p2 = (p2 == 0) ? p1 : p2;
            p3 = SA[i + 3]; SA[i + 3] = p3 & (p2 | INT_MAX); p3 = (p3 == 0) ? p2 : p3;
        }

        for (j += 3; i < j; i += 1)
        {
            p2 = p3; p3 = SA[i]; SA[i] = p3 & (p2 | INT_MAX); p3 = (p3 == 0) ? p2 : p3;
        }
    }

    return name - 1;
}

static void libsais_reconstruct_lms_suffixes(int * RESTRICT SA, int n, int m)
{
    const ptrdiff_t prefetch_distance = 32;

    const int * RESTRICT SAnm = &SA[n - m];

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)m - prefetch_distance - 3; i < j; i += 4)
    { 
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        libsais_prefetch(&SAnm[SA[i + prefetch_distance + 0]]); 
        libsais_prefetch(&SAnm[SA[i + prefetch_distance + 1]]);
        libsais_prefetch(&SAnm[SA[i + prefetch_distance + 2]]);
        libsais_prefetch(&SAnm[SA[i + prefetch_distance + 3]]);

        SA[i + 0] = SAnm[SA[i + 0]];
        SA[i + 1] = SAnm[SA[i + 1]];
        SA[i + 2] = SAnm[SA[i + 2]];
        SA[i + 3] = SAnm[SA[i + 3]];
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[i] = SAnm[SA[i]];
    }
}

static void libsais_place_lms_suffixes_interval_8u(int * RESTRICT SA, int n, int m, const int * RESTRICT buckets)
{
    const int * RESTRICT bucket_end = &buckets[7 * ALPHABET_SIZE];

    ptrdiff_t c, j = n;
    for (c = UCHAR_MAX - 1; c >= 0; --c)
    {
        ptrdiff_t l = (ptrdiff_t)buckets[BUCKETS_INDEX2(c, 1) + BUCKETS_INDEX2(1, 0)] - (ptrdiff_t)buckets[BUCKETS_INDEX2(c, 1)];
        if (l > 0)
        {
            ptrdiff_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(int));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (int)l], (size_t)l * sizeof(int));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(int));
}

static void libsais_place_lms_suffixes_interval_32s_4k(int * RESTRICT SA, int n, int k, int m, const int * RESTRICT buckets)
{
    const int * RESTRICT bucket_end = &buckets[3 * k];

    ptrdiff_t c, j = n;
    for (c = (ptrdiff_t)k - 2; c >= 0; --c)
    {
        ptrdiff_t l = (ptrdiff_t)buckets[BUCKETS_INDEX2(c, 1) + BUCKETS_INDEX2(1, 0)] - (ptrdiff_t)buckets[BUCKETS_INDEX2(c, 1)];
        if (l > 0)
        {
            ptrdiff_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(int));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (int)l], (size_t)l * sizeof(int));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(int));
}

static void libsais_place_lms_suffixes_interval_32s_2k(int * RESTRICT SA, int n, int k, int m, const int * RESTRICT buckets)
{
    ptrdiff_t c, j = n;
    for (c = BUCKETS_INDEX2((ptrdiff_t)k - 2, 0); c >= BUCKETS_INDEX2(0, 0); c -= BUCKETS_INDEX2(1, 0))
    {
        ptrdiff_t l = (ptrdiff_t)buckets[c + BUCKETS_INDEX2(1, 1)] - (ptrdiff_t)buckets[c + BUCKETS_INDEX2(0, 1)];
        if (l > 0)
        {
            ptrdiff_t i = buckets[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(int));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (int)l], (size_t)l * sizeof(int));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(int));
}

static void libsais_place_lms_suffixes_interval_32s_1k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int m, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int c = k - 1; ptrdiff_t i, l = buckets[c];
    for (i = (ptrdiff_t)m - 1; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetch(&SA[i - 2 * prefetch_distance]);

        libsais_prefetch(&T[SA[i - prefetch_distance - 0]]);
        libsais_prefetch(&T[SA[i - prefetch_distance - 1]]);
        libsais_prefetch(&T[SA[i - prefetch_distance - 2]]);
        libsais_prefetch(&T[SA[i - prefetch_distance - 3]]);

        int p0 = SA[i - 0]; if (T[p0] != c) { c = T[p0]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(int)); l = buckets[c]; } SA[--l] = p0;
        int p1 = SA[i - 1]; if (T[p1] != c) { c = T[p1]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(int)); l = buckets[c]; } SA[--l] = p1;
        int p2 = SA[i - 2]; if (T[p2] != c) { c = T[p2]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(int)); l = buckets[c]; } SA[--l] = p2;
        int p3 = SA[i - 3]; if (T[p3] != c) { c = T[p3]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(int)); l = buckets[c]; } SA[--l] = p3;
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i]; if (T[p] != c) { c = T[p]; memset(&SA[buckets[c]], 0, (size_t)(l - buckets[c]) * sizeof(int)); l = buckets[c]; } SA[--l] = p;
    }

    memset(&SA[0], 0, (size_t)l * sizeof(int));
}

static void libsais_place_lms_suffixes_histogram_32s_6k(int * RESTRICT SA, int n, int k, int m, const int * RESTRICT buckets)
{
    const int * RESTRICT bucket_end = &buckets[5 * k];

    ptrdiff_t c, j = n;
    for (c = (ptrdiff_t)k - 2; c >= 0; --c)
    {
        ptrdiff_t l = (ptrdiff_t)buckets[BUCKETS_INDEX4(c, 1)];
        if (l > 0)
        {
            ptrdiff_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(int));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (int)l], (size_t)l * sizeof(int));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(int));
}

static void libsais_place_lms_suffixes_histogram_32s_4k(int * RESTRICT SA, int n, int k, int m, const int * RESTRICT buckets)
{
    const int * RESTRICT bucket_end = &buckets[3 * k];

    ptrdiff_t c, j = n;
    for (c = (ptrdiff_t)k - 2; c >= 0; --c)
    {
        ptrdiff_t l = (ptrdiff_t)buckets[BUCKETS_INDEX2(c, 1)];
        if (l > 0)
        {
            ptrdiff_t i = bucket_end[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(int));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (int)l], (size_t)l * sizeof(int));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(int));
}

static void libsais_place_lms_suffixes_histogram_32s_2k(int * RESTRICT SA, int n, int k, int m, const int * RESTRICT buckets)
{
    ptrdiff_t c, j = n;
    for (c = BUCKETS_INDEX2((ptrdiff_t)k - 2, 0); c >= BUCKETS_INDEX2(0, 0); c -= BUCKETS_INDEX2(1, 0))
    {
        ptrdiff_t l = (ptrdiff_t)buckets[c + BUCKETS_INDEX2(0, 1)];
        if (l > 0)
        {
            ptrdiff_t i = buckets[c];
            if (j - i > 0)
            {
                memset(&SA[i], 0, (size_t)(j - i) * sizeof(int));
            }

            memmove(&SA[j = (i - l)], &SA[m -= (int)l], (size_t)l * sizeof(int));
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(int));
}

static void libsais_final_bwt_scan_left_to_right_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[6 * ALPHABET_SIZE];

    SA[induction_bucket[T[n - 1]]++] = (n - 1) | ((T[n - 2] < T[n - 1]) << (INT_BIT - 1));

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        int s0 = SA[i + prefetch_distance + 0]; const unsigned char * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL); Ts0--; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i + prefetch_distance + 1]; const unsigned char * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL); Ts1--; libsais_prefetch(s1 > 0 ? Ts1 : NULL);

        int p0 = SA[i + 0]; SA[i + 0] = p0 & INT_MAX; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | INT_MIN; SA[induction_bucket[T[p0]]++] = p0 | (((T[p0 - (p0 > 0)] < T[p0])) << (INT_BIT - 1)); }
        int p1 = SA[i + 1]; SA[i + 1] = p1 & INT_MAX; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | INT_MIN; SA[induction_bucket[T[p1]]++] = p1 | (((T[p1 - (p1 > 0)] < T[p1])) << (INT_BIT - 1)); }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; SA[i] = p & INT_MAX; if (p > 0) { p--; SA[i] = T[p] | INT_MIN; SA[induction_bucket[T[p]]++] = p | (((T[p - (p > 0)] < T[p])) << (INT_BIT - 1)); }
    }
}

static void libsais_final_sorting_scan_left_to_right_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[6 * ALPHABET_SIZE];

    SA[induction_bucket[T[n - 1]]++] = (n - 1) | ((T[n - 2] < T[n - 1]) << (INT_BIT - 1));

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        int s0 = SA[i + prefetch_distance + 0]; const unsigned char * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL); Ts0--; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i + prefetch_distance + 1]; const unsigned char * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL); Ts1--; libsais_prefetch(s1 > 0 ? Ts1 : NULL);

        int p0 = SA[i + 0]; SA[i + 0] = p0 ^ INT_MIN; if (p0 > 0) { p0--; SA[induction_bucket[T[p0]]++] = p0 | (((T[p0 - (p0 > 0)] < T[p0])) << (INT_BIT - 1)); }
        int p1 = SA[i + 1]; SA[i + 1] = p1 ^ INT_MIN; if (p1 > 0) { p1--; SA[induction_bucket[T[p1]]++] = p1 | (((T[p1 - (p1 > 0)] < T[p1])) << (INT_BIT - 1)); }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; SA[i] = p ^ INT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | (((T[p - (p > 0)] < T[p])) << (INT_BIT - 1)); }
    }
}

static void libsais_final_sorting_scan_left_to_right_32s(const int * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT induction_bucket)
{
    const ptrdiff_t prefetch_distance = 32;

    SA[induction_bucket[T[n - 1]]++] = (n - 1) | ((T[n - 2] < T[n - 1]) << (INT_BIT - 1));

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 3 * prefetch_distance]);

        int s0 = SA[i + 2 * prefetch_distance + 0]; const int * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i + 2 * prefetch_distance + 1]; const int * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL);
        int s2 = SA[i + 1 * prefetch_distance + 0]; if (s2 > 0) { libsais_prefetchw(&induction_bucket[T[s2 - 1]]); libsais_prefetch(&T[s2] - 2); }
        int s3 = SA[i + 1 * prefetch_distance + 1]; if (s3 > 0) { libsais_prefetchw(&induction_bucket[T[s3 - 1]]); libsais_prefetch(&T[s3] - 2); }

        int p0 = SA[i + 0]; SA[i + 0] = p0 ^ INT_MIN; if (p0 > 0) { p0--; SA[induction_bucket[T[p0]]++] = p0 | (((T[p0 - (p0 > 0)] < T[p0])) << (INT_BIT - 1)); }
        int p1 = SA[i + 1]; SA[i + 1] = p1 ^ INT_MIN; if (p1 > 0) { p1--; SA[induction_bucket[T[p1]]++] = p1 | (((T[p1 - (p1 > 0)] < T[p1])) << (INT_BIT - 1)); }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        int p = SA[i]; SA[i] = p ^ INT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | (((T[p - (p > 0)] < T[p])) << (INT_BIT - 1)); }
    }
}

static int libsais_final_bwt_scan_right_to_left_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[7 * ALPHABET_SIZE];

    ptrdiff_t i; int index = -1;
    for (i = (ptrdiff_t)n - 1; i >= prefetch_distance + 1; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        int s0 = SA[i - prefetch_distance - 0]; const unsigned char * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL); Ts0--; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i - prefetch_distance - 1]; const unsigned char * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL); Ts1--; libsais_prefetch(s1 > 0 ? Ts1 : NULL);

        int p0 = SA[i - 0]; index = (p0 == 0) ? (int)(i - 0) : index;
        SA[i - 0] = p0 & INT_MAX; if (p0 > 0) { p0--; unsigned char c0 = T[p0 - (p0 > 0)], c1 = T[p0]; SA[i - 0] = c1; int t = c0 | INT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p0 : t; }

        int p1 = SA[i - 1]; index = (p1 == 0) ? (int)(i - 1) : index;
        SA[i - 1] = p1 & INT_MAX; if (p1 > 0) { p1--; unsigned char c0 = T[p1 - (p1 > 0)], c1 = T[p1]; SA[i - 1] = c1; int t = c0 | INT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p1 : t; }
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i]; index = (p == 0) ? (int)i : index;
        SA[i] = p & INT_MAX; if (p > 0) { p--; unsigned char c0 = T[p - (p > 0)], c1 = T[p]; SA[i] = c1; int t = c0 | INT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; }
    }

    return index;
}

static void libsais_final_sorting_scan_right_to_left_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT buckets)
{
    const ptrdiff_t prefetch_distance = 32;

    int * RESTRICT induction_bucket = &buckets[7 * ALPHABET_SIZE];

    ptrdiff_t i;
    for (i = (ptrdiff_t)n - 1; i >= prefetch_distance + 1; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        int s0 = SA[i - prefetch_distance - 0]; const unsigned char * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL); Ts0--; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i - prefetch_distance - 1]; const unsigned char * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL); Ts1--; libsais_prefetch(s1 > 0 ? Ts1 : NULL);

        int p0 = SA[i - 0]; SA[i - 0] = p0 & INT_MAX; if (p0 > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | (((T[p0 - (p0 > 0)] > T[p0])) << (INT_BIT - 1)); }
        int p1 = SA[i - 1]; SA[i - 1] = p1 & INT_MAX; if (p1 > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | (((T[p1 - (p1 > 0)] > T[p1])) << (INT_BIT - 1)); }
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i]; SA[i] = p & INT_MAX; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | (((T[p - (p > 0)] > T[p])) << (INT_BIT - 1)); }
    }
}

static void libsais_final_sorting_scan_right_to_left_32s(const int * RESTRICT T, int * RESTRICT SA, int n, int * RESTRICT induction_bucket)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i;
    for (i = (ptrdiff_t)n - 1; i >= 2 * prefetch_distance + 1; i -= 2)
    {
        libsais_prefetchw(&SA[i - 3 * prefetch_distance]);

        int s0 = SA[i - 2 * prefetch_distance - 0]; const int * Ts0 = &T[s0] - 1; libsais_prefetch(s0 > 0 ? Ts0 : NULL);
        int s1 = SA[i - 2 * prefetch_distance - 1]; const int * Ts1 = &T[s1] - 1; libsais_prefetch(s1 > 0 ? Ts1 : NULL);
        int s2 = SA[i - 1 * prefetch_distance - 0]; if (s2 > 0) { libsais_prefetchw(&induction_bucket[T[s2 - 1]]); libsais_prefetch(&T[s2] - 2); }
        int s3 = SA[i - 1 * prefetch_distance - 1]; if (s3 > 0) { libsais_prefetchw(&induction_bucket[T[s3 - 1]]); libsais_prefetch(&T[s3] - 2); }

        int p0 = SA[i - 0]; SA[i - 0] = p0 & INT_MAX; if (p0 > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | (((T[p0 - (p0 > 0)] > T[p0])) << (INT_BIT - 1)); }
        int p1 = SA[i - 1]; SA[i - 1] = p1 & INT_MAX; if (p1 > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | (((T[p1 - (p1 > 0)] > T[p1])) << (INT_BIT - 1)); }
    }

    for (; i >= 0; i -= 1)
    {
        int p = SA[i]; SA[i] = p & INT_MAX; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | (((T[p - (p > 0)] > T[p])) << (INT_BIT - 1)); }
    }
}

static int libsais_induce_final_order_8u(const unsigned char * RESTRICT T, int * RESTRICT SA, int n, int bwt, int * RESTRICT buckets)
{
    if (bwt)
    {
        libsais_final_bwt_scan_left_to_right_8u(T, SA, n, buckets);
        return libsais_final_bwt_scan_right_to_left_8u(T, SA, n, buckets);
    }
    else
    {
        libsais_final_sorting_scan_left_to_right_8u(T, SA, n, buckets);
        libsais_final_sorting_scan_right_to_left_8u(T, SA, n, buckets);
        return 0;
    }
}

static void libsais_induce_final_order_32s_6k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    libsais_final_sorting_scan_left_to_right_32s(T, SA, n, &buckets[4 * k]);
    libsais_final_sorting_scan_right_to_left_32s(T, SA, n, &buckets[5 * k]);
}

static void libsais_induce_final_order_32s_4k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    libsais_final_sorting_scan_left_to_right_32s(T, SA, n, &buckets[2 * k]);
    libsais_final_sorting_scan_right_to_left_32s(T, SA, n, &buckets[3 * k]);
}

static void libsais_induce_final_order_32s_2k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    libsais_final_sorting_scan_left_to_right_32s(T, SA, n, &buckets[1 * k]);
    libsais_final_sorting_scan_right_to_left_32s(T, SA, n, &buckets[0 * k]);
}

static void libsais_induce_final_order_32s_1k(const int * RESTRICT T, int * RESTRICT SA, int n, int k, int * RESTRICT buckets)
{
    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_start_32s_1k(k, buckets);
    libsais_final_sorting_scan_left_to_right_32s(T, SA, n, buckets);

    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_end_32s_1k(k, buckets);
    libsais_final_sorting_scan_right_to_left_32s(T, SA, n, buckets);
}

static int libsais_compact_lms_suffixes_32s(int * RESTRICT T, int * RESTRICT SA, int n, int m, int fs)
{
    const ptrdiff_t prefetch_distance = 32;

    int f = 0;

    {
        int * RESTRICT SAm = &SA[m];

        int i, j;
        for (i = 0, j = m - 2 * (int)prefetch_distance - 3; i < j; i += 4)
        {
            libsais_prefetch(&SA[i + 3 * prefetch_distance]);

            libsais_prefetchw(&SAm[((unsigned int)SA[i + 2 * prefetch_distance + 0]) >> 1]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + 2 * prefetch_distance + 1]) >> 1]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + 2 * prefetch_distance + 2]) >> 1]);
            libsais_prefetchw(&SAm[((unsigned int)SA[i + 2 * prefetch_distance + 3]) >> 1]);

            unsigned int q0 = (unsigned int)SA[i + prefetch_distance + 0]; const int * Tq0 = &T[q0]; libsais_prefetchw(SAm[q0 >> 1] < 0 ? Tq0 : NULL);
            unsigned int q1 = (unsigned int)SA[i + prefetch_distance + 1]; const int * Tq1 = &T[q1]; libsais_prefetchw(SAm[q1 >> 1] < 0 ? Tq1 : NULL);
            unsigned int q2 = (unsigned int)SA[i + prefetch_distance + 2]; const int * Tq2 = &T[q2]; libsais_prefetchw(SAm[q2 >> 1] < 0 ? Tq2 : NULL);
            unsigned int q3 = (unsigned int)SA[i + prefetch_distance + 3]; const int * Tq3 = &T[q3]; libsais_prefetchw(SAm[q3 >> 1] < 0 ? Tq3 : NULL);

            unsigned int p0 = (unsigned int)SA[i + 0]; int s0 = SAm[p0 >> 1]; if (s0 < 0) { T[p0] |= INT_MIN; f++; s0 = i + 0 + INT_MIN + f; } SAm[p0 >> 1] = s0 - f;
            unsigned int p1 = (unsigned int)SA[i + 1]; int s1 = SAm[p1 >> 1]; if (s1 < 0) { T[p1] |= INT_MIN; f++; s1 = i + 1 + INT_MIN + f; } SAm[p1 >> 1] = s1 - f;
            unsigned int p2 = (unsigned int)SA[i + 2]; int s2 = SAm[p2 >> 1]; if (s2 < 0) { T[p2] |= INT_MIN; f++; s2 = i + 2 + INT_MIN + f; } SAm[p2 >> 1] = s2 - f;
            unsigned int p3 = (unsigned int)SA[i + 3]; int s3 = SAm[p3 >> 1]; if (s3 < 0) { T[p3] |= INT_MIN; f++; s3 = i + 3 + INT_MIN + f; } SAm[p3 >> 1] = s3 - f;
        }

        for (j += 2 * (int)prefetch_distance + 3; i < j; i += 1)
        {
            unsigned int p = (unsigned int)SA[i]; int s = SAm[p >> 1]; if (s < 0) { T[p] |= INT_MIN; f++; s = i + INT_MIN + f; } SAm[p >> 1] = s - f;
        }
    }

    {
        int * RESTRICT SAl = &SA[0];
        int * RESTRICT SAr = &SA[0];

        ptrdiff_t i, j, l = (ptrdiff_t)m - 1, r = (ptrdiff_t)n + (ptrdiff_t)fs - 1;
        for (i = (ptrdiff_t)m + ((ptrdiff_t)n >> 1) - 1, j = (ptrdiff_t)m + 3; i >= j; i -= 4)
        {
            libsais_prefetch(&SA[i - prefetch_distance]);

            int p0 = SA[i - 0]; SAl[l] = p0 & INT_MAX; l -= p0 < 0; SAr[r] = p0 - 1; r -= p0 > 0;
            int p1 = SA[i - 1]; SAl[l] = p1 & INT_MAX; l -= p1 < 0; SAr[r] = p1 - 1; r -= p1 > 0;
            int p2 = SA[i - 2]; SAl[l] = p2 & INT_MAX; l -= p2 < 0; SAr[r] = p2 - 1; r -= p2 > 0;
            int p3 = SA[i - 3]; SAl[l] = p3 & INT_MAX; l -= p3 < 0; SAr[r] = p3 - 1; r -= p3 > 0;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            int p = SA[i]; SAl[l] = p & INT_MAX; l -= p < 0; SAr[r] = p - 1; r -= p > 0;
        }

        memcpy(&SA[(ptrdiff_t)n + (ptrdiff_t)fs - (ptrdiff_t)m], &SA[(ptrdiff_t)l + 1], (size_t)f * sizeof(int));
    }

    return f;
}

static void libsais_merge_compacted_lms_suffixes_32s(int * RESTRICT T, int * RESTRICT SA, int n, int m, int f)
{
    const ptrdiff_t prefetch_distance = 32;

    const int * RESTRICT SAnm = &SA[n - m - 1];

    {
        int i, j, l = 0, tmp = SAnm[l];
        for (i = 0, j = n - 6; i < j; i += 4)
        {
            libsais_prefetch(&T[i + prefetch_distance]);

            int c0 = T[i + 0]; if (c0 < 0) { T[i + 0] = c0 & INT_MAX; SA[tmp] = i + 0; i++; tmp = SAnm[++l]; }
            int c1 = T[i + 1]; if (c1 < 0) { T[i + 1] = c1 & INT_MAX; SA[tmp] = i + 1; i++; tmp = SAnm[++l]; }
            int c2 = T[i + 2]; if (c2 < 0) { T[i + 2] = c2 & INT_MAX; SA[tmp] = i + 2; i++; tmp = SAnm[++l]; }
            int c3 = T[i + 3]; if (c3 < 0) { T[i + 3] = c3 & INT_MAX; SA[tmp] = i + 3; i++; tmp = SAnm[++l]; }
        }

        for (j += 6; i < j; i += 1)
        {
            int c0 = T[i]; if (c0 < 0) { T[i] = c0 & INT_MAX; SA[tmp] = i; i++; tmp = SAnm[++l]; }
        }
    }

    {
        ptrdiff_t i, j, l = f; int tmp = SAnm[l];
        for (i = 0, j = (ptrdiff_t)m - 3; i < j; i += 4)
        {
            libsais_prefetch(&SA[i + prefetch_distance]);

            if (SA[i + 0] == 0) { SA[i + 0] = tmp; tmp = SAnm[++l]; }
            if (SA[i + 1] == 0) { SA[i + 1] = tmp; tmp = SAnm[++l]; }
            if (SA[i + 2] == 0) { SA[i + 2] = tmp; tmp = SAnm[++l]; }
            if (SA[i + 3] == 0) { SA[i + 3] = tmp; tmp = SAnm[++l]; }
        }

        for (j += 3; i < j; i += 1)
        {
            if (SA[i] == 0) { SA[i] = tmp; tmp = SAnm[++l]; }
        }
    }
}

static void libsais_reconstruct_compacted_lms_suffixes_32s_2k(int * RESTRICT T, int * RESTRICT SA, int n, int k, int m, int fs, int f, int * RESTRICT buckets)
{
    if (f > 0)
    {
        memcpy(&SA[n - m - 1], &SA[n + fs - m], (size_t)f * sizeof(int));

        libsais_count_and_gather_compacted_lms_suffixes_32s_2k(T, SA, n, k, buckets);
        libsais_reconstruct_lms_suffixes(SA, n, m - f);

        memcpy(&SA[n - m - 1 + f], &SA[0], ((size_t)m - (size_t)f) * sizeof(int));
        memset(&SA[0], 0, (size_t)m * sizeof(int));

        libsais_merge_compacted_lms_suffixes_32s(T, SA, n, m, f);
    }
    else
    {
        libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets);
        libsais_reconstruct_lms_suffixes(SA, n, m);
    }
}

static void libsais_reconstruct_compacted_lms_suffixes_32s_1k(int * RESTRICT T, int * RESTRICT SA, int n, int k, int m, int fs, int f)
{
    if (f > 0)
    {
        memmove(&SA[n - m - 1], &SA[n + fs - m], (size_t)f * sizeof(int));

        libsais_gather_compacted_lms_suffixes_32s(T, SA, n);
        libsais_reconstruct_lms_suffixes(SA, n, m - f);

        memcpy(&SA[n - m - 1 + f], &SA[0], ((size_t)m - (size_t)f) * sizeof(int));
        memset(&SA[0], 0, (size_t)m * sizeof(int));

        libsais_merge_compacted_lms_suffixes_32s(T, SA, n, m, f);
    }
    else
    {
        libsais_gather_lms_suffixes_32s(T, SA, n);
        libsais_reconstruct_lms_suffixes(SA, n, m);
    }
}

static int libsais_main_32s(int * RESTRICT T, int * RESTRICT SA, int n, int k, int fs)
{
    if (k > 0 && fs / k >= 6)
    {
        int alignment = (fs - 1024) / k >= 6 ? 1024 : 16;
        int * RESTRICT buckets = (fs - alignment) / k >= 6 ? (int *)libsais_align_up(&SA[n + fs - 6 * k - alignment], (size_t)alignment * sizeof(int)) : &SA[n + fs - 6 * k];

        int m = libsais_count_and_gather_lms_suffixes_32s_4k(T, SA, n, k, buckets);
        if (m > 1)
        {
            memset(SA, 0, ((size_t)n - (size_t)m) * sizeof(int));

            int first_lms_suffix    = SA[n - m];
            int left_suffixes_count = libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(T, k, buckets, first_lms_suffix);

            libsais_radix_sort_lms_suffixes_32s_2k(T, SA, n, m, &buckets[4 * k]);
            libsais_radix_sort_set_markers_32s(SA, k, &buckets[4 * k], INT_MIN);

            libsais_initialize_buckets_for_partial_sorting_32s_6k(T, k, buckets, first_lms_suffix, left_suffixes_count);
            libsais_induce_partial_order_32s_6k(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count);

            int names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_4k(SA, n, m);
            if (names < m)
            {
                int f = libsais_compact_lms_suffixes_32s(T, SA, n, m, fs);

                if (libsais_main_32s(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_2k(T, SA, n, k, m, fs, f, buckets);
            }
            else
            {
                libsais_count_lms_suffixes_32s_2k(T, n, k, buckets);
            }

            libsais_initialize_buckets_start_and_end_32s_4k(k, buckets);
            libsais_place_lms_suffixes_histogram_32s_4k(SA, n, k, m, buckets);
            libsais_induce_final_order_32s_4k(T, SA, n, k, buckets);
        }
        else
        {
            SA[0] = SA[n - 1];

            libsais_initialize_buckets_start_and_end_32s_6k(k, buckets);
            libsais_place_lms_suffixes_histogram_32s_6k(SA, n, k, m, buckets);
            libsais_induce_final_order_32s_6k(T, SA, n, k, buckets);
        }

        return 0;
    }
    else if (k > 0 && fs / k >= 4)
    {
        int alignment = (fs - 1024) / k >= 4 ? 1024 : 16;
        int * RESTRICT buckets = (fs - alignment) / k >= 4 ? (int *)libsais_align_up(&SA[n + fs - 4 * k - alignment], (size_t)alignment * sizeof(int)) : &SA[n + fs - 4 * k];

        int m = libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets);
        if (m > 1)
        {
            libsais_initialize_buckets_for_radix_and_partial_sorting_32s_4k(T, k, buckets, SA[n - m]);

            libsais_radix_sort_lms_suffixes_32s_2k(T, SA, n, m, &buckets[1]);
            libsais_radix_sort_set_markers_32s(SA, k, &buckets[1], SUFFIX_GROUP_MARKER);
            
            libsais_place_lms_suffixes_interval_32s_4k(SA, n, k, m - 1, buckets);
            libsais_induce_partial_order_32s_4k(T, SA, n, k, buckets);

            int names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_4k(SA, n, m);
            if (names < m)
            {
                int f = libsais_compact_lms_suffixes_32s(T, SA, n, m, fs);

                if (libsais_main_32s(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_2k(T, SA, n, k, m, fs, f, buckets);
            }
            else
            {
                libsais_count_lms_suffixes_32s_2k(T, n, k, buckets);
            }
        }
        else
        {
            SA[0] = SA[n - 1];
        }

        libsais_initialize_buckets_start_and_end_32s_4k(k, buckets);
        libsais_place_lms_suffixes_histogram_32s_4k(SA, n, k, m, buckets);
        libsais_induce_final_order_32s_4k(T, SA, n, k, buckets);

        return 0;
    }
    else if (k > 0 && fs / k >= 2)
    {
        int alignment = (fs - 1024) / k >= 2 ? 1024 : 16;
        int * RESTRICT buckets = (fs - alignment) / k >= 2 ? (int *)libsais_align_up(&SA[n + fs - 2 * k - alignment], (size_t)alignment * sizeof(int)) : &SA[n + fs - 2 * k];

        int m = libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets);
        if (m > 1)
        {
            libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_2k(T, k, buckets, SA[n - m]);

            libsais_radix_sort_lms_suffixes_32s_2k(T, SA, n, m, &buckets[1]);
            libsais_place_lms_suffixes_interval_32s_2k(SA, n, k, m - 1, buckets);

            libsais_initialize_buckets_start_and_end_32s_2k(k, buckets);
            libsais_induce_partial_order_32s_2k(T, SA, n, k, buckets);

            int names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_1k(T, SA, n, m);
            if (names < m)
            {
                int f = libsais_compact_lms_suffixes_32s(T, SA, n, m, fs);

                if (libsais_main_32s(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_2k(T, SA, n, k, m, fs, f, buckets);
            }
            else
            {
                libsais_count_lms_suffixes_32s_2k(T, n, k, buckets);
            }
        }
        else
        {
            SA[0] = SA[n - 1];
        }

        libsais_initialize_buckets_end_32s_2k(k, buckets);
        libsais_place_lms_suffixes_histogram_32s_2k(SA, n, k, m, buckets);

        libsais_initialize_buckets_start_and_end_32s_2k(k, buckets);
        libsais_induce_final_order_32s_2k(T, SA, n, k, buckets);

        return 0;
    }
    else
    {
        int * buffer = fs < k ? (int *)libsais_aligned_malloc((size_t)k * sizeof(int), 4096) : (int *)NULL;

        int alignment = fs - 1024 >= k ? 1024 : 16;
        int * RESTRICT buckets = fs - alignment >= k ? (int *)libsais_align_up(&SA[n + fs - k - alignment], (size_t)alignment * sizeof(int)) : fs >= k ? &SA[n + fs - k] : buffer;

        if (buckets == NULL) { return -2; }

        memset(SA, 0, (size_t)n * sizeof(int));

        libsais_count_suffixes_32s(T, n, k, buckets); 
        libsais_initialize_buckets_end_32s_1k(k, buckets);

        int m = libsais_radix_sort_lms_suffixes_32s_1k(T, SA, n, buckets);
        if (m > 1)
        {
            libsais_induce_partial_order_32s_1k(T, SA, n, k, buckets);

            int names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_1k(T, SA, n, m);
            if (names < m)
            {
                if (buffer != NULL) { libsais_aligned_free(buffer); buckets = NULL; }

                int f = libsais_compact_lms_suffixes_32s(T, SA, n, m, fs);

                if (libsais_main_32s(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_1k(T, SA, n, k, m, fs, f);

                if (buckets == NULL) { buckets = buffer = (int *)libsais_aligned_malloc((size_t)k * sizeof(int), 4096); }
                if (buckets == NULL) { return -2; }
            }
            
            libsais_count_suffixes_32s(T, n, k, buckets);
            libsais_initialize_buckets_end_32s_1k(k, buckets);
            libsais_place_lms_suffixes_interval_32s_1k(T, SA, n, k, m, buckets);
        }

        libsais_induce_final_order_32s_1k(T, SA, n, k, buckets);
        libsais_aligned_free(buffer);

        return 0;
    }
}

static int libsais_main_8u(const unsigned char * T, int * SA, int n, int bwt)
{
    int * RESTRICT buckets = (int *)libsais_aligned_malloc(8 * ALPHABET_SIZE * sizeof(int), 4096);

    if (buckets != NULL)
    {
        int m = libsais_count_and_gather_lms_suffixes_8u(T, SA, n, buckets);

        libsais_initialize_buckets_start_and_end_8u(buckets);

        if (m > 0)
        {
            int first_lms_suffix    = SA[n - m];
            int left_suffixes_count = libsais_initialize_buckets_for_lms_suffixes_radix_sort_8u(T, buckets, first_lms_suffix);

            libsais_radix_sort_lms_suffixes_8u(T, SA, n, m, buckets);
            libsais_initialize_buckets_for_partial_sorting_8u(T, buckets, first_lms_suffix, left_suffixes_count);
            libsais_induce_partial_order_8u(T, SA, n, buckets, first_lms_suffix, left_suffixes_count);

            int names = libsais_renumber_and_gather_lms_suffixes_8u(SA, n, m);
            if (names < m)
            {
                if (libsais_main_32s(SA + n - m, SA, m, names, n - 2 * m) != 0)
                {
                    libsais_aligned_free(buckets);
                    return -2;
                }

                libsais_gather_lms_suffixes_8u(T, SA, n);
                libsais_reconstruct_lms_suffixes(SA, n, m);
            }

            libsais_place_lms_suffixes_interval_8u(SA, n, m, buckets);
        }
        else
        {
            memset(SA, 0, (size_t)n * sizeof(int));
        }

        int index = libsais_induce_final_order_8u(T, SA, n, bwt, buckets);

        libsais_aligned_free(buckets);
        return index;
    }

    return -2;
}

static void libsais_bwt_copy_8u(unsigned char * RESTRICT U, int * RESTRICT A, int n)
{
    const ptrdiff_t prefetch_distance = 32;

    ptrdiff_t i, j;
    for (i = 0, j = (ptrdiff_t)n - 7; i < j; i += 8)
    {
        libsais_prefetch(&A[i + prefetch_distance]);

        U[i + 0] = (unsigned char)A[i + 0];
        U[i + 1] = (unsigned char)A[i + 1];
        U[i + 2] = (unsigned char)A[i + 2];
        U[i + 3] = (unsigned char)A[i + 3];
        U[i + 4] = (unsigned char)A[i + 4];
        U[i + 5] = (unsigned char)A[i + 5];
        U[i + 6] = (unsigned char)A[i + 6];
        U[i + 7] = (unsigned char)A[i + 7];
    }

    for (j += 7; i < j; i += 1)
    {
        U[i] = (unsigned char)A[i];
    }
}

int libsais(const unsigned char * T, int * SA, int n)
{
    if ((T == NULL) || (SA == NULL) || (n < 0))
    {
        return -1;
    }
    else if (n < 2)
    {
        if (n == 1) { SA[0] = 0; }
        return 0;
    }

    return libsais_main_8u(T, SA, n, 0);
}

int libsais_bwt(const unsigned char * T, unsigned char * U, int * A, int n)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0)) 
    { 
        return -1; 
    }
    else if (n <= 1) 
    { 
        if (n == 1) { U[0] = T[0]; }
        return n; 
    }

    int index = libsais_main_8u(T, A, n, 1);
    if (index >= 0) 
    { 
        U[0] = T[n - 1];
        libsais_bwt_copy_8u(U + 1, A, index);
        libsais_bwt_copy_8u(U + 1 + index, A + 1 + index, n - index - 1);

        index++;
    }

    return index;
}
