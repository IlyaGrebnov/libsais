# libsais

The libsais is a library for fast (see [Benchmarks](#benchmarks) below) linear time suffix array and Burrows-Wheeler transform construction based on induced sorting algorithm described in the following paper: Ge Nong, Sen Zhang and Wai Hong Chan, *Two Efficient Algorithms for Linear Suffix Array Construction*, 2008.

Copyright (c) 2021 Ilya Grebnov <ilya.grebnov@gmail.com>

>The libsais is inspired by [libdivsufsort](https://github.com/y-256/libdivsufsort) and [sais](https://sites.google.com/site/yuta256/sais) libraries by Yuta Mori.

## Introduction
The libsais provides simple C99 API to construct suffix array and Burrows-Wheeler transformed string from a given string over constant-size alphabet. The algorithm runs in a linear time using typically only ~16KB of extra memory (with 2n bytes as absolute worst-case; where n is the length of the string). OpenMP acceleration uses 132KB of addition memory per thread.

> * The libsais works with compilers from GNU, Microsoft and Intel, but I recommend Clang for best performance.
> * The libsais is sensitive to fast memory and software prefetching and might not be suitable for some workloads. Please benchmark yourself.

## License
The libsais is released under the [Apache License Version 2.0](LICENSE "Apache license")

## Changes
* April 4, 2021 (2.0.0)
  * OpenMP acceleration. 
* February 23, 2021 (1.0.0)
  * Initial release.

## API
```c
    /**
    * Constructs the suffix array of a given string.
    * @param T [0..n-1] The input string.
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the given string.
    * @param fs Extra space available at the end of SA array (can be 0).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    int libsais(const unsigned char * T, int * SA, int n, int fs);

#if defined(_OPENMP)
    /**
    * Constructs the suffix array of a given string in parallel using OpenMP.
    * @param T [0..n-1] The input string.
    * @param SA [0..n-1+fs] The output array of suffixes.
    * @param n The length of the given string.
    * @param fs Extra space available at the end of SA array (can be 0).
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    int libsais_omp(const unsigned char * T, int * SA, int n, int fs, int threads);
#endif

    /**
    * Constructs the burrows-wheeler transformed string of a given string.
    * @param T [0..n-1] The input string.
    * @param U [0..n-1] The output string (can be T).
    * @param A [0..n-1+fs] The temporary array.
    * @param n The length of the given string.
    * @param fs Extra space available at the end of A array (can be 0).
    * @return The primary index if no error occurred, -1 or -2 otherwise.
    */
    int libsais_bwt(const unsigned char * T, unsigned char * U, int * A, int n, int fs);

#if defined(_OPENMP)
    /**
    * Constructs the burrows-wheeler transformed string of a given string in parallel using OpenMP.
    * @param T [0..n-1] The input string.
    * @param U [0..n-1] The output string (can be T).
    * @param A [0..n-1+fs] The temporary array.
    * @param n The length of the given string.
    * @param fs Extra space available at the end of A array (can be 0).
    * @param threads The number of OpenMP threads to use (can be 0 for OpenMP default).
    * @return The primary index if no error occurred, -1 or -2 otherwise.
    */
    int libsais_bwt_omp(const unsigned char * T, unsigned char * U, int * A, int n, int fs, int threads);
#endif
```

---

# Benchmarks

Benchmarks are moved to own [Benchmarks.md](Benchmarks.md) file.

## Large page and multi-core systems support

Large-pages and OpenMP improves the libsais performance, but multi-core scalability is limited to 4-5 cores. Here is an example of such improvements on Manzini Corpus.

> * Multi-core scalability is limited by RAM performance. I hope, DDR5 and more RAM channels improve this in future.

| file            |      size | baseline| LP      | LP w 2c | LP w 3c | LP w 4c | LP w 5c | LP w 6c | LP w 7c | LP w 8c  |
|:---------------:|:---------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:--------:|
| chr22.dna       | 34553758  |42.90MB/s|49.52MB/s|60.71MB/s|69.95MB/s|73.70MB/s|75.44MB/s|75.87MB/s|76.06MB/s|75.37 MB/s|
| etext99         | 105277340 |31.98MB/s|40.20MB/s|47.54MB/s|52.90MB/s|55.72MB/s|56.57MB/s|56.67MB/s|56.63MB/s|56.18 MB/s|
| gcc-3.0.tar     | 86630400  |43.95MB/s|49.56MB/s|57.72MB/s|65.56MB/s|69.17MB/s|70.71MB/s|70.77MB/s|70.99MB/s|70.72 MB/s|
| howto           | 39422105  |42.41MB/s|47.49MB/s|55.89MB/s|62.71MB/s|66.62MB/s|67.55MB/s|68.15MB/s|68.60MB/s|68.21 MB/s|
| jdk13c          | 69728899  |42.11MB/s|47.10MB/s|53.92MB/s|61.51MB/s|65.88MB/s|67.02MB/s|67.65MB/s|67.54MB/s|67.20 MB/s|
| linux-2.4.5.tar | 116254720 |41.76MB/s|48.55MB/s|56.71MB/s|63.68MB/s|67.61MB/s|68.56MB/s|69.21MB/s|69.26MB/s|68.76 MB/s|
| rctail96        | 114711151 |35.81MB/s|42.74MB/s|50.13MB/s|56.83MB/s|59.61MB/s|60.57MB/s|60.79MB/s|60.59MB/s|60.26 MB/s|
| rfc             | 116421901 |39.47MB/s|46.36MB/s|54.09MB/s|61.37MB/s|65.34MB/s|66.09MB/s|66.35MB/s|66.55MB/s|66.11 MB/s|
| sprot34.dat     | 109617186 |35.29MB/s|44.19MB/s|52.25MB/s|58.77MB/s|61.43MB/s|63.07MB/s|62.86MB/s|62.87MB/s|62.59 MB/s|
| w3c2            | 104201579 |41.32MB/s|46.22MB/s|53.67MB/s|60.97MB/s|63.47MB/s|65.35MB/s|65.39MB/s|65.16MB/s|65.05 MB/s|

## Additional memory

The libsais reuses free space allocated for suffix array during induction. Sometimes this space might not be sufficient for most efficient algorithm and libsais will need to fallback to less efficient one (libsais has 4 algorithms at different break points). You could avoid this fallbacks and improve performance by allocating additional space at the end of suffix array.

> * All other files from [Benchmarks](#benchmarks) above do not suffer from this fallbacks.

|  file           |    size     |     libsais + O(n)  (ST)   |     libsais + O(1) (ST)    |speedup (ST)|    libsais + O(n)  (MT)    |     libsais + O(1) (ST)    |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            osdb |    10085684 |   0.224 sec (  44.99 MB/s) |   0.234 sec (  43.15 MB/s) |  **+4.26%**|   0.155 sec (  64.98 MB/s) |   0.166 sec (  60.60 MB/s) |  **+7.23%**|
|           x-ray |     8474240 |   0.192 sec (  44.25 MB/s) |   0.223 sec (  37.92 MB/s) | **+16.67%**|   0.125 sec (  67.67 MB/s) |   0.159 sec (  53.26 MB/s) | **+27.07%**|
|             sao |     7251944 |   0.178 sec (  40.82 MB/s) |   0.185 sec (  39.20 MB/s) |  **+4.12%**|   0.138 sec (  52.60 MB/s) |   0.150 sec (  48.44 MB/s) |  **+8.58%**|
|         ooffice |     6152192 |   0.113 sec (  54.51 MB/s) |   0.119 sec (  51.51 MB/s) |  **+5.81%**|   0.083 sec (  74.44 MB/s) |   0.090 sec (  68.27 MB/s) |  **+9.03%**|
|            abac |      200000 |   0.002 sec (  83.78 MB/s) |   0.003 sec (  73.56 MB/s) | **+13.89%**|   0.002 sec ( 114.44 MB/s) |   0.002 sec (  94.68 MB/s) | **+20.87%**|
|           test3 |     2097088 |   0.035 sec (  60.13 MB/s) |   0.038 sec (  55.44 MB/s) |  **+8.46%**|   0.027 sec (  76.72 MB/s) |   0.031 sec (  67.60 MB/s) | **+13.49%**|
