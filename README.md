# libsais

libsais is a library for fast (see [Benchmarks](#benchmarks) below) linear time suffix array and Burrows-Wheeler transform construction based on induced sorting algorithm described in the following paper: Ge Nong, Sen Zhang and Wai Hong Chan, *Two Efficient Algorithms for Linear Suffix Array Construction*, 2008.

Copyright (c) 2021 Ilya Grebnov <ilya.grebnov@gmail.com>

>This program is based on the work of [sais-lite](https://sites.google.com/site/yuta256/sais) library by Yuta Mori.

## Introduction
libsais library provides simple C99 API to construct suffix array and Burrows-Wheeler transformed string from a given string over constant-size alphabet. The algorithm runs in a linear time using typically only ~12KB of extra memory with 2n bytes as absolute worst-case extra working space, where n is the length of the string.

> Note, libsais library is heavily dependent on fast memory and software prefetching, so it might not be suitable for embeded workloads or heavily concurrent workloads which could saturate dual channel RAM.

## License
libsais library is released under the [Apache License Version 2.0](LICENSE "Apache license")

## Changes
* February 23, 2021
  * Initial release.

## API
```c
    /**
    * Constructs the suffix array of a given string.
    * @param T [0..n-1] The input string.
    * @param SA [0..n-1] The output array of suffixes.
    * @param n The length of the given string.
    * @return 0 if no error occurred, -1 or -2 otherwise.
    */
    int libsais(const unsigned char * T, int * SA, int n);

    /**
    * Constructs the burrows-wheeler transformed string of a given string.
    * @param T [0..n-1] The input string.
    * @param U [0..n-1] The output string. (can be T)
    * @param A [0..n-1] The temporary array.
    * @param n The length of the given string.
    * @return The primary index if no error occurred, -1 or -2 otherwise.
    */
    int libsais_bwt(const unsigned char * T, unsigned char * U, int * A, int n);
```

---

# Benchmarks
  * CPU: Intel Core i7-9700K Processor (12M Cache, 5 GHz OC)
  * RAM: 16 GB Dual Channel 4133 MHz (17-17-17-37-400-2T OC)
  * OS: Microsoft Windows 10 Pro 64 Bit
  * Compiler: Clang 10.0.0 '-O3 -DNDEBUG'

The times are the minimum of five runs measuring single-threaded performance of suffix array construction.

### Silesia Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|         dickens |   10192446      |   0.202 sec (  50.43 MB/s) |   0.448 sec (  22.78 MB/s) |**+121.43%**|
|         mozilla |   51220480      |   1.120 sec (  45.73 MB/s) |   1.803 sec (  28.41 MB/s) | **+60.93%**|
|              mr |    9970564      |   0.185 sec (  53.90 MB/s) |   0.374 sec (  26.66 MB/s) |**+102.23%**|
|             nci |   33553445      |   0.625 sec (  53.70 MB/s) |   1.184 sec (  28.35 MB/s) | **+89.42%**|
|         ooffice |    6152192      |   0.116 sec (  53.03 MB/s) |   0.178 sec (  34.61 MB/s) | **+53.25%**|
|            osdb |   10085684      |   0.228 sec (  44.31 MB/s) |   0.368 sec (  27.38 MB/s) | **+61.83%**|
|         reymont |    6627202      |   0.123 sec (  53.78 MB/s) |   0.253 sec (  26.15 MB/s) |**+105.65%**|
|           samba |   21606400      |   0.407 sec (  53.06 MB/s) |   0.663 sec (  32.57 MB/s) | **+62.93%**|
|             sao |    7251944      |   0.181 sec (  40.04 MB/s) |   0.232 sec (  31.29 MB/s) | **+27.99%**|
|         webster |   41458703      |   1.019 sec (  40.70 MB/s) |   2.195 sec (  18.89 MB/s) |**+115.49%**|
|             xml |    5345280      |   0.079 sec (  68.05 MB/s) |   0.139 sec (  38.36 MB/s) | **+77.39%**|

### Large Canterbury Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|       bible.txt |    4047392      |   0.066 sec (  61.53 MB/s) |   0.142 sec (  28.48 MB/s) |**+116.00%**|
|          E.coli |    4638690      |   0.078 sec (  59.12 MB/s) |   0.203 sec (  22.83 MB/s) |**+158.97%**|
|    world192.txt |    2473400      |   0.037 sec (  67.32 MB/s) |   0.075 sec (  32.92 MB/s) |**+104.48%**|

### Manzini's Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|       chr22.dna |   34553758      |   0.772 sec (  44.76 MB/s) |   2.050 sec (  16.86 MB/s) |**+165.50%**|
|         etext99 |  105277340      |   3.190 sec (  33.01 MB/s) |   7.091 sec (  14.85 MB/s) |**+122.30%**|
|           howto |   39422105      |   0.926 sec (  42.58 MB/s) |   1.890 sec (  20.86 MB/s) |**+104.08%**|
|          jdk13c |   69728899      |   1.641 sec (  42.50 MB/s) |   3.086 sec (  22.60 MB/s) | **+88.08%**|
| linux-2.4.5.tar |  116254720      |   2.744 sec (  42.36 MB/s) |   5.073 sec (  22.92 MB/s) | **+84.85%**|
|        rctail96 |  114711151      |   3.183 sec (  36.04 MB/s) |   6.512 sec (  17.62 MB/s) |**+104.56%**|
|             rfc |  116421901      |   2.919 sec (  39.88 MB/s) |   5.710 sec (  20.39 MB/s) | **+95.60%**|
|     sprot34.dat |  109617186      |   3.010 sec (  36.42 MB/s) |   6.408 sec (  17.11 MB/s) |**+112.90%**|
|            w3c2 |  104201579      |   2.524 sec (  41.29 MB/s) |   4.612 sec (  22.59 MB/s) | **+82.72%**|

### Large Text Compression Benchmark Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|          enwik8 |  100000000      |   2.878 sec (  34.74 MB/s) |   6.216 sec (  16.09 MB/s) |**+115.96%**|
|          enwik9 | 1000000000      |  36.786 sec (  27.18 MB/s) |  77.817 sec (  12.85 MB/s) |**+111.54%**|

### The Gauntlet Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|            abac |     200000      |   0.003 sec (  75.94 MB/s) |   0.002 sec ( 129.08 MB/s) |    -41.17% |
|            abba |   10500596      |   0.167 sec (  62.77 MB/s) |   0.444 sec (  23.66 MB/s) |**+165.35%**|
|        book1x20 |   15375420      |   0.323 sec (  47.55 MB/s) |   0.719 sec (  21.38 MB/s) |**+122.37%**|
|   fib_s14930352 |   14930352      |   0.344 sec (  43.38 MB/s) |   1.119 sec (  13.34 MB/s) |**+225.11%**|
|           fss10 |   12078908      |   0.240 sec (  50.22 MB/s) |   0.850 sec (  14.20 MB/s) |**+253.60%**|
|            fss9 |    2851443      |   0.041 sec (  69.53 MB/s) |   0.128 sec (  22.25 MB/s) |**+212.42%**|
|         houston |    3839141      |   0.037 sec ( 104.15 MB/s) |   0.022 sec ( 177.07 MB/s) |    -41.18% |
|       paper5x80 |     956322      |   0.014 sec (  69.76 MB/s) |   0.024 sec (  40.22 MB/s) | **+73.45%**|
|           test1 |    2097152      |   0.041 sec (  51.06 MB/s) |   0.073 sec (  28.75 MB/s) | **+77.63%**|
|           test2 |    2097152      |   0.039 sec (  54.12 MB/s) |   0.049 sec (  42.85 MB/s) | **+26.29%**|
|           test3 |    2097088      |   0.036 sec (  58.05 MB/s) |   0.045 sec (  47.04 MB/s) | **+23.41%**|

### Pizza & Chilli Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|        dblp.xml |  296135874      |   8.628 sec (  34.32 MB/s) |  16.145 sec (  18.34 MB/s) | **+87.11%**|
|             dna |  403927746      |  14.983 sec (  26.96 MB/s) |  37.108 sec (  10.89 MB/s) |**+147.67%**|
|  english.1024MB | 1073741824      |  46.777 sec (  22.95 MB/s) |  99.729 sec (  10.77 MB/s) |**+113.20%**|
|         pitches |   55832855      |   1.253 sec (  44.56 MB/s) |   2.390 sec (  23.36 MB/s) | **+90.74%**|
|        proteins | 1184051855      |  59.239 sec (  19.99 MB/s) | 120.934 sec (   9.79 MB/s) |**+104.15%**|
|         sources |  210866607      |   5.407 sec (  39.00 MB/s) |  10.227 sec (  20.62 MB/s) | **+89.13%**|

### Pizza & Chilli Repetitive Corpus ###

|  file           |      size       |       libsais 1.0.0        |       divsufsort 2.0.2     | improvement |
|:---------------:|:---------------:|:--------------------------:|:--------------------------:|:-----------:|
|            cere |  461286644      |  16.447 sec (  28.05 MB/s) |  35.326 sec (  13.06 MB/s) |**+114.79%**|
|       coreutils |  205281778      |   5.260 sec (  39.03 MB/s) |  10.842 sec (  18.93 MB/s) |**+106.12%**|
| einstein.de.txt |   92758441      |   2.493 sec (  37.21 MB/s) |   4.696 sec (  19.75 MB/s) | **+88.41%**|
| einstein.en.txt |  467626544      |  15.487 sec (  30.19 MB/s) |  32.412 sec (  14.43 MB/s) |**+109.28%**|
|Escherichia_Coli |  112689515      |   3.381 sec (  33.33 MB/s) |   7.690 sec (  14.65 MB/s) |**+127.44%**|
|       influenza |  154808555      |   3.843 sec (  40.28 MB/s) |   9.292 sec (  16.66 MB/s) |**+141.76%**|
|          kernel |  257961616      |   7.039 sec (  36.65 MB/s) |  13.718 sec (  18.80 MB/s) | **+94.88%**|
|            para |  429265758      |  14.429 sec (  29.75 MB/s) |  33.080 sec (  12.98 MB/s) |**+129.26%**|
|   world_leaders |   46968181      |   0.881 sec (  53.34 MB/s) |   1.313 sec (  35.76 MB/s) | **+49.17%**|



