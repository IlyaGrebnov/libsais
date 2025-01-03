# Specifications
  * CPU: Intel Core i7-9700K Processor (12M Cache, 5 GHz all cores)
  * RAM: 2x8 GB dual-channel DDR4 (4133 MHz, 17-17-17-37-400-2T with tight subtimings)
  * OS: Microsoft Windows 10 Pro 64 Bit (with MEM_LARGE_PAGES support enabled)
  * Compiler: Clang 11.0.0 '-Ofast -march=skylake -fopenmp -DNDEBUG'

The times are the minimum of five runs measuring **single-threaded (ST)** and **multi-threaded (MT)** performance of suffix array construction.

### [Silesia Corpus](https://www.data-compression.info/Corpora/SilesiaCorpus/index.html) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|         dickens |    10192446 |   0.204 sec (  50.00 MB/s) |   0.446 sec (  22.87 MB/s) |**+118.64%**|   0.120 sec (  84.98 MB/s) |   0.296 sec (  34.42 MB/s) |**+146.88%**|
|         mozilla |    51220480 |   1.105 sec (  46.34 MB/s) |   1.819 sec (  28.15 MB/s) | **+64.60%**|   0.712 sec (  71.96 MB/s) |   1.206 sec (  42.47 MB/s) | **+69.44%**|
|              mr |     9970564 |   0.187 sec (  53.21 MB/s) |   0.370 sec (  26.92 MB/s) | **+97.66%**|   0.131 sec (  76.03 MB/s) |   0.261 sec (  38.22 MB/s) | **+98.95%**|
|             nci |    33553445 |   0.641 sec (  52.38 MB/s) |   1.151 sec (  29.16 MB/s) | **+79.59%**|   0.364 sec (  92.20 MB/s) |   1.036 sec (  32.39 MB/s) |**+184.65%**|
|         ooffice |     6152192 |   0.116 sec (  53.23 MB/s) |   0.180 sec (  34.26 MB/s) | **+55.38%**|   0.082 sec (  75.06 MB/s) |   0.108 sec (  57.08 MB/s) | **+31.51%**|
|            osdb |    10085684 |   0.224 sec (  45.02 MB/s) |   0.361 sec (  27.94 MB/s) | **+61.12%**|   0.154 sec (  65.61 MB/s) |   0.272 sec (  37.07 MB/s) | **+77.00%**|
|         reymont |     6627202 |   0.125 sec (  53.09 MB/s) |   0.249 sec (  26.61 MB/s) | **+99.49%**|   0.076 sec (  87.44 MB/s) |   0.183 sec (  36.22 MB/s) |**+141.39%**|
|           samba |    21606400 |   0.401 sec (  53.82 MB/s) |   0.661 sec (  32.70 MB/s) | **+64.59%**|   0.274 sec (  79.00 MB/s) |   0.467 sec (  46.27 MB/s) | **+70.74%**|
|             sao |     7251944 |   0.181 sec (  40.02 MB/s) |   0.232 sec (  31.24 MB/s) | **+28.10%**|   0.132 sec (  55.02 MB/s) |   0.143 sec (  50.63 MB/s) |  **+8.67%**|
|         webster |    41458703 |   1.020 sec (  40.64 MB/s) |   2.190 sec (  18.93 MB/s) |**+114.70%**|   0.595 sec (  69.72 MB/s) |   1.653 sec (  25.08 MB/s) |**+177.98%**|
|           x-ray |     8474240 |   0.215 sec (  39.35 MB/s) |   0.337 sec (  25.13 MB/s) | **+56.58%**|   0.153 sec (  55.40 MB/s) |   0.184 sec (  46.10 MB/s) | **+20.19%**|
|             xml |     5345280 |   0.079 sec (  67.54 MB/s) |   0.138 sec (  38.87 MB/s) | **+73.78%**|   0.050 sec ( 106.50 MB/s) |   0.108 sec (  49.52 MB/s) |**+115.06%**|



### [Large Canterbury Corpus](https://www.data-compression.info/Corpora/CanterburyCorpus/) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       bible.txt |     4047392 |   0.066 sec (  61.40 MB/s) |   0.140 sec (  28.95 MB/s) |**+112.14%**|   0.038 sec ( 106.74 MB/s) |   0.100 sec (  40.58 MB/s) |**+163.00%**|
|          E.coli |     4638690 |   0.079 sec (  58.74 MB/s) |   0.199 sec (  23.35 MB/s) |**+151.59%**|   0.038 sec ( 121.18 MB/s) |   0.154 sec (  30.15 MB/s) |**+301.95%**|
|    world192.txt |     2473400 |   0.037 sec (  66.72 MB/s) |   0.075 sec (  32.92 MB/s) |**+102.65%**|   0.023 sec ( 108.02 MB/s) |   0.050 sec (  49.58 MB/s) |**+117.85%**|



### [Manzini Corpus](https://people.unipmn.it/manzini/lightweight/corpus/) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       chr22.dna |    34553758 |   0.778 sec (  44.43 MB/s) |   2.005 sec (  17.23 MB/s) |**+157.85%**|   0.436 sec (  79.27 MB/s) |   1.676 sec (  20.61 MB/s) |**+284.61%**|
|         etext99 |   105277340 |   3.174 sec (  33.17 MB/s) |   7.083 sec (  14.86 MB/s) |**+123.15%**|   1.739 sec (  60.53 MB/s) |   5.345 sec (  19.70 MB/s) |**+207.33%**|
|     gcc-3.0.tar |    86630400 |   1.930 sec (  44.89 MB/s) |   3.458 sec (  25.05 MB/s) | **+79.18%**|   1.188 sec (  72.90 MB/s) |   2.508 sec (  34.54 MB/s) |**+111.04%**|
|           howto |    39422105 |   0.912 sec (  43.24 MB/s) |   1.884 sec (  20.93 MB/s) |**+106.60%**|   0.557 sec (  70.80 MB/s) |   1.279 sec (  30.83 MB/s) |**+129.64%**|
|          jdk13c |    69728899 |   1.630 sec (  42.79 MB/s) |   3.019 sec (  23.10 MB/s) | **+85.25%**|   1.008 sec (  69.20 MB/s) |   2.490 sec (  28.00 MB/s) |**+147.11%**|
| linux-2.4.5.tar |   116254720 |   2.688 sec (  43.25 MB/s) |   4.980 sec (  23.34 MB/s) | **+85.27%**|   1.651 sec (  70.40 MB/s) |   3.540 sec (  32.84 MB/s) |**+114.34%**|
|        rctail96 |   114711151 |   3.162 sec (  36.27 MB/s) |   6.370 sec (  18.01 MB/s) |**+101.41%**|   1.813 sec (  63.27 MB/s) |   5.009 sec (  22.90 MB/s) |**+176.26%**|
|             rfc |   116421901 |   2.929 sec (  39.75 MB/s) |   5.689 sec (  20.46 MB/s) | **+94.23%**|   1.681 sec (  69.24 MB/s) |   4.182 sec (  27.84 MB/s) |**+148.70%**|
|     sprot34.dat |   109617186 |   2.983 sec (  36.75 MB/s) |   6.324 sec (  17.33 MB/s) |**+112.03%**|   1.712 sec (  64.03 MB/s) |   4.509 sec (  24.31 MB/s) |**+163.38%**|
|            w3c2 |   104201579 |   2.450 sec (  42.53 MB/s) |   4.426 sec (  23.54 MB/s) | **+80.67%**|   1.565 sec (  66.57 MB/s) |   3.570 sec (  29.19 MB/s) |**+128.10%**|



### [Large Text Compression Benchmark Corpus](https://www.mattmahoney.net/dc/textdata.html) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|          enwik8 |   100000000 |   2.841 sec (  35.20 MB/s) |   6.209 sec (  16.11 MB/s) |**+118.58%**|   1.604 sec (  62.36 MB/s) |   4.446 sec (  22.49 MB/s) |**+177.24%**|
|          enwik9 |  1000000000 |  40.718 sec (  24.56 MB/s) |  82.407 sec (  12.13 MB/s) |**+102.39%**|  19.138 sec (  52.25 MB/s) |  63.373 sec (  15.78 MB/s) |**+231.14%**|



### [The Gauntlet Corpus](https://github.com/michaelmaniscalco/gauntlet_corpus) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            abac |      200000 |   0.003 sec (  73.23 MB/s) |   0.002 sec ( 130.08 MB/s) |   -43.71%  |   0.002 sec (  94.65 MB/s) |   0.002 sec ( 125.78 MB/s) |   -24.75%  |
|            abba |    10500596 |   0.169 sec (  62.02 MB/s) |   0.431 sec (  24.36 MB/s) |**+154.64%**|   0.090 sec ( 116.53 MB/s) |   0.442 sec (  23.76 MB/s) |**+390.38%**|
|        book1x20 |    15375420 |   0.328 sec (  46.85 MB/s) |   0.720 sec (  21.36 MB/s) |**+119.34%**|   0.195 sec (  79.02 MB/s) |   0.505 sec (  30.47 MB/s) |**+159.34%**|
|   fib_s14930352 |    14930352 |   0.341 sec (  43.80 MB/s) |   1.070 sec (  13.95 MB/s) |**+213.96%**|   0.173 sec (  86.06 MB/s) |   1.059 sec (  14.10 MB/s) |**+510.23%**|
|           fss10 |    12078908 |   0.239 sec (  50.54 MB/s) |   0.811 sec (  14.89 MB/s) |**+239.49%**|   0.135 sec (  89.17 MB/s) |   0.793 sec (  15.24 MB/s) |**+485.07%**|
|            fss9 |     2851443 |   0.041 sec (  68.79 MB/s) |   0.123 sec (  23.12 MB/s) |**+197.52%**|   0.025 sec ( 116.35 MB/s) |   0.121 sec (  23.59 MB/s) |**+393.28%**|
|         houston |     3839141 |   0.037 sec ( 103.92 MB/s) |   0.024 sec ( 162.18 MB/s) |   -35.92%  |   0.018 sec ( 209.58 MB/s) |   0.023 sec ( 164.72 MB/s) | **+27.24%**|
|       paper5x80 |      956322 |   0.014 sec (  69.13 MB/s) |   0.024 sec (  40.00 MB/s) | **+72.82%**|   0.008 sec ( 112.87 MB/s) |   0.020 sec (  47.87 MB/s) |**+135.81%**|
|           test1 |     2097152 |   0.039 sec (  53.45 MB/s) |   0.084 sec (  24.82 MB/s) |**+115.35%**|   0.020 sec ( 106.42 MB/s) |   0.080 sec (  26.25 MB/s) |**+305.39%**|
|           test2 |     2097152 |   0.039 sec (  53.87 MB/s) |   0.058 sec (  36.02 MB/s) | **+49.56%**|   0.019 sec ( 109.39 MB/s) |   0.052 sec (  40.23 MB/s) |**+171.92%**|
|           test3 |     2097088 |   0.037 sec (  56.56 MB/s) |   0.046 sec (  45.93 MB/s) | **+23.13%**|   0.031 sec (  67.13 MB/s) |   0.047 sec (  44.81 MB/s) | **+49.79%**|



### [Pizza & Chilli Corpus](https://pizzachili.dcc.uchile.cl/texts.html) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|        dblp.xml |   296135874 |   8.428 sec (  35.14 MB/s) |  16.001 sec (  18.51 MB/s) | **+89.85%**|   4.749 sec (  62.36 MB/s) |  12.932 sec (  22.90 MB/s) |**+172.29%**|
|             dna |   403927746 |  14.257 sec (  28.33 MB/s) |  36.002 sec (  11.22 MB/s) |**+152.53%**|   6.462 sec (  62.51 MB/s) |  32.014 sec (  12.62 MB/s) |**+395.42%**|
|  english.1024MB |  1073741824 |  51.568 sec (  20.82 MB/s) | 104.407 sec (  10.28 MB/s) |**+102.47%**|  22.777 sec (  47.14 MB/s) |  84.357 sec (  12.73 MB/s) |**+270.37%**|
|         pitches |    55832855 |   1.222 sec (  45.70 MB/s) |   2.380 sec (  23.46 MB/s) | **+94.78%**|   0.820 sec (  68.10 MB/s) |   1.422 sec (  39.26 MB/s) | **+73.45%**|
|        proteins |  1184051855 |  60.105 sec (  19.70 MB/s) | 124.688 sec (   9.50 MB/s) |**+107.45%**|  25.981 sec (  45.57 MB/s) |  79.149 sec (  14.96 MB/s) |**+204.65%**|
|         sources |   210866607 |   5.381 sec (  39.19 MB/s) |  10.210 sec (  20.65 MB/s) | **+89.74%**|   3.223 sec (  65.42 MB/s) |   7.350 sec (  28.69 MB/s) |**+128.03%**|



### [Pizza & Chilli Repetitive Corpus](https://pizzachili.dcc.uchile.cl/repcorpus.html) ###

|  file           |    size     |     libsais 2.1.0  (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.1.0  (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            cere |   461286644 |  16.184 sec (  28.50 MB/s) |  34.422 sec (  13.40 MB/s) |**+112.69%**|   7.349 sec (  62.77 MB/s) |  30.502 sec (  15.12 MB/s) |**+315.04%**|
|       coreutils |   205281778 |   5.264 sec (  39.00 MB/s) |  10.506 sec (  19.54 MB/s) | **+99.58%**|   3.090 sec (  66.43 MB/s) |   8.170 sec (  25.13 MB/s) |**+164.39%**|
| einstein.de.txt |    92758441 |   2.449 sec (  37.87 MB/s) |   4.582 sec (  20.25 MB/s) | **+87.06%**|   1.381 sec (  67.19 MB/s) |   3.744 sec (  24.78 MB/s) |**+171.18%**|
| einstein.en.txt |   467626544 |  15.279 sec (  30.61 MB/s) |  31.504 sec (  14.84 MB/s) |**+106.20%**|   7.902 sec (  59.18 MB/s) |  27.152 sec (  17.22 MB/s) |**+243.60%**|
|Escherichia_Coli |   112689515 |   3.360 sec (  33.54 MB/s) |   7.495 sec (  15.03 MB/s) |**+123.06%**|   1.683 sec (  66.94 MB/s) |   6.283 sec (  17.94 MB/s) |**+273.25%**|
|       influenza |   154808555 |   3.816 sec (  40.57 MB/s) |   8.982 sec (  17.24 MB/s) |**+135.38%**|   2.209 sec (  70.07 MB/s) |   7.757 sec (  19.96 MB/s) |**+251.11%**|
|          kernel |   257961616 |   7.091 sec (  36.38 MB/s) |  13.511 sec (  19.09 MB/s) | **+90.54%**|   3.866 sec (  66.72 MB/s) |  10.158 sec (  25.40 MB/s) |**+162.74%**|
|            para |   429265758 |  15.407 sec (  27.86 MB/s) |  32.821 sec (  13.08 MB/s) |**+113.03%**|   6.905 sec (  62.17 MB/s) |  29.143 sec (  14.73 MB/s) |**+322.06%**|
|   world_leaders |    46968181 |   0.880 sec (  53.35 MB/s) |   1.301 sec (  36.11 MB/s) | **+47.73%**|   0.460 sec ( 102.16 MB/s) |   1.054 sec (  44.58 MB/s) |**+129.18%**|
|dblp.xml.00001.1 |   104857600 |   4.445 sec (  23.59 MB/s) |   5.913 sec (  17.73 MB/s) | **+33.03%**|   1.976 sec (  53.07 MB/s) |   5.130 sec (  20.44 MB/s) |**+159.62%**|
|dblp.xml.00001.2 |   104857600 |   4.423 sec (  23.71 MB/s) |   6.050 sec (  17.33 MB/s) | **+36.78%**|   1.967 sec (  53.32 MB/s) |   5.217 sec (  20.10 MB/s) |**+165.27%**|
| dblp.xml.0001.1 |   104857600 |   3.062 sec (  34.24 MB/s) |   5.534 sec (  18.95 MB/s) | **+80.70%**|   1.427 sec (  73.50 MB/s) |   4.664 sec (  22.48 MB/s) |**+226.91%**|
| dblp.xml.0001.2 |   104857600 |   3.118 sec (  33.63 MB/s) |   5.592 sec (  18.75 MB/s) | **+79.32%**|   1.664 sec (  63.00 MB/s) |   4.763 sec (  22.01 MB/s) |**+186.21%**|
|       dna.001.1 |   104857600 |   3.090 sec (  33.94 MB/s) |   6.368 sec (  16.47 MB/s) |**+106.10%**|   1.539 sec (  68.13 MB/s) |   5.338 sec (  19.64 MB/s) |**+246.87%**|
|   english.001.2 |   104857600 |   3.445 sec (  30.44 MB/s) |   5.920 sec (  17.71 MB/s) | **+71.87%**|   1.745 sec (  60.08 MB/s) |   4.268 sec (  24.57 MB/s) |**+144.55%**|
|  proteins.001.1 |   104857600 |   3.406 sec (  30.79 MB/s) |   5.884 sec (  17.82 MB/s) | **+72.76%**|   1.797 sec (  58.36 MB/s) |   3.667 sec (  28.59 MB/s) |**+104.12%**|
|   sources.001.2 |   104857600 |   2.916 sec (  35.96 MB/s) |   4.828 sec (  21.72 MB/s) | **+65.55%**|   1.599 sec (  65.58 MB/s) |   3.695 sec (  28.38 MB/s) |**+131.08%**|
|           fib41 |   267914296 |   8.390 sec (  31.93 MB/s) |  37.397 sec (   7.16 MB/s) |**+345.72%**|   3.794 sec (  70.61 MB/s) |  37.158 sec (   7.21 MB/s) |**+879.30%**|
|           rs.13 |   216747218 |   6.425 sec (  33.74 MB/s) |  28.994 sec (   7.48 MB/s) |**+351.30%**|   2.807 sec (  77.21 MB/s) |  31.609 sec (   6.86 MB/s) |**+1025.97%**|
|            tm29 |   268435456 |   8.931 sec (  30.06 MB/s) |  31.267 sec (   8.59 MB/s) |**+250.09%**|   3.729 sec (  71.99 MB/s) |  31.384 sec (   8.55 MB/s) |**+741.69%**|



## Large pages and multi-core systems support

Large-pages and OpenMP improves the libsais performance. Here is an example of such improvements on Manzini Corpus.

| file            |      size | baseline| LP      | LP w 2c | LP w 3c | LP w 4c | LP w 5c | LP w 6c | LP w 7c | LP w 8c |
|:---------------:|:---------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
|       chr22.dna |  34553758 |43.50MB/s|50.18MB/s|61.20MB/s|73.66MB/s|78.91MB/s|81.20MB/s|81.49MB/s|81.52MB/s|80.42MB/s|
|         etext99 | 105277340 |32.96MB/s|40.73MB/s|50.19MB/s|59.34MB/s|62.97MB/s|64.06MB/s|62.83MB/s|63.08MB/s|62.49MB/s|
|     gcc-3.0.tar |  86630400 |44.32MB/s|50.13MB/s|58.51MB/s|68.85MB/s|73.82MB/s|75.76MB/s|76.14MB/s|75.85MB/s|75.24MB/s|
|           howto |  39422105 |42.78MB/s|48.10MB/s|57.38MB/s|67.75MB/s|71.91MB/s|73.67MB/s|73.61MB/s|73.17MB/s|72.38MB/s|
|          jdk13c |  69728899 |42.70MB/s|47.77MB/s|54.50MB/s|64.85MB/s|69.63MB/s|71.66MB/s|72.15MB/s|71.96MB/s|71.24MB/s|
| linux-2.4.5.tar | 116254720 |42.46MB/s|48.85MB/s|57.60MB/s|67.92MB/s|72.29MB/s|73.88MB/s|74.11MB/s|73.59MB/s|73.27MB/s|
|        rctail96 | 114711151 |36.39MB/s|43.19MB/s|50.96MB/s|60.60MB/s|64.33MB/s|65.43MB/s|65.79MB/s|65.78MB/s|65.18MB/s|
|             rfc | 116421901 |39.81MB/s|46.76MB/s|55.92MB/s|66.48MB/s|70.79MB/s|71.68MB/s|72.21MB/s|71.92MB/s|71.06MB/s|
|     sprot34.dat | 109617186 |36.09MB/s|45.06MB/s|53.26MB/s|61.60MB/s|59.69MB/s|62.25MB/s|67.20MB/s|66.84MB/s|66.38MB/s|
|            w3c2 | 104201579 |42.97MB/s|47.09MB/s|54.01MB/s|63.79MB/s|67.67MB/s|69.84MB/s|69.94MB/s|69.65MB/s|68.86MB/s|

Note, multi-core scalability is limited by RAM bandwidth and adding more RAM channels improves performance:
![enwik9 BWT throughput in MB/s on Azure DS14 v2 (Intel Xeon Platinum 8171M)](Azure_enwik9_benchmark.png?raw=true "enwik9 BWT throughput in MB/s on Azure DS14 v2 (Intel Xeon Platinum 8171M)")

## libsais64 for inputs larger than 2GB

Starting from version 2.2.0 libsais64 could process inputs larger than 2GB.

The times below are the minimum of five runs measuring **multi-threaded (MT)** performance of suffix array construction on Azure DS14 v2 (Intel Xeon Platinum 8171M).

|  file           |    size     |    libsais64 2.2.0  (MT)   |   divsufsort64 2.0.2 (MT)  |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|
|         english |  2210395553 |  61.499 sec (  34.28 MB/s) | 435.199 sec (   4.84 MB/s) |**+607.65%**|
|   GRCh38.p13.fa |  3321586957 |  84.068 sec (  37.68 MB/s) | 782.938 sec (   4.05 MB/s) |**+831.32%**|
|         enwik10 | 10000000000 | 303.542 sec (  31.42 MB/s) |1927.351 sec (   4.95 MB/s) |**+534.95%**|

## Additional memory

The libsais reuses space allocated for suffix array during construction. Sometimes this free space is not sufficient for most optimal algorithm (this is uncommon) and libsais will need to fallback to less efficient one (libsais has 4 algorithms at different break-points point: 6k, 4k, 2k and 1k; where k is alphabet size). To improve performance for those cases you could allocating additional space at the end of suffix array.

|  file           |    size     |     libsais + O(n)  (ST)   |     libsais + O(1) (ST)    |speedup (ST)|    libsais + O(n)  (MT)    |     libsais + O(1) (ST)    |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            osdb |    10085684 |   0.222 sec (  45.52 MB/s) |   0.228 sec (  44.20 MB/s) |  **+2.97%**|   0.150 sec (  67.30 MB/s) |   0.162 sec (  62.25 MB/s) |  **+8.11%**|
|           x-ray |     8474240 |   0.190 sec (  44.52 MB/s) |   0.217 sec (  39.11 MB/s) | **+13.82%**|   0.122 sec (  69.46 MB/s) |   0.156 sec (  54.16 MB/s) | **+28.25%**|
|             sao |     7251944 |   0.175 sec (  41.48 MB/s) |   0.182 sec (  39.75 MB/s) |  **+4.37%**|   0.127 sec (  57.26 MB/s) |   0.140 sec (  51.87 MB/s) | **+10.39%**|
|         ooffice |     6152192 |   0.113 sec (  54.55 MB/s) |   0.117 sec (  52.45 MB/s) |  **+4.01%**|   0.081 sec (  76.38 MB/s) |   0.088 sec (  70.30 MB/s) |  **+8.65%**|
|            abac |      200000 |   0.002 sec (  84.36 MB/s) |   0.003 sec (  73.63 MB/s) | **+14.56%**|   0.002 sec ( 105.08 MB/s) |   0.002 sec (  86.64 MB/s) | **+21.27%**|
|           test3 |     2097088 |   0.034 sec (  61.54 MB/s) |   0.037 sec (  56.45 MB/s) |  **+9.03%**|   0.028 sec (  75.76 MB/s) |   0.032 sec (  64.93 MB/s) | **+16.68%**|

> * All other files from [Benchmarks](#benchmarks) above do not suffer from this fallbacks.
