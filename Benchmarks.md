# Specifications
  * OS: Windows 11 Pro 25H2 (64-bit)
  * CPU: AMD Ryzen 9 9950X3D (16C / 32T, 128MB L3 cache, PBO +200Mhz, CO -10 all-core)
  * RAM: 2 x 48 GB DDR5-6000 (28-36-36-60)
  * Compiler: clang 19.1.5 '-O3 -march=znver5 -fopenmp -DLIBSAIS_OPENMP -DNDEBUG'

The times reflect the median of 15 runs for suffix array construction, reported in both **single-threaded (ST)** and **multi-threaded (MT)** modes. Both libraries use 16 threads in MT mode with default OpenMP settings; for each file all ST runs come first, and every run is preceded by a 250 ms pause.

### [Silesia Corpus](https://www.data-compression.info/Corpora/SilesiaCorpus/index.html) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|         dickens |    10192446 |   0.085 sec ( 120.57 MB/s) |   0.310 sec (  32.84 MB/s) |**+267.13%**|   0.038 sec ( 267.59 MB/s) |   0.219 sec (  46.50 MB/s) |**+475.52%**|
|         mozilla |    51220480 |   0.471 sec ( 108.83 MB/s) |   1.219 sec (  42.00 MB/s) |**+159.10%**|   0.299 sec ( 171.13 MB/s) |   0.786 sec (  65.20 MB/s) |**+162.48%**|
|              mr |     9970564 |   0.072 sec ( 138.03 MB/s) |   0.244 sec (  40.84 MB/s) |**+238.00%**|   0.055 sec ( 180.30 MB/s) |   0.190 sec (  52.35 MB/s) |**+244.40%**|
|             nci |    33553445 |   0.242 sec ( 138.55 MB/s) |   0.695 sec (  48.29 MB/s) |**+186.94%**|   0.135 sec ( 247.96 MB/s) |   0.678 sec (  49.52 MB/s) |**+400.78%**|
|         ooffice |     6152192 |   0.058 sec ( 106.81 MB/s) |   0.131 sec (  46.84 MB/s) |**+128.01%**|   0.042 sec ( 147.05 MB/s) |   0.089 sec (  69.40 MB/s) |**+111.87%**|
|            osdb |    10085684 |   0.094 sec ( 106.82 MB/s) |   0.220 sec (  45.93 MB/s) |**+132.59%**|   0.069 sec ( 146.81 MB/s) |   0.193 sec (  52.21 MB/s) |**+181.19%**|
|         reymont |     6627202 |   0.055 sec ( 121.53 MB/s) |   0.175 sec (  37.96 MB/s) |**+220.18%**|   0.032 sec ( 209.28 MB/s) |   0.146 sec (  45.41 MB/s) |**+360.90%**|
|           samba |    21606400 |   0.177 sec ( 122.25 MB/s) |   0.444 sec (  48.71 MB/s) |**+150.95%**|   0.118 sec ( 183.82 MB/s) |   0.337 sec (  64.07 MB/s) |**+186.91%**|
|             sao |     7251944 |   0.075 sec (  97.16 MB/s) |   0.155 sec (  46.72 MB/s) |**+107.94%**|   0.034 sec ( 211.66 MB/s) |   0.105 sec (  69.24 MB/s) |**+205.70%**|
|         webster |    41458703 |   0.383 sec ( 108.29 MB/s) |   1.345 sec (  30.83 MB/s) |**+251.23%**|   0.220 sec ( 188.51 MB/s) |   0.997 sec (  41.57 MB/s) |**+353.44%**|
|           x-ray |     8474240 |   0.108 sec (  78.65 MB/s) |   0.224 sec (  37.87 MB/s) |**+107.66%**|   0.084 sec ( 100.47 MB/s) |   0.123 sec (  68.83 MB/s) | **+45.97%**|
|             xml |     5345280 |   0.036 sec ( 148.52 MB/s) |   0.096 sec (  55.59 MB/s) |**+167.16%**|   0.023 sec ( 228.26 MB/s) |   0.089 sec (  60.01 MB/s) |**+280.36%**|



### [Large Canterbury Corpus](https://www.data-compression.info/Corpora/CanterburyCorpus/) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       bible.txt |     4047392 |   0.032 sec ( 125.91 MB/s) |   0.106 sec (  38.36 MB/s) |**+228.21%**|   0.014 sec ( 294.15 MB/s) |   0.084 sec (  48.05 MB/s) |**+512.18%**|
|          E.coli |     4638690 |   0.040 sec ( 115.22 MB/s) |   0.157 sec (  29.48 MB/s) |**+290.79%**|   0.015 sec ( 315.22 MB/s) |   0.137 sec (  33.86 MB/s) |**+830.86%**|
|    world192.txt |     2473400 |   0.019 sec ( 130.69 MB/s) |   0.056 sec (  44.07 MB/s) |**+196.52%**|   0.012 sec ( 200.26 MB/s) |   0.042 sec (  59.36 MB/s) |**+237.34%**|



### [Manzini Corpus](https://people.unipmn.it/manzini/lightweight/corpus/) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       chr22.dna |    34553758 |   0.328 sec ( 105.45 MB/s) |   1.268 sec (  27.24 MB/s) |**+287.07%**|   0.133 sec ( 259.58 MB/s) |   1.077 sec (  32.08 MB/s) |**+709.06%**|
|         etext99 |   105277340 |   1.185 sec (  88.84 MB/s) |   4.216 sec (  24.97 MB/s) |**+255.79%**|   0.744 sec ( 141.53 MB/s) |   2.958 sec (  35.60 MB/s) |**+297.61%**|
|     gcc-3.0.tar |    86630400 |   0.796 sec ( 108.90 MB/s) |   2.361 sec (  36.70 MB/s) |**+196.76%**|   0.523 sec ( 165.76 MB/s) |   1.705 sec (  50.82 MB/s) |**+226.16%**|
|           howto |    39422105 |   0.376 sec ( 104.72 MB/s) |   1.254 sec (  31.45 MB/s) |**+233.00%**|   0.216 sec ( 182.62 MB/s) |   0.852 sec (  46.28 MB/s) |**+294.58%**|
|          jdk13c |    69728899 |   0.621 sec ( 112.34 MB/s) |   1.883 sec (  37.04 MB/s) |**+203.30%**|   0.416 sec ( 167.49 MB/s) |   1.574 sec (  44.31 MB/s) |**+278.00%**|
| linux-2.4.5.tar |   116254720 |   1.122 sec ( 103.66 MB/s) |   3.407 sec (  34.13 MB/s) |**+203.75%**|   0.738 sec ( 157.51 MB/s) |   2.390 sec (  48.64 MB/s) |**+223.84%**|
|        rctail96 |   114711151 |   1.170 sec (  98.03 MB/s) |   3.831 sec (  29.94 MB/s) |**+227.39%**|   0.781 sec ( 146.96 MB/s) |   2.924 sec (  39.24 MB/s) |**+274.54%**|
|             rfc |   116421901 |   1.112 sec ( 104.68 MB/s) |   3.651 sec (  31.88 MB/s) |**+228.30%**|   0.717 sec ( 162.31 MB/s) |   2.625 sec (  44.35 MB/s) |**+266.02%**|
|     sprot34.dat |   109617186 |   1.113 sec (  98.47 MB/s) |   3.793 sec (  28.90 MB/s) |**+240.76%**|   0.726 sec ( 150.90 MB/s) |   2.523 sec (  43.45 MB/s) |**+247.33%**|
|            w3c2 |   104201579 |   0.975 sec ( 106.83 MB/s) |   2.884 sec (  36.12 MB/s) |**+195.73%**|   0.691 sec ( 150.72 MB/s) |   2.290 sec (  45.50 MB/s) |**+231.25%**|



### [Large Text Compression Benchmark Corpus](https://www.mattmahoney.net/dc/textdata.html) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|          enwik8 |   100000000 |   1.065 sec (  93.92 MB/s) |   3.841 sec (  26.03 MB/s) |**+260.77%**|   0.589 sec ( 169.76 MB/s) |   2.562 sec (  39.04 MB/s) |**+334.83%**|
|          enwik9 |  1000000000 |  12.227 sec (  81.79 MB/s) |  43.415 sec (  23.03 MB/s) |**+255.08%**|   7.905 sec ( 126.50 MB/s) |  29.591 sec (  33.79 MB/s) |**+274.31%**|



### [The Gauntlet Corpus](https://github.com/michaelmaniscalco/gauntlet_corpus) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            abac |      200000 |   0.001 sec ( 159.53 MB/s) |   0.001 sec ( 169.18 MB/s) |    -5.70%  |   0.001 sec ( 154.94 MB/s) |   0.001 sec ( 137.23 MB/s) | **+12.91%**|
|            abba |    10500596 |   0.059 sec ( 176.59 MB/s) |   0.267 sec (  39.35 MB/s) |**+348.78%**|   0.035 sec ( 299.05 MB/s) |   0.296 sec (  35.50 MB/s) |**+742.37%**|
|        book1x20 |    15375420 |   0.124 sec ( 124.12 MB/s) |   0.484 sec (  31.80 MB/s) |**+290.34%**|   0.078 sec ( 198.05 MB/s) |   0.327 sec (  46.96 MB/s) |**+321.73%**|
|   fib_s14930352 |    14930352 |   0.109 sec ( 136.94 MB/s) |   0.538 sec (  27.77 MB/s) |**+393.19%**|   0.070 sec ( 214.32 MB/s) |   0.593 sec (  25.17 MB/s) |**+751.61%**|
|           fss10 |    12078908 |   0.080 sec ( 151.15 MB/s) |   0.411 sec (  29.37 MB/s) |**+414.65%**|   0.049 sec ( 246.23 MB/s) |   0.465 sec (  25.95 MB/s) |**+848.68%**|
|            fss9 |     2851443 |   0.018 sec ( 154.68 MB/s) |   0.084 sec (  33.79 MB/s) |**+357.78%**|   0.011 sec ( 252.23 MB/s) |   0.104 sec (  27.53 MB/s) |**+816.34%**|
|         houston |     3839141 |   0.015 sec ( 260.64 MB/s) |   0.013 sec ( 290.51 MB/s) |   -10.28%  |   0.007 sec ( 531.25 MB/s) |   0.017 sec ( 221.21 MB/s) |**+140.15%**|
|       paper5x80 |      956322 |   0.006 sec ( 148.20 MB/s) |   0.015 sec (  62.75 MB/s) |**+136.17%**|   0.005 sec ( 195.75 MB/s) |   0.014 sec (  68.40 MB/s) |**+186.17%**|
|           test1 |     2097152 |   0.020 sec ( 105.79 MB/s) |   0.034 sec (  61.60 MB/s) | **+71.75%**|   0.009 sec ( 235.00 MB/s) |   0.036 sec (  58.00 MB/s) |**+305.14%**|
|           test2 |     2097152 |   0.020 sec ( 104.97 MB/s) |   0.026 sec (  81.29 MB/s) | **+29.13%**|   0.009 sec ( 226.48 MB/s) |   0.027 sec (  77.08 MB/s) |**+193.84%**|
|           test3 |     2097088 |   0.016 sec ( 131.64 MB/s) |   0.026 sec (  81.07 MB/s) | **+62.39%**|   0.015 sec ( 140.44 MB/s) |   0.033 sec (  63.29 MB/s) |**+121.90%**|



### [Pizza & Chilli Corpus](https://pizzachili.dcc.uchile.cl/texts.html) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|        dblp.xml |   296135874 |   3.097 sec (  95.64 MB/s) |   9.362 sec (  31.63 MB/s) |**+202.34%**|   2.199 sec ( 134.67 MB/s) |   7.221 sec (  41.01 MB/s) |**+228.40%**|
|             dna |   403927746 |   4.978 sec (  81.14 MB/s) |  19.232 sec (  21.00 MB/s) |**+286.33%**|   2.757 sec ( 146.51 MB/s) |  16.477 sec (  24.52 MB/s) |**+497.65%**|
|  english.1024MB |  1073741824 |  14.517 sec (  73.97 MB/s) |  52.624 sec (  20.40 MB/s) |**+262.51%**|   9.449 sec ( 113.64 MB/s) |  37.603 sec (  28.55 MB/s) |**+297.97%**|
|         pitches |    55832855 |   0.592 sec (  94.35 MB/s) |   1.681 sec (  33.21 MB/s) |**+184.11%**|   0.377 sec ( 148.14 MB/s) |   1.009 sec (  55.32 MB/s) |**+167.77%**|
|        proteins |  1184051855 |  16.767 sec (  70.62 MB/s) |  67.793 sec (  17.47 MB/s) |**+304.33%**|  10.965 sec ( 107.98 MB/s) |  36.295 sec (  32.62 MB/s) |**+231.00%**|
|         sources |   210866607 |   2.185 sec (  96.52 MB/s) |   6.819 sec (  30.92 MB/s) |**+212.11%**|   1.485 sec ( 142.02 MB/s) |   4.766 sec (  44.24 MB/s) |**+221.03%**|



### [Pizza & Chilli Repetitive Corpus](https://pizzachili.dcc.uchile.cl/repcorpus.html) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            cere |   461286644 |   4.951 sec (  93.17 MB/s) |  19.371 sec (  23.81 MB/s) |**+291.25%**|   3.189 sec ( 144.64 MB/s) |  16.574 sec (  27.83 MB/s) |**+419.69%**|
|       coreutils |   205281778 |   1.991 sec ( 103.08 MB/s) |   6.924 sec (  29.65 MB/s) |**+247.70%**|   1.412 sec ( 145.39 MB/s) |   5.246 sec (  39.13 MB/s) |**+271.54%**|
| einstein.de.txt |    92758441 |   0.843 sec ( 110.05 MB/s) |   2.705 sec (  34.29 MB/s) |**+220.99%**|   0.546 sec ( 169.98 MB/s) |   2.152 sec (  43.10 MB/s) |**+294.35%**|
| einstein.en.txt |   467626544 |   5.087 sec (  91.93 MB/s) |  15.996 sec (  29.23 MB/s) |**+214.46%**|   3.377 sec ( 138.47 MB/s) |  12.687 sec (  36.86 MB/s) |**+275.69%**|
|Escherichia_Coli |   112689515 |   1.152 sec (  97.84 MB/s) |   4.617 sec (  24.41 MB/s) |**+300.81%**|   0.724 sec ( 155.63 MB/s) |   3.833 sec (  29.40 MB/s) |**+429.35%**|
|       influenza |   154808555 |   1.442 sec ( 107.39 MB/s) |   5.396 sec (  28.69 MB/s) |**+274.35%**|   0.986 sec ( 156.99 MB/s) |   4.620 sec (  33.51 MB/s) |**+368.45%**|
|          kernel |   257961616 |   2.526 sec ( 102.14 MB/s) |   9.089 sec (  28.38 MB/s) |**+259.86%**|   1.750 sec ( 147.42 MB/s) |   6.640 sec (  38.85 MB/s) |**+279.47%**|
|            para |   429265758 |   4.738 sec (  90.61 MB/s) |  18.522 sec (  23.18 MB/s) |**+290.97%**|   3.066 sec ( 140.03 MB/s) |  15.780 sec (  27.20 MB/s) |**+414.74%**|
|   world_leaders |    46968181 |   0.286 sec ( 164.42 MB/s) |   0.777 sec (  60.45 MB/s) |**+171.98%**|   0.151 sec ( 311.67 MB/s) |   0.647 sec (  72.63 MB/s) |**+329.14%**|
|dblp.xml.00001.1 |   104857600 |   1.833 sec (  57.19 MB/s) |   3.849 sec (  27.24 MB/s) |**+109.94%**|   0.690 sec ( 152.02 MB/s) |   3.305 sec (  31.73 MB/s) |**+379.14%**|
|dblp.xml.00001.2 |   104857600 |   1.838 sec (  57.04 MB/s) |   3.900 sec (  26.88 MB/s) |**+112.16%**|   0.689 sec ( 152.12 MB/s) |   3.361 sec (  31.20 MB/s) |**+387.63%**|
| dblp.xml.0001.1 |   104857600 |   1.822 sec (  57.57 MB/s) |   3.841 sec (  27.30 MB/s) |**+110.89%**|   0.689 sec ( 152.18 MB/s) |   3.294 sec (  31.83 MB/s) |**+378.05%**|
| dblp.xml.0001.2 |   104857600 |   1.824 sec (  57.48 MB/s) |   3.870 sec (  27.10 MB/s) |**+112.15%**|   0.699 sec ( 150.01 MB/s) |   3.304 sec (  31.74 MB/s) |**+372.63%**|
|       dna.001.1 |   104857600 |   1.916 sec (  54.73 MB/s) |   4.579 sec (  22.90 MB/s) |**+139.00%**|   0.655 sec ( 160.12 MB/s) |   3.879 sec (  27.03 MB/s) |**+492.39%**|
|   english.001.2 |   104857600 |   1.993 sec (  52.62 MB/s) |   4.478 sec (  23.42 MB/s) |**+124.72%**|   0.766 sec ( 136.89 MB/s) |   3.297 sec (  31.81 MB/s) |**+330.36%**|
|  proteins.001.1 |   104857600 |   1.963 sec (  53.42 MB/s) |   4.591 sec (  22.84 MB/s) |**+133.90%**|   0.740 sec ( 141.76 MB/s) |   2.923 sec (  35.87 MB/s) |**+295.16%**|
|   sources.001.2 |   104857600 |   1.762 sec (  59.50 MB/s) |   3.769 sec (  27.82 MB/s) |**+113.85%**|   0.711 sec ( 147.43 MB/s) |   3.009 sec (  34.85 MB/s) |**+323.04%**|
|           fib41 |   267914296 |   3.083 sec (  86.89 MB/s) |  16.813 sec (  15.93 MB/s) |**+445.31%**|   1.799 sec ( 148.95 MB/s) |  16.923 sec (  15.83 MB/s) |**+840.86%**|
|           rs.13 |   216747218 |   2.129 sec ( 101.78 MB/s) |  13.135 sec (  16.50 MB/s) |**+516.81%**|   1.457 sec ( 148.74 MB/s) |  13.262 sec (  16.34 MB/s) |**+810.11%**|
|            tm29 |   268435456 |   8.500 sec (  31.58 MB/s) |  19.786 sec (  13.57 MB/s) |**+132.79%**|   1.986 sec ( 135.17 MB/s) |  20.266 sec (  13.25 MB/s) |**+920.48%**|



### [Skyline Corpus](http://panthema.net/2012/1119-eSAIS-Inducing-Suffix-and-LCP-Arrays-in-External-Memory/eSAIS-DC3-LCP-0.5.0/src/input/skyline.h.html) ###

|  file           |    size     |     libsais 2.12.0 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.12.0 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|   skyline20.txt |     1048575 |   0.013 sec (  82.76 MB/s) |   0.051 sec (  20.59 MB/s) |**+301.83%**|   0.009 sec ( 118.47 MB/s) |   0.054 sec (  19.40 MB/s) |**+510.71%**|
|   skyline22.txt |     4194303 |   0.060 sec (  69.54 MB/s) |   0.262 sec (  15.99 MB/s) |**+334.88%**|   0.032 sec ( 131.99 MB/s) |   0.273 sec (  15.36 MB/s) |**+759.21%**|
|   skyline24.txt |    16777215 |   0.349 sec (  48.03 MB/s) |   1.369 sec (  12.25 MB/s) |**+292.05%**|   0.129 sec ( 129.99 MB/s) |   1.378 sec (  12.17 MB/s) |**+967.87%**|
|   skyline26.txt |    67108863 |   2.804 sec (  23.94 MB/s) |   7.611 sec (   8.82 MB/s) |**+171.46%**|   0.702 sec (  95.55 MB/s) |   7.608 sec (   8.82 MB/s) |**+983.29%**|
|   skyline28.txt |   268435455 |  12.350 sec (  21.74 MB/s) |  37.239 sec (   7.21 MB/s) |**+201.53%**|   3.382 sec (  79.38 MB/s) |  37.102 sec (   7.24 MB/s) |**+997.09%**|

## Additional memory

The libsais reuses the space allocated for the suffix array during construction. In rare cases, this space is not large enough for the fastest algorithm and libsais will need to fallback to less efficient one (libsais has 4 algorithms at different break-points point: 6k, 4k, 2k and 1k; where k is alphabet size). To improve performance for those cases you could allocating additional space at the end of suffix array.

|  file           |    size     |     libsais + O(n)  (ST)   |     libsais + O(1) (ST)    |speedup (ST)|    libsais + O(n)  (MT)    |     libsais + O(1) (MT)    |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            osdb |    10085684 |   0.086 sec ( 116.90 MB/s) |   0.094 sec ( 107.55 MB/s) |  **+8.70%**|   0.062 sec ( 162.77 MB/s) |   0.068 sec ( 147.58 MB/s) | **+10.30%**|
|           x-ray |     8474240 |   0.081 sec ( 105.25 MB/s) |   0.107 sec (  79.36 MB/s) | **+32.62%**|   0.061 sec ( 138.99 MB/s) |   0.084 sec ( 100.63 MB/s) | **+38.13%**|
|         ooffice |     6152192 |   0.053 sec ( 117.14 MB/s) |   0.058 sec ( 106.83 MB/s) |  **+9.66%**|   0.037 sec ( 166.66 MB/s) |   0.045 sec ( 136.95 MB/s) | **+21.69%**|
|           test3 |     2097088 |   0.015 sec ( 138.95 MB/s) |   0.016 sec ( 131.89 MB/s) |  **+5.35%**|   0.014 sec ( 149.02 MB/s) |   0.015 sec ( 141.20 MB/s) |  **+5.54%**|

> * All other files from [Benchmarks](#benchmarks) above do not suffer from this fallbacks.
