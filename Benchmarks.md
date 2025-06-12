# Specifications
  * OS: Windows 10 Pro (64-bit)
  * CPU: Intel Core i7-9700K (12 MB cache, 5.0 GHz on all cores)
  * RAM: 4 x 16 GB DDR4-3600 MHz (14-14-14-30-300-2T with tight subtimings)
  * Compiler: clang 19.1.5 '-O3 -march=skylake -fopenmp -DLIBSAIS_OPENMP -DNDEBUG'

The times reflect the median of 11 runs for suffix-array construction, reported in both **single-threaded (ST)** and **multi-threaded (MT)** modes.

### [Silesia Corpus](https://www.data-compression.info/Corpora/SilesiaCorpus/index.html) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|         dickens |    10192446 |   0.191 sec (  53.50 MB/s) |   0.452 sec (  22.56 MB/s) |**+137.16%**|   0.101 sec ( 100.64 MB/s) |   0.303 sec (  33.59 MB/s) |**+199.59%**|
|         mozilla |    51220480 |   1.053 sec (  48.66 MB/s) |   1.829 sec (  28.00 MB/s) | **+73.80%**|   0.653 sec (  78.47 MB/s) |   1.228 sec (  41.70 MB/s) | **+88.18%**|
|              mr |     9970564 |   0.172 sec (  58.13 MB/s) |   0.395 sec (  25.27 MB/s) |**+130.06%**|   0.136 sec (  73.45 MB/s) |   0.290 sec (  34.39 MB/s) |**+113.55%**|
|             nci |    33553445 |   0.635 sec (  52.83 MB/s) |   1.205 sec (  27.83 MB/s) | **+89.82%**|   0.403 sec (  83.27 MB/s) |   1.094 sec (  30.66 MB/s) |**+171.59%**|
|         ooffice |     6152192 |   0.108 sec (  56.92 MB/s) |   0.184 sec (  33.50 MB/s) | **+69.88%**|   0.085 sec (  72.11 MB/s) |   0.120 sec (  51.13 MB/s) | **+41.04%**|
|            osdb |    10085684 |   0.214 sec (  47.18 MB/s) |   0.378 sec (  26.70 MB/s) | **+76.74%**|   0.155 sec (  64.94 MB/s) |   0.296 sec (  34.02 MB/s) | **+90.90%**|
|         reymont |     6627202 |   0.115 sec (  57.71 MB/s) |   0.259 sec (  25.59 MB/s) |**+125.53%**|   0.078 sec (  84.73 MB/s) |   0.196 sec (  33.78 MB/s) |**+150.85%**|
|           samba |    21606400 |   0.378 sec (  57.10 MB/s) |   0.680 sec (  31.75 MB/s) | **+79.82%**|   0.281 sec (  76.83 MB/s) |   0.502 sec (  43.02 MB/s) | **+78.61%**|
|             sao |     7251944 |   0.176 sec (  41.28 MB/s) |   0.238 sec (  30.51 MB/s) | **+35.32%**|   0.141 sec (  51.45 MB/s) |   0.159 sec (  45.66 MB/s) | **+12.67%**|
|         webster |    41458703 |   0.993 sec (  41.75 MB/s) |   2.221 sec (  18.67 MB/s) |**+123.67%**|   0.555 sec (  74.66 MB/s) |   1.697 sec (  24.44 MB/s) |**+205.52%**|
|           x-ray |     8474240 |   0.203 sec (  41.84 MB/s) |   0.335 sec (  25.30 MB/s) | **+65.38%**|   0.134 sec (  63.26 MB/s) |   0.181 sec (  46.82 MB/s) | **+35.13%**|
|             xml |     5345280 |   0.070 sec (  75.84 MB/s) |   0.139 sec (  38.40 MB/s) | **+97.51%**|   0.046 sec ( 116.44 MB/s) |   0.111 sec (  48.14 MB/s) |**+141.89%**|



### [Large Canterbury Corpus](https://www.data-compression.info/Corpora/CanterburyCorpus/) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       bible.txt |     4047392 |   0.060 sec (  66.93 MB/s) |   0.143 sec (  28.37 MB/s) |**+135.95%**|   0.035 sec ( 115.09 MB/s) |   0.103 sec (  39.42 MB/s) |**+191.95%**|
|          E.coli |     4638690 |   0.074 sec (  62.70 MB/s) |   0.206 sec (  22.50 MB/s) |**+178.71%**|   0.036 sec ( 128.18 MB/s) |   0.161 sec (  28.86 MB/s) |**+344.18%**|
|    world192.txt |     2473400 |   0.034 sec (  73.36 MB/s) |   0.076 sec (  32.40 MB/s) |**+126.43%**|   0.022 sec ( 111.69 MB/s) |   0.051 sec (  48.46 MB/s) |**+130.47%**|



### [Manzini Corpus](https://people.unipmn.it/manzini/lightweight/corpus/) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       chr22.dna |    34553758 |   0.744 sec (  46.42 MB/s) |   2.047 sec (  16.88 MB/s) |**+175.02%**|   0.401 sec (  86.22 MB/s) |   1.720 sec (  20.09 MB/s) |**+329.15%**|
|         etext99 |   105277340 |   3.201 sec (  32.89 MB/s) |   7.217 sec (  14.59 MB/s) |**+125.51%**|   1.837 sec (  57.32 MB/s) |   5.444 sec (  19.34 MB/s) |**+196.40%**|
|     gcc-3.0.tar |    86630400 |   1.872 sec (  46.28 MB/s) |   3.560 sec (  24.34 MB/s) | **+90.16%**|   1.269 sec (  68.25 MB/s) |   2.596 sec (  33.37 MB/s) |**+104.52%**|
|           howto |    39422105 |   0.894 sec (  44.08 MB/s) |   1.934 sec (  20.38 MB/s) |**+116.27%**|   0.597 sec (  66.08 MB/s) |   1.342 sec (  29.37 MB/s) |**+125.02%**|
|          jdk13c |    69728899 |   1.577 sec (  44.20 MB/s) |   3.088 sec (  22.58 MB/s) | **+95.79%**|   1.070 sec (  65.16 MB/s) |   2.600 sec (  26.82 MB/s) |**+142.97%**|
| linux-2.4.5.tar |   116254720 |   2.629 sec (  44.21 MB/s) |   5.164 sec (  22.51 MB/s) | **+96.40%**|   1.745 sec (  66.61 MB/s) |   3.670 sec (  31.67 MB/s) |**+110.29%**|
|        rctail96 |   114711151 |   3.039 sec (  37.74 MB/s) |   6.413 sec (  17.89 MB/s) |**+111.00%**|   1.697 sec (  67.61 MB/s) |   5.130 sec (  22.36 MB/s) |**+202.34%**|
|             rfc |   116421901 |   2.813 sec (  41.39 MB/s) |   5.693 sec (  20.45 MB/s) |**+102.42%**|   1.574 sec (  73.95 MB/s) |   4.261 sec (  27.32 MB/s) |**+170.66%**|
|     sprot34.dat |   109617186 |   2.926 sec (  37.46 MB/s) |   6.357 sec (  17.24 MB/s) |**+117.25%**|   1.590 sec (  68.92 MB/s) |   4.549 sec (  24.09 MB/s) |**+186.06%**|
|            w3c2 |   104201579 |   2.337 sec (  44.59 MB/s) |   4.500 sec (  23.15 MB/s) | **+92.58%**|   1.449 sec (  71.89 MB/s) |   3.635 sec (  28.67 MB/s) |**+150.77%**|



### [Large Text Compression Benchmark Corpus](https://www.mattmahoney.net/dc/textdata.html) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|          enwik8 |   100000000 |   2.812 sec (  35.56 MB/s) |   6.238 sec (  16.03 MB/s) |**+121.80%**|   1.500 sec (  66.67 MB/s) |   4.437 sec (  22.54 MB/s) |**+195.82%**|
|          enwik9 |  1000000000 |  40.587 sec (  24.64 MB/s) |  82.544 sec (  12.11 MB/s) |**+103.38%**|  18.280 sec (  54.70 MB/s) |  63.639 sec (  15.71 MB/s) |**+248.13%**|



### [The Gauntlet Corpus](https://github.com/michaelmaniscalco/gauntlet_corpus) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            abac |      200000 |   0.002 sec (  94.13 MB/s) |   0.002 sec ( 125.20 MB/s) |   -24.81%  |   0.002 sec ( 124.73 MB/s) |   0.002 sec ( 121.79 MB/s) |  **+2.42%**|
|            abba |    10500596 |   0.148 sec (  71.12 MB/s) |   0.444 sec (  23.66 MB/s) |**+200.60%**|   0.078 sec ( 134.74 MB/s) |   0.462 sec (  22.71 MB/s) |**+493.23%**|
|        book1x20 |    15375420 |   0.307 sec (  50.14 MB/s) |   0.714 sec (  21.54 MB/s) |**+132.71%**|   0.171 sec (  89.74 MB/s) |   0.504 sec (  30.51 MB/s) |**+194.12%**|
|   fib_s14930352 |    14930352 |   0.324 sec (  46.04 MB/s) |   1.095 sec (  13.63 MB/s) |**+237.83%**|   0.162 sec (  92.11 MB/s) |   1.095 sec (  13.63 MB/s) |**+575.66%**|
|           fss10 |    12078908 |   0.217 sec (  55.58 MB/s) |   0.834 sec (  14.48 MB/s) |**+283.93%**|   0.120 sec ( 101.03 MB/s) |   0.827 sec (  14.60 MB/s) |**+592.10%**|
|            fss9 |     2851443 |   0.036 sec (  79.09 MB/s) |   0.129 sec (  22.12 MB/s) |**+257.48%**|   0.022 sec ( 128.85 MB/s) |   0.129 sec (  22.11 MB/s) |**+482.74%**|
|         houston |     3839141 |   0.033 sec ( 118.12 MB/s) |   0.024 sec ( 162.30 MB/s) |   -27.22%  |   0.017 sec ( 219.93 MB/s) |   0.023 sec ( 165.40 MB/s) | **+32.97%**|
|       paper5x80 |      956322 |   0.012 sec (  79.14 MB/s) |   0.024 sec (  39.34 MB/s) |**+101.17%**|   0.008 sec ( 116.66 MB/s) |   0.020 sec (  47.66 MB/s) |**+144.77%**|
|           test1 |     2097152 |   0.036 sec (  57.74 MB/s) |   0.083 sec (  25.28 MB/s) |**+128.44%**|   0.019 sec ( 108.61 MB/s) |   0.077 sec (  27.23 MB/s) |**+298.88%**|
|           test2 |     2097152 |   0.036 sec (  58.60 MB/s) |   0.055 sec (  38.21 MB/s) | **+53.36%**|   0.019 sec ( 109.98 MB/s) |   0.052 sec (  40.69 MB/s) |**+170.28%**|
|           test3 |     2097088 |   0.035 sec (  60.63 MB/s) |   0.044 sec (  47.71 MB/s) | **+27.07%**|   0.029 sec (  71.76 MB/s) |   0.045 sec (  46.51 MB/s) | **+54.31%**|



### [Pizza & Chilli Corpus](https://pizzachili.dcc.uchile.cl/texts.html) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|        dblp.xml |   296135874 |   8.311 sec (  35.63 MB/s) |  16.167 sec (  18.32 MB/s) | **+94.52%**|   4.486 sec (  66.01 MB/s) |  13.193 sec (  22.45 MB/s) |**+194.08%**|
|             dna |   403927746 |  14.293 sec (  28.26 MB/s) |  36.903 sec (  10.95 MB/s) |**+158.19%**|   6.159 sec (  65.58 MB/s) |  32.486 sec (  12.43 MB/s) |**+427.42%**|
|  english.1024MB |  1073741824 |  50.606 sec (  21.22 MB/s) | 104.811 sec (  10.24 MB/s) |**+107.11%**|  21.722 sec (  49.43 MB/s) |  85.056 sec (  12.62 MB/s) |**+291.57%**|
|         pitches |    55832855 |   1.154 sec (  48.40 MB/s) |   2.416 sec (  23.11 MB/s) |**+109.43%**|   0.740 sec (  75.41 MB/s) |   1.450 sec (  38.50 MB/s) | **+95.86%**|
|        proteins |  1184051855 |  59.536 sec (  19.89 MB/s) | 124.421 sec (   9.52 MB/s) |**+108.98%**|  24.944 sec (  47.47 MB/s) |  80.608 sec (  14.69 MB/s) |**+223.16%**|
|         sources |   210866607 |   5.212 sec (  40.46 MB/s) |  10.395 sec (  20.29 MB/s) | **+99.46%**|   3.018 sec (  69.88 MB/s) |   7.483 sec (  28.18 MB/s) |**+147.98%**|



### [Pizza & Chilli Repetitive Corpus](https://pizzachili.dcc.uchile.cl/repcorpus.html) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            cere |   461286644 |  16.361 sec (  28.19 MB/s) |  35.093 sec (  13.14 MB/s) |**+114.50%**|   6.935 sec (  66.52 MB/s) |  31.072 sec (  14.85 MB/s) |**+348.08%**|
|       coreutils |   205281778 |   5.047 sec (  40.67 MB/s) |  10.652 sec (  19.27 MB/s) |**+111.04%**|   2.862 sec (  71.73 MB/s) |   8.372 sec (  24.52 MB/s) |**+192.55%**|
| einstein.de.txt |    92758441 |   2.386 sec (  38.88 MB/s) |   4.628 sec (  20.04 MB/s) | **+93.98%**|   1.272 sec (  72.91 MB/s) |   3.782 sec (  24.53 MB/s) |**+197.24%**|
| einstein.en.txt |   467626544 |  14.944 sec (  31.29 MB/s) |  32.041 sec (  14.59 MB/s) |**+114.40%**|   7.439 sec (  62.86 MB/s) |  27.360 sec (  17.09 MB/s) |**+267.80%**|
|Escherichia_Coli |   112689515 |   3.289 sec (  34.26 MB/s) |   7.607 sec (  14.81 MB/s) |**+131.25%**|   1.563 sec (  72.08 MB/s) |   6.424 sec (  17.54 MB/s) |**+310.91%**|
|       influenza |   154808555 |   3.702 sec (  41.82 MB/s) |   9.195 sec (  16.84 MB/s) |**+148.38%**|   2.049 sec (  75.57 MB/s) |   7.957 sec (  19.46 MB/s) |**+288.40%**|
|          kernel |   257961616 |   6.782 sec (  38.04 MB/s) |  13.568 sec (  19.01 MB/s) |**+100.06%**|   3.617 sec (  71.32 MB/s) |  10.282 sec (  25.09 MB/s) |**+184.27%**|
|            para |   429265758 |  15.425 sec (  27.83 MB/s) |  33.343 sec (  12.87 MB/s) |**+116.16%**|   6.518 sec (  65.85 MB/s) |  29.341 sec (  14.63 MB/s) |**+350.12%**|
|   world_leaders |    46968181 |   0.825 sec (  56.91 MB/s) |   1.307 sec (  35.93 MB/s) | **+58.42%**|   0.421 sec ( 111.45 MB/s) |   1.063 sec (  44.18 MB/s) |**+152.25%**|
|dblp.xml.00001.1 |   104857600 |   2.961 sec (  35.42 MB/s) |   5.682 sec (  18.45 MB/s) | **+91.94%**|   1.506 sec (  69.63 MB/s) |   4.897 sec (  21.41 MB/s) |**+225.16%**|
|dblp.xml.00001.2 |   104857600 |   2.982 sec (  35.17 MB/s) |   5.748 sec (  18.24 MB/s) | **+92.79%**|   1.507 sec (  69.57 MB/s) |   4.944 sec (  21.21 MB/s) |**+228.03%**|
| dblp.xml.0001.1 |   104857600 |   3.006 sec (  34.88 MB/s) |   5.678 sec (  18.47 MB/s) | **+88.88%**|   1.512 sec (  69.36 MB/s) |   4.898 sec (  21.41 MB/s) |**+223.98%**|
| dblp.xml.0001.2 |   104857600 |   2.996 sec (  35.00 MB/s) |   5.696 sec (  18.41 MB/s) | **+90.10%**|   1.510 sec (  69.45 MB/s) |   4.846 sec (  21.64 MB/s) |**+220.94%**|
|       dna.001.1 |   104857600 |   2.955 sec (  35.49 MB/s) |   6.481 sec (  16.18 MB/s) |**+119.33%**|   1.369 sec (  76.57 MB/s) |   5.446 sec (  19.25 MB/s) |**+297.68%**|
|   english.001.2 |   104857600 |   3.285 sec (  31.92 MB/s) |   5.959 sec (  17.60 MB/s) | **+81.40%**|   1.579 sec (  66.39 MB/s) |   4.350 sec (  24.10 MB/s) |**+175.45%**|
|  proteins.001.1 |   104857600 |   3.245 sec (  32.31 MB/s) |   5.943 sec (  17.65 MB/s) | **+83.11%**|   1.617 sec (  64.85 MB/s) |   3.724 sec (  28.15 MB/s) |**+130.33%**|
|   sources.001.2 |   104857600 |   2.815 sec (  37.25 MB/s) |   4.917 sec (  21.32 MB/s) | **+74.69%**|   1.465 sec (  71.59 MB/s) |   3.783 sec (  27.72 MB/s) |**+158.30%**|
|           fib41 |   267914296 |   8.610 sec (  31.12 MB/s) |  38.252 sec (   7.00 MB/s) |**+344.28%**|   3.497 sec (  76.62 MB/s) |  37.812 sec (   7.09 MB/s) |**+981.35%**|
|           rs.13 |   216747218 |   6.559 sec (  33.04 MB/s) |  29.787 sec (   7.28 MB/s) |**+354.13%**|   2.710 sec (  79.99 MB/s) |  31.887 sec (   6.80 MB/s) |**+1076.78%**|
|            tm29 |   268435456 |   8.104 sec (  33.12 MB/s) |  32.866 sec (   8.17 MB/s) |**+305.55%**|   3.171 sec (  84.65 MB/s) |  33.149 sec (   8.10 MB/s) |**+945.39%**|



### [Skyline Corpus](http://panthema.net/2012/1119-eSAIS-Inducing-Suffix-and-LCP-Arrays-in-External-Memory/eSAIS-DC3-LCP-0.5.0/src/input/skyline.h.html) ###

|  file           |    size     |     libsais 2.10.2 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.2 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|   skyline20.txt |     1048575 |   0.022 sec (  46.76 MB/s) |   0.065 sec (  16.11 MB/s) |**+190.35%**|   0.014 sec (  77.26 MB/s) |   0.064 sec (  16.32 MB/s) |**+373.42%**|
|   skyline22.txt |     4194303 |   0.103 sec (  40.92 MB/s) |   0.310 sec (  13.51 MB/s) |**+202.87%**|   0.053 sec (  78.80 MB/s) |   0.307 sec (  13.66 MB/s) |**+476.65%**|
|   skyline24.txt |    16777215 |   0.544 sec (  30.87 MB/s) |   1.778 sec (   9.44 MB/s) |**+227.11%**|   0.272 sec (  61.69 MB/s) |   1.753 sec (   9.57 MB/s) |**+544.68%**|
|   skyline26.txt |    67108863 |   2.469 sec (  27.18 MB/s) |   9.129 sec (   7.35 MB/s) |**+269.71%**|   1.147 sec (  58.52 MB/s) |   9.002 sec (   7.46 MB/s) |**+684.96%**|
|   skyline28.txt |   268435455 |  12.014 sec (  22.34 MB/s) |  47.934 sec (   5.60 MB/s) |**+298.98%**|   4.815 sec (  55.75 MB/s) |  47.442 sec (   5.66 MB/s) |**+885.34%**|

## Large pages and multi-core systems support

Using large pages and OpenMP multithreading boosts the libsais performance.  Below is an example of the speedup achieved on the Manzini corpus:

| file            |   size    | 1 core  | 1c + LP | 2c + LP | 3c + LP | 4c + LP | 5c + LP | 6c + LP | 7c + LP | 8c + LP |
|:---------------:|:---------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
|       chr22.dna |  34553758 |45.33MB/s|54.64MB/s|65.41MB/s|79.72MB/s|85.76MB/s|88.05MB/s|89.53MB/s|88.90MB/s|88.84MB/s|
|         etext99 | 105277340 |32.90MB/s|42.24MB/s|51.63MB/s|60.49MB/s|64.86MB/s|66.15MB/s|66.69MB/s|66.82MB/s|66.20MB/s|
|     gcc-3.0.tar |  86630400 |45.68MB/s|53.35MB/s|60.52MB/s|70.95MB/s|76.31MB/s|78.74MB/s|79.50MB/s|79.79MB/s|79.30MB/s|
|           howto |  39422105 |43.64MB/s|51.01MB/s|61.01MB/s|70.94MB/s|75.75MB/s|78.36MB/s|78.86MB/s|78.61MB/s|78.14MB/s|
|          jdk13c |  69728899 |43.97MB/s|50.99MB/s|56.30MB/s|66.75MB/s|72.18MB/s|74.07MB/s|75.27MB/s|75.81MB/s|75.41MB/s|
| linux-2.4.5.tar | 116254720 |43.40MB/s|52.17MB/s|59.87MB/s|70.44MB/s|75.25MB/s|77.20MB/s|78.02MB/s|78.07MB/s|77.83MB/s|
|        rctail96 | 114711151 |36.68MB/s|44.93MB/s|52.76MB/s|62.16MB/s|66.52MB/s|68.58MB/s|68.99MB/s|69.00MB/s|68.44MB/s|
|             rfc | 116421901 |39.97MB/s|48.86MB/s|57.77MB/s|68.04MB/s|73.59MB/s|75.18MB/s|75.48MB/s|75.75MB/s|75.45MB/s|
|     sprot34.dat | 109617186 |36.88MB/s|46.68MB/s|55.01MB/s|64.42MB/s|69.12MB/s|70.42MB/s|71.22MB/s|71.12MB/s|70.59MB/s|
|            w3c2 | 104201579 |43.50MB/s|49.50MB/s|55.46MB/s|65.22MB/s|70.28MB/s|72.56MB/s|73.15MB/s|73.30MB/s|72.92MB/s|

## Additional memory

The libsais reuses the space allocated for the suffix array during construction. In rare cases, this space is not large enough for the fastest algorithm and libsais will need to fallback to less efficient one (libsais has 4 algorithms at different break-points point: 6k, 4k, 2k and 1k; where k is alphabet size). To improve performance for those cases you could allocating additional space at the end of suffix array.

|  file           |    size     |     libsais + O(n)  (ST)   |     libsais + O(1) (ST)    |speedup (ST)|    libsais + O(n)  (MT)    |     libsais + O(1) (MT)    |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            osdb |    10085684 |   0.201 sec (  50.18 MB/s) |   0.206 sec (  49.06 MB/s) |  **+2.29%**|   0.123 sec (  82.19 MB/s) |   0.133 sec (  75.84 MB/s) |  **+8.37%**|
|           x-ray |     8474240 |   0.175 sec (  48.36 MB/s) |   0.201 sec (  42.21 MB/s) | **+14.59%**|   0.103 sec (  82.37 MB/s) |   0.134 sec (  63.24 MB/s) | **+30.26%**|
|             sao |     7251944 |   0.160 sec (  45.40 MB/s) |   0.169 sec (  42.98 MB/s) |  **+5.63%**|   0.109 sec (  66.37 MB/s) |   0.125 sec (  58.17 MB/s) | **+14.09%**|
|         ooffice |     6152192 |   0.101 sec (  60.73 MB/s) |   0.106 sec (  58.27 MB/s) |  **+4.21%**|   0.069 sec (  88.67 MB/s) |   0.075 sec (  81.62 MB/s) |  **+8.63%**|
|           test3 |     2097088 |   0.031 sec (  67.91 MB/s) |   0.034 sec (  62.22 MB/s) |  **+9.16%**|   0.025 sec (  82.32 MB/s) |   0.030 sec (  70.56 MB/s) | **+16.66%**|

> * All other files from [Benchmarks](#benchmarks) above do not suffer from this fallbacks.

## libsais64 for inputs larger than 2GB

Starting from version 2.2.0 libsais64 could process inputs larger than 2GB.

The times below are the minimum of five runs measuring **multi-threaded (MT)** performance of suffix array construction on Azure DS14 v2 (Intel Xeon Platinum 8171M).

|  file           |    size     |    libsais64 2.2.0  (MT)   |   divsufsort64 2.0.2 (MT)  |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|
|         english |  2210395553 |  61.499 sec (  34.28 MB/s) | 435.199 sec (   4.84 MB/s) |**+607.65%**|
|   GRCh38.p13.fa |  3321586957 |  84.068 sec (  37.68 MB/s) | 782.938 sec (   4.05 MB/s) |**+831.32%**|
|         enwik10 | 10000000000 | 303.542 sec (  31.42 MB/s) |1927.351 sec (   4.95 MB/s) |**+534.95%**|
