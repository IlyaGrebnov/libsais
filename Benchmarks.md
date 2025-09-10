# Specifications
  * OS: Windows 11 Pro (64-bit)
  * CPU: AMD Ryzen 9 9950X3D (16C / 32T, 128MB L3 cache, PBO +200Mhz, CO -10 all-core)
  * RAM: 2 x 48 GB DDR5-6000 (28-36-36-60)
  * Compiler: clang 19.1.5 '-O3 -march=znver5 -fopenmp -DLIBSAIS_OPENMP -DNDEBUG'

The times reflect the median of 15 runs for suffix array construction, reported in both **single-threaded (ST)** and **multi-threaded (MT)** modes. For optimal performance, libsais was limited to 12 threads, while divsufsort was limited to 16 threads.

### [Silesia Corpus](https://www.data-compression.info/Corpora/SilesiaCorpus/index.html) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|         dickens |    10192446 |   0.108 sec (  94.52 MB/s) |   0.297 sec (  34.33 MB/s) |**+175.35%**|   0.059 sec ( 173.57 MB/s) |   0.195 sec (  52.21 MB/s) |**+232.45%**|
|         mozilla |    51220480 |   0.542 sec (  94.55 MB/s) |   1.163 sec (  44.04 MB/s) |**+114.70%**|   0.318 sec ( 161.13 MB/s) |   0.725 sec (  70.69 MB/s) |**+127.94%**|
|              mr |     9970564 |   0.095 sec ( 104.55 MB/s) |   0.239 sec (  41.66 MB/s) |**+150.96%**|   0.069 sec ( 145.10 MB/s) |   0.163 sec (  61.28 MB/s) |**+136.81%**|
|             nci |    33553445 |   0.274 sec ( 122.48 MB/s) |   0.682 sec (  49.20 MB/s) |**+148.97%**|   0.142 sec ( 235.87 MB/s) |   0.605 sec (  55.43 MB/s) |**+325.56%**|
|         ooffice |     6152192 |   0.072 sec (  85.36 MB/s) |   0.130 sec (  47.18 MB/s) | **+80.94%**|   0.052 sec ( 119.34 MB/s) |   0.081 sec (  76.34 MB/s) | **+56.32%**|
|            osdb |    10085684 |   0.110 sec (  91.44 MB/s) |   0.214 sec (  47.04 MB/s) | **+94.40%**|   0.073 sec ( 138.57 MB/s) |   0.168 sec (  60.01 MB/s) |**+130.91%**|
|         reymont |     6627202 |   0.069 sec (  95.98 MB/s) |   0.173 sec (  38.36 MB/s) |**+150.23%**|   0.042 sec ( 157.05 MB/s) |   0.129 sec (  51.38 MB/s) |**+205.66%**|
|           samba |    21606400 |   0.208 sec ( 103.82 MB/s) |   0.440 sec (  49.15 MB/s) |**+111.24%**|   0.133 sec ( 163.02 MB/s) |   0.304 sec (  71.01 MB/s) |**+129.56%**|
|             sao |     7251944 |   0.107 sec (  67.75 MB/s) |   0.154 sec (  47.08 MB/s) | **+43.89%**|   0.079 sec (  91.84 MB/s) |   0.095 sec (  75.99 MB/s) | **+20.85%**|
|         webster |    41458703 |   0.445 sec (  93.24 MB/s) |   1.335 sec (  31.04 MB/s) |**+200.33%**|   0.243 sec ( 170.50 MB/s) |   0.944 sec (  43.93 MB/s) |**+288.10%**|
|           x-ray |     8474240 |   0.139 sec (  60.85 MB/s) |   0.219 sec (  38.61 MB/s) | **+57.62%**|   0.093 sec (  91.56 MB/s) |   0.106 sec (  80.02 MB/s) | **+14.43%**|
|             xml |     5345280 |   0.044 sec ( 122.71 MB/s) |   0.097 sec (  55.15 MB/s) |**+122.52%**|   0.029 sec ( 182.47 MB/s) |   0.077 sec (  69.57 MB/s) |**+162.29%**|



### [Large Canterbury Corpus](https://www.data-compression.info/Corpora/CanterburyCorpus/) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       bible.txt |     4047392 |   0.040 sec ( 100.99 MB/s) |   0.106 sec (  38.15 MB/s) |**+164.69%**|   0.026 sec ( 156.84 MB/s) |   0.075 sec (  53.85 MB/s) |**+191.28%**|
|          E.coli |     4638690 |   0.052 sec (  89.05 MB/s) |   0.156 sec (  29.70 MB/s) |**+199.79%**|   0.024 sec ( 192.42 MB/s) |   0.123 sec (  37.71 MB/s) |**+410.27%**|
|    world192.txt |     2473400 |   0.023 sec ( 108.73 MB/s) |   0.056 sec (  43.92 MB/s) |**+147.54%**|   0.017 sec ( 149.17 MB/s) |   0.037 sec (  66.84 MB/s) |**+123.17%**|



### [Manzini Corpus](https://people.unipmn.it/manzini/lightweight/corpus/) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|       chr22.dna |    34553758 |   0.400 sec (  86.39 MB/s) |   1.252 sec (  27.61 MB/s) |**+212.96%**|   0.177 sec ( 195.37 MB/s) |   1.024 sec (  33.73 MB/s) |**+479.24%**|
|         etext99 |   105277340 |   1.341 sec (  78.51 MB/s) |   4.171 sec (  25.24 MB/s) |**+211.02%**|   0.776 sec ( 135.69 MB/s) |   2.855 sec (  36.87 MB/s) |**+268.00%**|
|     gcc-3.0.tar |    86630400 |   0.909 sec (  95.30 MB/s) |   2.332 sec (  37.15 MB/s) |**+156.53%**|   0.529 sec ( 163.65 MB/s) |   1.614 sec (  53.67 MB/s) |**+204.93%**|
|           howto |    39422105 |   0.435 sec (  90.63 MB/s) |   1.235 sec (  31.91 MB/s) |**+183.99%**|   0.237 sec ( 166.57 MB/s) |   0.794 sec (  49.66 MB/s) |**+235.45%**|
|          jdk13c |    69728899 |   0.692 sec ( 100.74 MB/s) |   1.864 sec (  37.40 MB/s) |**+169.37%**|   0.412 sec ( 169.13 MB/s) |   1.495 sec (  46.63 MB/s) |**+262.72%**|
| linux-2.4.5.tar |   116254720 |   1.281 sec (  90.74 MB/s) |   3.380 sec (  34.39 MB/s) |**+163.81%**|   0.898 sec ( 129.48 MB/s) |   2.409 sec (  48.26 MB/s) |**+168.31%**|
|        rctail96 |   114711151 |   1.314 sec (  87.30 MB/s) |   3.790 sec (  30.27 MB/s) |**+188.41%**|   0.935 sec ( 122.65 MB/s) |   2.945 sec (  38.95 MB/s) |**+214.89%**|
|             rfc |   116421901 |   1.269 sec (  91.75 MB/s) |   3.601 sec (  32.33 MB/s) |**+183.84%**|   0.887 sec ( 131.32 MB/s) |   2.626 sec (  44.33 MB/s) |**+196.21%**|
|     sprot34.dat |   109617186 |   1.242 sec (  88.26 MB/s) |   3.739 sec (  29.32 MB/s) |**+201.01%**|   0.848 sec ( 129.27 MB/s) |   2.532 sec (  43.30 MB/s) |**+198.56%**|
|            w3c2 |   104201579 |   1.083 sec (  96.21 MB/s) |   2.848 sec (  36.59 MB/s) |**+162.95%**|   0.809 sec ( 128.80 MB/s) |   2.326 sec (  44.80 MB/s) |**+187.50%**|



### [Large Text Compression Benchmark Corpus](https://www.mattmahoney.net/dc/textdata.html) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|          enwik8 |   100000000 |   1.224 sec (  81.69 MB/s) |   3.785 sec (  26.42 MB/s) |**+209.19%**|   0.827 sec ( 120.98 MB/s) |   2.568 sec (  38.94 MB/s) |**+210.67%**|
|          enwik9 |  1000000000 |  13.633 sec (  73.35 MB/s) |  42.792 sec (  23.37 MB/s) |**+213.88%**|   9.410 sec ( 106.27 MB/s) |  29.223 sec (  34.22 MB/s) |**+210.56%**|



### [The Gauntlet Corpus](https://github.com/michaelmaniscalco/gauntlet_corpus) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            abac |      200000 |   0.002 sec ( 120.88 MB/s) |   0.001 sec ( 156.91 MB/s) |   -22.97%  |   0.001 sec ( 193.69 MB/s) |   0.001 sec ( 151.70 MB/s) | **+27.68%**|
|            abba |    10500596 |   0.070 sec ( 150.70 MB/s) |   0.267 sec (  39.27 MB/s) |**+283.78%**|   0.046 sec ( 228.50 MB/s) |   0.297 sec (  35.37 MB/s) |**+546.11%**|
|        book1x20 |    15375420 |   0.138 sec ( 111.81 MB/s) |   0.477 sec (  32.21 MB/s) |**+247.16%**|   0.086 sec ( 178.20 MB/s) |   0.330 sec (  46.62 MB/s) |**+282.23%**|
|   fib_s14930352 |    14930352 |   0.122 sec ( 122.45 MB/s) |   0.537 sec (  27.79 MB/s) |**+340.60%**|   0.082 sec ( 183.15 MB/s) |   0.594 sec (  25.14 MB/s) |**+628.64%**|
|           fss10 |    12078908 |   0.091 sec ( 133.24 MB/s) |   0.406 sec (  29.73 MB/s) |**+348.11%**|   0.063 sec ( 192.55 MB/s) |   0.463 sec (  26.10 MB/s) |**+637.75%**|
|            fss9 |     2851443 |   0.026 sec ( 109.23 MB/s) |   0.098 sec (  29.10 MB/s) |**+275.34%**|   0.015 sec ( 192.40 MB/s) |   0.104 sec (  27.47 MB/s) |**+600.49%**|
|         houston |     3839141 |   0.023 sec ( 165.51 MB/s) |   0.018 sec ( 212.46 MB/s) |   -22.10%  |   0.012 sec ( 333.62 MB/s) |   0.018 sec ( 219.29 MB/s) | **+52.14%**|
|       paper5x80 |      956322 |   0.009 sec ( 105.78 MB/s) |   0.018 sec (  53.26 MB/s) | **+98.60%**|   0.006 sec ( 173.22 MB/s) |   0.014 sec (  67.59 MB/s) |**+156.27%**|
|           test1 |     2097152 |   0.026 sec (  81.67 MB/s) |   0.042 sec (  50.29 MB/s) | **+62.39%**|   0.014 sec ( 147.02 MB/s) |   0.041 sec (  51.48 MB/s) |**+185.60%**|
|           test2 |     2097152 |   0.026 sec (  81.80 MB/s) |   0.032 sec (  66.13 MB/s) | **+23.70%**|   0.014 sec ( 149.18 MB/s) |   0.032 sec (  66.40 MB/s) |**+124.69%**|
|           test3 |     2097088 |   0.023 sec (  91.81 MB/s) |   0.034 sec (  62.20 MB/s) | **+47.61%**|   0.024 sec (  87.59 MB/s) |   0.038 sec (  55.21 MB/s) | **+58.65%**|



### [Pizza & Chilli Corpus](https://pizzachili.dcc.uchile.cl/texts.html) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|        dblp.xml |   296135874 |   3.422 sec (  86.54 MB/s) |   9.243 sec (  32.04 MB/s) |**+170.10%**|   2.568 sec ( 115.33 MB/s) |   7.178 sec (  41.26 MB/s) |**+179.53%**|
|             dna |   403927746 |   5.691 sec (  70.98 MB/s) |  19.028 sec (  21.23 MB/s) |**+234.35%**|   3.476 sec ( 116.21 MB/s) |  16.365 sec (  24.68 MB/s) |**+370.82%**|
|  english.1024MB |  1073741824 |  16.110 sec (  66.65 MB/s) |  52.004 sec (  20.65 MB/s) |**+222.81%**|  11.219 sec (  95.71 MB/s) |  37.177 sec (  28.88 MB/s) |**+231.37%**|
|         pitches |    55832855 |   0.690 sec (  80.96 MB/s) |   1.656 sec (  33.71 MB/s) |**+140.12%**|   0.462 sec ( 120.98 MB/s) |   1.024 sec (  54.51 MB/s) |**+121.96%**|
|        proteins |  1184051855 |  18.368 sec (  64.46 MB/s) |  66.699 sec (  17.75 MB/s) |**+263.13%**|  12.717 sec (  93.10 MB/s) |  35.746 sec (  33.12 MB/s) |**+181.08%**|
|         sources |   210866607 |   2.469 sec (  85.40 MB/s) |   6.702 sec (  31.46 MB/s) |**+171.45%**|   1.780 sec ( 118.43 MB/s) |   4.741 sec (  44.48 MB/s) |**+166.27%**|



### [Pizza & Chilli Repetitive Corpus](https://pizzachili.dcc.uchile.cl/repcorpus.html) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            cere |   461286644 |   5.356 sec (  86.12 MB/s) |  19.154 sec (  24.08 MB/s) |**+257.61%**|   3.597 sec ( 128.23 MB/s) |  16.417 sec (  28.10 MB/s) |**+356.39%**|
|       coreutils |   205281778 |   2.181 sec (  94.11 MB/s) |   6.814 sec (  30.13 MB/s) |**+212.39%**|   1.626 sec ( 126.26 MB/s) |   5.216 sec (  39.36 MB/s) |**+220.81%**|
| einstein.de.txt |    92758441 |   0.926 sec ( 100.22 MB/s) |   2.676 sec (  34.66 MB/s) |**+189.17%**|   0.635 sec ( 146.08 MB/s) |   2.167 sec (  42.81 MB/s) |**+241.27%**|
| einstein.en.txt |   467626544 |   5.480 sec (  85.33 MB/s) |  15.822 sec (  29.56 MB/s) |**+188.72%**|   3.885 sec ( 120.36 MB/s) |  12.593 sec (  37.13 MB/s) |**+224.13%**|
|Escherichia_Coli |   112689515 |   1.272 sec (  88.61 MB/s) |   4.562 sec (  24.70 MB/s) |**+258.72%**|   0.864 sec ( 130.37 MB/s) |   3.827 sec (  29.44 MB/s) |**+342.80%**|
|       influenza |   154808555 |   1.588 sec (  97.51 MB/s) |   5.325 sec (  29.07 MB/s) |**+235.42%**|   1.159 sec ( 133.60 MB/s) |   4.619 sec (  33.51 MB/s) |**+298.67%**|
|          kernel |   257961616 |   2.747 sec (  93.90 MB/s) |   8.993 sec (  28.69 MB/s) |**+227.33%**|   2.030 sec ( 127.09 MB/s) |   6.614 sec (  39.00 MB/s) |**+225.88%**|
|            para |   429265758 |   5.129 sec (  83.69 MB/s) |  18.302 sec (  23.45 MB/s) |**+256.81%**|   3.444 sec ( 124.63 MB/s) |  15.618 sec (  27.49 MB/s) |**+353.45%**|
|   world_leaders |    46968181 |   0.338 sec ( 138.85 MB/s) |   0.758 sec (  62.00 MB/s) |**+123.96%**|   0.202 sec ( 232.62 MB/s) |   0.637 sec (  73.76 MB/s) |**+215.37%**|
|dblp.xml.00001.1 |   104857600 |   1.944 sec (  53.94 MB/s) |   3.809 sec (  27.53 MB/s) | **+95.95%**|   0.848 sec ( 123.66 MB/s) |   3.339 sec (  31.40 MB/s) |**+293.78%**|
|dblp.xml.00001.2 |   104857600 |   1.922 sec (  54.57 MB/s) |   3.825 sec (  27.42 MB/s) | **+99.02%**|   0.845 sec ( 124.11 MB/s) |   3.358 sec (  31.22 MB/s) |**+297.49%**|
| dblp.xml.0001.1 |   104857600 |   1.876 sec (  55.89 MB/s) |   3.753 sec (  27.94 MB/s) |**+100.03%**|   0.845 sec ( 124.08 MB/s) |   3.272 sec (  32.05 MB/s) |**+287.18%**|
| dblp.xml.0001.2 |   104857600 |   1.840 sec (  56.98 MB/s) |   3.743 sec (  28.01 MB/s) |**+103.42%**|   0.844 sec ( 124.26 MB/s) |   3.242 sec (  32.35 MB/s) |**+284.16%**|
|       dna.001.1 |   104857600 |   1.957 sec (  53.59 MB/s) |   4.475 sec (  23.43 MB/s) |**+128.73%**|   0.785 sec ( 133.57 MB/s) |   3.835 sec (  27.34 MB/s) |**+388.49%**|
|   english.001.2 |   104857600 |   2.035 sec (  51.53 MB/s) |   4.360 sec (  24.05 MB/s) |**+114.26%**|   0.914 sec ( 114.76 MB/s) |   3.242 sec (  32.35 MB/s) |**+254.80%**|
|  proteins.001.1 |   104857600 |   2.001 sec (  52.41 MB/s) |   4.486 sec (  23.37 MB/s) |**+124.24%**|   0.871 sec ( 120.43 MB/s) |   2.880 sec (  36.40 MB/s) |**+230.82%**|
|   sources.001.2 |   104857600 |   1.777 sec (  59.02 MB/s) |   3.649 sec (  28.74 MB/s) |**+105.40%**|   0.847 sec ( 123.79 MB/s) |   2.940 sec (  35.67 MB/s) |**+247.04%**|
|           fib41 |   267914296 |   3.362 sec (  79.68 MB/s) |  16.837 sec (  15.91 MB/s) |**+400.76%**|   2.082 sec ( 128.66 MB/s) |  17.004 sec (  15.76 MB/s) |**+716.58%**|
|           rs.13 |   216747218 |   2.292 sec (  94.57 MB/s) |  13.105 sec (  16.54 MB/s) |**+471.82%**|   1.671 sec ( 129.71 MB/s) |  13.311 sec (  16.28 MB/s) |**+696.59%**|
|            tm29 |   268435456 |   8.414 sec (  31.90 MB/s) |  19.460 sec (  13.79 MB/s) |**+131.29%**|   2.407 sec ( 111.51 MB/s) |  20.004 sec (  13.42 MB/s) |**+730.96%**|



### [Skyline Corpus](http://panthema.net/2012/1119-eSAIS-Inducing-Suffix-and-LCP-Arrays-in-External-Memory/eSAIS-DC3-LCP-0.5.0/src/input/skyline.h.html) ###

|  file           |    size     |     libsais 2.10.4 (ST)    |    divsufsort 2.0.2 (ST)   |speedup (ST)|     libsais 2.10.4 (MT)    |    divsufsort 2.0.2 (MT)   |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|   skyline20.txt |     1048575 |   0.016 sec (  64.17 MB/s) |   0.055 sec (  18.92 MB/s) |**+239.18%**|   0.011 sec (  97.18 MB/s) |   0.055 sec (  19.15 MB/s) |**+407.38%**|
|   skyline22.txt |     4194303 |   0.068 sec (  61.81 MB/s) |   0.260 sec (  16.16 MB/s) |**+282.46%**|   0.033 sec ( 127.33 MB/s) |   0.276 sec (  15.20 MB/s) |**+737.65%**|
|   skyline24.txt |    16777215 |   0.375 sec (  44.76 MB/s) |   1.355 sec (  12.38 MB/s) |**+261.57%**|   0.140 sec ( 119.44 MB/s) |   1.376 sec (  12.19 MB/s) |**+879.45%**|
|   skyline26.txt |    67108863 |   2.862 sec (  23.45 MB/s) |   7.519 sec (   8.93 MB/s) |**+162.71%**|   0.772 sec (  86.97 MB/s) |   7.566 sec (   8.87 MB/s) |**+880.50%**|
|   skyline28.txt |   268435455 |  12.159 sec (  22.08 MB/s) |  36.970 sec (   7.26 MB/s) |**+204.06%**|   3.827 sec (  70.15 MB/s) |  36.876 sec (   7.28 MB/s) |**+863.64%**|

## Additional memory

The libsais reuses the space allocated for the suffix array during construction. In rare cases, this space is not large enough for the fastest algorithm and libsais will need to fallback to less efficient one (libsais has 4 algorithms at different break-points point: 6k, 4k, 2k and 1k; where k is alphabet size). To improve performance for those cases you could allocating additional space at the end of suffix array.

|  file           |    size     |     libsais + O(n)  (ST)   |     libsais + O(1) (ST)    |speedup (ST)|    libsais + O(n)  (MT)    |     libsais + O(1) (MT)    |speedup (MT)|
|:---------------:|:-----------:|:--------------------------:|:--------------------------:|:----------:|:--------------------------:|:--------------------------:|:----------:|
|            osdb |    10085684 |   0.100 sec ( 100.68 MB/s) |   0.107 sec (  94.43 MB/s) |  **+6.61%**|   0.071 sec ( 142.67 MB/s) |   0.078 sec ( 128.57 MB/s) | **+10.96%**|
|           x-ray |     8474240 |   0.112 sec (  75.75 MB/s) |   0.138 sec (  61.20 MB/s) | **+23.77%**|   0.070 sec ( 121.65 MB/s) |   0.093 sec (  91.46 MB/s) | **+33.01%**|
|             sao |     7251944 |   0.091 sec (  80.02 MB/s) |   0.106 sec (  68.19 MB/s) | **+17.35%**|   0.068 sec ( 107.05 MB/s) |   0.080 sec (  90.48 MB/s) | **+18.32%**|
|         ooffice |     6152192 |   0.065 sec (  94.39 MB/s) |   0.071 sec (  86.20 MB/s) |  **+9.49%**|   0.048 sec ( 127.93 MB/s) |   0.055 sec ( 111.39 MB/s) | **+14.85%**|
|           test3 |     2097088 |   0.017 sec ( 125.08 MB/s) |   0.019 sec ( 109.88 MB/s) | **+13.83%**|   0.015 sec ( 135.61 MB/s) |   0.018 sec ( 114.18 MB/s) | **+18.76%**|

> * All other files from [Benchmarks](#benchmarks) above do not suffer from this fallbacks.
