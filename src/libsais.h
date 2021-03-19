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

#ifndef LIBSAIS_H
#define LIBSAIS_H 1

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif
