/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <cstdint>

#include "PCCMath.h"

namespace pcc {

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformGbrToYCbCrBt709(T<Tv>& gbr)
{
  const Tv g = gbr[0];
  const Tv b = gbr[1];
  const Tv r = gbr[2];
  const double y =
    PCCClip(std::round(0.212600 * r + 0.715200 * g + 0.072200 * b), 0., 255.);
  const double u = PCCClip(
    std::round(-0.114572 * r - 0.385428 * g + 0.5 * b + 128.0), 0., 255.);
  const double v = PCCClip(
    std::round(0.5 * r - 0.454153 * g - 0.045847 * b + 128.0), 0., 255.);
  return {Tv(y), Tv(u), Tv(v)};
}

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformYCbCrBt709ToGbr(T<Tv>& ycbcr)
{
  const double y1 = ycbcr[0];
  const double u1 = ycbcr[1] - 128.0;
  const double v1 = ycbcr[2] - 128.0;
  const double r =
    PCCClip(round(y1 /*- 0.00000 * u1*/ + 1.57480 * v1), 0.0, 255.0);
  const double g =
    PCCClip(round(y1 - 0.18733 * u1 - 0.46813 * v1), 0.0, 255.0);
  const double b =
    PCCClip(round(y1 + 1.85563 * u1 /*+ 0.00000 * v1*/), 0.0, 255.0);
  return {Tv(g), Tv(b), Tv(r)};
}

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformGbrToYCgCoR(int bitDepth, T<Tv>& gbr)
{
  int g = gbr[0];
  int b = gbr[1];
  int r = gbr[2];

  int co = r - b;
  int t = b + (co >> 1);
  int cg = g - t;
  int y = t + (cg >> 1);

  int offset = 1 << bitDepth;

  // NB: YCgCoR needs extra 1-bit for chroma
  return {Tv(y), Tv(cg + offset), Tv(co + offset)};
}

//============================================================================

template<template<typename> class T, typename Tv>
T<Tv>
transformYCgCoRToGbr(int bitDepth, T<Tv>& ycgco)
{
  int offset = 1 << bitDepth;
  int y0 = ycgco[0];
  int cg = ycgco[1] - offset;
  int co = ycgco[2] - offset;

  int t = y0 - (cg >> 1);

  int g = cg + t;
  int b = t - (co >> 1);
  int r = co + b;

  int maxVal = (1 << bitDepth) - 1;
  g = PCCClip(g, 0, maxVal);
  b = PCCClip(b, 0, maxVal);
  r = PCCClip(r, 0, maxVal);

  return {Tv(g), Tv(b), Tv(r)};
}

//============================================================================

template<class T> 
inline static void convertRgbToYCbCrBt709(T (&src)[3]) { 
    T R = src[0]; T G = src[1]; T B = src[2]; 
    src[0] = static_cast<T>(+0.212600 * R + 0.715200 * G + 0.072200 * B); 
    src[1] = static_cast<T>(-0.114572 * R - 0.385428 * G + 0.500000 * B); 
    src[2] = static_cast<T>(+0.500000 * R - 0.454153 * G - 0.045847 * B); 
} 

template<class T> 
inline static void convertYCbCrToRgbBt709(T (&src)[3]) { 
    T Y = src[0]; T Cb = src[1]; T Cr = src[2]; 
    src[0] = static_cast<T>(Y + 1.574800 * Cr); 
    src[1] = static_cast<T>(Y - 0.187324 * Cb - 0.468124 * Cr); 
    src[2] = static_cast<T>(Y + 1.855600 * Cb); 
} 

template<class T> 
inline static void convertRgbToYCbCrBt601(T (&src)[3]) { 
    T R = src[0]; T G = src[1]; T B = src[2]; 
    src[0] = static_cast<T>(+0.299000 * R + 0.587000 * G + 0.114000 * B); 
    src[1] = static_cast<T>(-0.168736 * R - 0.331264 * G + 0.500000 * B); 
    src[2] = static_cast<T>(+0.500000 * R - 0.418688 * G - 0.081312 * B); 
} 

template<class T> 
inline static void convertYCbCrToRgbBt601(T (&src)[3]) { 
    T Y = src[0]; T Cb = src[1]; T Cr = src[2]; 
    src[0] = static_cast<T>(Y + 1.402000 * Cr); 
    src[1] = static_cast<T>(Y - 0.344136 * Cb - 0.714136 * Cr); 
    src[2] = static_cast<T>(Y + 1.772000 * Cb); 
} 

template<class T> 
inline void convertRgbToYcgco(T (&src)[3]) { 
    T R = src[0]; T G = src[1]; T B = src[2]; 
    src[0] = static_cast<T>( 0.25 * R + 0.5 * G + 0.25 * B); 
    src[1] = static_cast<T>(-0.25 * R + 0.5 * G - 0.25 * B); 
    src[2] = static_cast<T>( 0.50 * R - 0.5 * B); 
} 

template<class T> 
inline void convertYcgcoToRgb(T (&src)[3]) { 
    T Y = src[0]; T Cg = src[1]; T Co = src[2]; 
    src[0] = static_cast<T>(Y + Co - Cg); 
    src[1] = static_cast<T>(Y + Cg); 
    src[2] = static_cast<T>(Y - Co - Cg); 
} 

template<class T> 
inline void convertFromRgb(const ColourMatrix cicp_matrix_coefficients_idx, T (&src)[3]) { 
    switch (cicp_matrix_coefficients_idx) { 
        case ColourMatrix::kBt709: convertRgbToYCbCrBt709(src); break; 
        case ColourMatrix::kBt601: convertRgbToYCbCrBt601(src); break; 
        case ColourMatrix::kYCgCo: convertRgbToYcgco(src);      break; 
        default: break; 
    } 
} 

template<class T> 
inline void convertToRgb(const ColourMatrix cicp_matrix_coefficients_idx, T (&src)[3]) { 
    switch (cicp_matrix_coefficients_idx) { 
        case ColourMatrix::kBt709: convertYCbCrToRgbBt709(src); break; 
        case ColourMatrix::kBt601: convertYCbCrToRgbBt601(src); break; 
        case ColourMatrix::kYCgCo: convertYcgcoToRgb(src);      break; 
        default: break; 
    } 
}

//============================================================================

}  // namespace pcc
