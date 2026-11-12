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

#include "ply.h"

#include "PCCMisc.h"
#include "PCCPointSet.h"
#include <string>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>

namespace pcc {

const AttributeDescription& getAttributeDescription(
  const std::vector<AttributeDescription>& attrDescs,
  KnownAttributeLabel attributeLabel, int index);

AttributeDescription& getAttributeDescription(
  std::vector<AttributeDescription>& attrDescs,
  KnownAttributeLabel attributeLabel, int index);

//============================================================================

int64_t adjust_frac_bits(int64_t val, int curr_frac_bits, int final_frac_bits)
{
    if (curr_frac_bits < final_frac_bits) {
        int shift = final_frac_bits - curr_frac_bits;
        val = val << shift;
    }
    else {
        int shift = curr_frac_bits - final_frac_bits;
        val = val >> shift;
    }
    return val;
}

//============================================================================
// Clamp used on either side of the opacity logit<->linear transform to avoid
// log(0)/log(inf) at the domain boundaries.
static constexpr double kOpacityLogitEpsilon = 1e-6;

static double opacityLogitToLinear(double logit)
{
  return 1.0 / (1.0 + std::exp(-logit));
}

static double opacityLinearToLogit(double linear)
{
  double clipped = PCCClip(linear, kOpacityLogitEpsilon, 1.0 - kOpacityLogitEpsilon);
  return std::log(clipped / (1.0 - clipped));
}

std::vector<int> getAttrIndices(
  const AttributeDescription& desc,
  const std::vector<ply::PropertyInfo>& properties);
static double computeShAcTop10Correlation(SHBufferView& shView, size_t pointCount) {
  int numSHCoeffs = shView.numBands();
  if (pointCount == 0 || numSHCoeffs <= 1) return 0.0;

  const int acCount = numSHCoeffs - 1;  // bands 1..15
  std::vector<double> acMean(acCount, 0.0);
  for (size_t p = 0; p < pointCount; ++p) {
    for (int l = 0; l < acCount; ++l) {
      acMean[l] += (shView.sh(p, l + 1, 0) + shView.sh(p, l + 1, 1) + shView.sh(p, l + 1, 2)) / 3.0;
    }
  }
  for (int l = 0; l < acCount; ++l) acMean[l] /= (double)pointCount;

  std::vector<double> var(acCount, 0.0);
  std::vector<double> cov(acCount * acCount, 0.0);
  for (size_t p = 0; p < pointCount; ++p) {
    std::vector<double> centered(acCount);
    for (int l = 0; l < acCount; ++l) {
      centered[l] = (shView.sh(p, l + 1, 0) + shView.sh(p, l + 1, 1) + shView.sh(p, l + 1, 2)) / 3.0 - acMean[l];
    }
    for (int r = 0; r < acCount; ++r) {
      var[r] += centered[r] * centered[r];
      for (int c = r + 1; c < acCount; ++c) cov[r * acCount + c] += centered[r] * centered[c];
    }
  }
  for (int r = 0; r < acCount; ++r) {
    var[r] /= (double)(pointCount > 1 ? pointCount - 1 : 1);
    for (int c = r + 1; c < acCount; ++c)
      cov[r * acCount + c] /= (double)(pointCount > 1 ? pointCount - 1 : 1);
  }

  std::vector<double> absCorrs;
  for (int r = 0; r < acCount; ++r) {
    for (int c = r + 1; c < acCount; ++c) {
      double denom = std::sqrt(std::max(1e-20, var[r] * var[c]));
      absCorrs.push_back(std::abs(cov[r * acCount + c] / denom));
    }
  }
  std::sort(absCorrs.rbegin(), absCorrs.rend());
  int topN = std::min(10, (int)absCorrs.size());
  double sum = 0.0;
  for (int i = 0; i < topN; ++i) sum += absCorrs[i];
  return topN > 0 ? sum / topN : 0.0;
}

//============================================================================
bool extractSHQuantRegions(bool withSHCoefficients, SHBufferView& shView,
    size_t pointCount, std::vector<AttributeDescription>& attrDescs) {  if (pointCount == 0 || !withSHCoefficients) return false;

  int numSHCoeffs = shView.numBands();
  auto& firstShParams =
    getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, 0).params;
  const bool shAcScalingEnabled = firstShParams.attr_ac_scaling_enabled_flag;

  const bool useAdaptiveThresholds = firstShParams.shAdaptiveTpOptEnabled_
                                   && firstShParams.dataAdaptiveShTransformEnabled_;

  int thresh1 = 15, thresh2 = 30, thresh3 = 40, thresh4 = 5, thresh5 = 10;
  const int idx3_arr[3] = { 11, 13, 12 };
  const int idx2_arr[3] = { 10, 14, 15 };

  if (useAdaptiveThresholds) {
    double sum1 = 0.0, sumSq1 = 0.0;
    double sum3 = 0.0, sumSq3 = 0.0;
    const size_t count1 = pointCount * 9 * 3;
    const size_t count3 = pointCount * 3 * 3;

    for (size_t p = 0; p < pointCount; ++p) {
      for (size_t l = 1; l <= 9; ++l) {
        for (int c = 0; c < 3; ++c) {
          double val = shView.sh(p, l, c);
          sum1 += val;
          sumSq1 += val * val;
        }
      }
      for (int i = 0; i < 3; ++i) {
        size_t l = idx3_arr[i];
        for (int c = 0; c < 3; ++c) {
          double val = shView.sh(p, l, c);
          sum3 += val;
          sumSq3 += val * val;
        }
      }
    }

    double mean1 = sum1 / count1;
    double mean3 = sum3 / count3;
    double var1 = (sumSq1 / count1) - (mean1 * mean1);
    double var3 = (sumSq3 / count3) - (mean3 * mean3);
    double std1 = std::sqrt(std::max(0.0, var1));
    double std3 = std::sqrt(std::max(0.0, var3));

    double stdRatio = std3 > 1e-10 ? std1 / std3 : 0.0;
    if (stdRatio > 2.0 && stdRatio < 10.0) {
      thresh1 = 10; thresh2 = 15; thresh3 = 20; thresh4 = 25; thresh5 = 30;
    }
  }

  // NQ skip: signal non-retained channels with M=0.
  std::vector<bool> shDrRetainedMask;
  const bool doNqSkip = attrDescs[0].params.enableShDimensionReduction_
    && attrDescs[0].params.enableShDimensionReductionNqSkipNonretained_
    && !attrDescs[0].params.shNumRetainComponents_.empty();
  if (doNqSkip)
    buildShDimRedRetainedMask(attrDescs[0].params, shDrRetainedMask);

  for (size_t l = 0; l < std::min(numSHCoeffs, 16); ++l) {
    auto& params =
      getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, l).params;
    const uint32_t bitDepth =
      (shAcScalingEnabled && l > 0) ? params.attr_frac_bits : firstShParams.attr_frac_bits;
    const uint32_t max_quant_val = (1u << bitDepth) - 1;

    for (int c = 0 ; c < 3 ; c++)
    {
      params.shNumRegionsMinus1_[c] = 2;
      params.shNumBitsForTransitionPoint_[c] = bitDepth;
      params.shTransitionPoints_[c].clear();
      params.shStepSizes_[c].clear();
      params.shStepSizesFxp_[c].clear();
      
      if (doNqSkip && !shDrRetainedMask[l * 3 + c]) {
        params.shNumRegionsMinus1_[c] = 0;
        params.min_offset_fxp[c] = 0; 
        double minLim = params.attr_offset;
        double maxLim = params.attr_offset + params.attr_scale_minus1 + 1;
        double d_code = std::round((maxLim - minLim) * (1 << bitDepth) /
                                   (double)(1ull << params.attr_frac_bits));
        double s_uniform = d_code / (double)(max_quant_val + 1);
        params.shStepSizes_[c].push_back(s_uniform);
        int64_t s_uniformFxp = (int64_t)round(s_uniform * pow(2.0, SS_SCALE));
        if (s_uniformFxp == 0) s_uniformFxp = 1;
        params.shStepSizesFxp_[c].push_back(s_uniformFxp);
        continue; 
      }

      std::vector<double> values;

      values.reserve(pointCount);
      for (size_t p = 0; p < pointCount; ++p) {
        values.push_back(shView.sh(p, l, c));
      }

      double tp_low, tp_high;

      if (useAdaptiveThresholds) {
        double perc_low, perc_high;
        if ( l == 0 ) {
          if ( c == 0 ) {
            perc_low  = thresh4;
            perc_high = 100 - thresh4;
          } else {
            perc_low  = thresh5;
            perc_high = 100 - thresh5;
          }
        } else if ( l >= 1 && l <= 9 ) {
          perc_low  = thresh1;
          perc_high = 100 - thresh1;
        } else if ( l == static_cast<size_t>( idx2_arr[0] ) || l == static_cast<size_t>( idx2_arr[1] ) ||
                    l == static_cast<size_t>( idx2_arr[2] ) ) {
          perc_low  = thresh2;
          perc_high = 100 - thresh2;
        } else if ( l == static_cast<size_t>( idx3_arr[0] ) || l == static_cast<size_t>( idx3_arr[1] ) ||
                    l == static_cast<size_t>( idx3_arr[2] ) ) {
          perc_low  = thresh3;
          perc_high = 100 - thresh3;
        } else {
          perc_low  = 25.0;
          perc_high = 75.0;
        }

        size_t lowIdx  = static_cast<size_t>( std::round( perc_low / 100.0 * pointCount ) );
        size_t highIdx = static_cast<size_t>( std::round( perc_high / 100.0 * pointCount ) );
        lowIdx  = std::min( lowIdx, pointCount > 0 ? pointCount - 1 : 0 );
        highIdx = std::min( highIdx, pointCount > 0 ? pointCount - 1 : 0 );

        std::nth_element( values.begin(), values.begin() + lowIdx, values.end() );
        tp_low = values[lowIdx];

        std::nth_element( values.begin(), values.begin() + highIdx, values.end() );
        tp_high = values[highIdx];
      } else {
        double ratio = 0.25;
        size_t first_transition_index = size_t(pointCount * ratio);
        size_t second_transition_index = size_t (pointCount * (1 - ratio));
        std::nth_element( values.begin(), values.begin() + first_transition_index, values.end() );
        tp_low = values[first_transition_index];

        std::nth_element( values.begin(), values.begin() + second_transition_index, values.end() );
        tp_high = values[second_transition_index];
      }

      double minLim = params.attr_offset;
      double maxLim = params.attr_offset + params.attr_scale_minus1 + 1;
      double Range = maxLim - minLim;
      double d_code = std::round(Range * (1 << bitDepth) / (double)(1ull << params.attr_frac_bits));

      double s_fine = d_code / (double)(max_quant_val+1);
      double s_coarse = d_code / (double)((max_quant_val + 1) / 2);

      int tp_low_quant  = static_cast<int>( std::round( ( tp_low - minLim ) / d_code * max_quant_val ) );
      int tp_high_quant = static_cast<int>( std::round( ( tp_high - minLim ) / d_code * max_quant_val ) );

      params.shTransitionPoints_[c].push_back( std::max( 0, std::min( static_cast<int>( max_quant_val - 1 ), tp_low_quant ) ) );
      params.shTransitionPoints_[c].push_back( std::max( tp_low_quant + 1, std::min( static_cast<int>( max_quant_val ), tp_high_quant ) ) );

      params.shStepSizes_[c].push_back(s_coarse);
      params.shStepSizes_[c].push_back(s_fine);
      params.shStepSizes_[c].push_back(s_coarse);
      const double ss_scale = pow(2.0, SS_SCALE);
      int64_t s_coarseFxp = (int64_t)round(s_coarse * ss_scale);
      int64_t s_fineFxp = (int64_t)round(s_fine * ss_scale);
      if (s_fineFxp == 0) s_fineFxp = 1;
      if (s_coarseFxp == 0) s_coarseFxp = 1;
      params.shStepSizesFxp_[c].push_back(s_coarseFxp);
      params.shStepSizesFxp_[c].push_back(s_fineFxp);
      params.shStepSizesFxp_[c].push_back(s_coarseFxp);
    }
  }

  return true;
}

static void
permut_02(
  const double quat_in[4],
  const float scale_in[3],
  double quat_out[4],
  float scale_out[3])
{
  scale_out[0] = scale_in[2];
  scale_out[1] = scale_in[1];
  scale_out[2] = scale_in[0];
  quat_out[0] = quat_in[0] - quat_in[2];
  quat_out[1] = quat_in[1] - quat_in[3];
  quat_out[2] = quat_in[2] + quat_in[0];
  quat_out[3] = quat_in[3] + quat_in[1];
}

static void
permut_12(
  const double quat_in[4],
  const float scale_in[3],
  double quat_out[4],
  float scale_out[3])
{
  scale_out[0] = scale_in[0];
  scale_out[1] = scale_in[2];
  scale_out[2] = scale_in[1];
  quat_out[0] = quat_in[0] - quat_in[1];
  quat_out[1] = quat_in[1] + quat_in[0];
  quat_out[2] = quat_in[2] + quat_in[3];
  quat_out[3] = quat_in[3] - quat_in[2];
}

static void
canonicalize(double quat[4], float scale[3])
{
  // Apply the 6 scales permutation, and search for quaternion with maximum
  // absolute value component
  double quats[6][4];
  float scales[6][3];

  std::copy_n(quat, 4, quats[0]);
  std::copy_n(scale, 3, scales[0]);
  permut_12(quats[0], scales[0], quats[1], scales[1]);
  permut_02(quats[1], scales[1], quats[2], scales[2]);
  permut_12(quats[2], scales[2], quats[3], scales[3]);
  permut_02(quats[3], scales[3], quats[4], scales[4]);
  permut_12(quats[4], scales[4], quats[5], scales[5]);

  double maxval = -1;
  int imax = 0;
  constexpr double mul[6] = {   // These values is to take into account that when
    1.00, 0.7071067811865476,   // the functions permut_02 or permut_12 are used,
    0.50, 0.3535533905932738,   // the resulting quaternion will have a norm equal
    0.25, 0.1767766952966369 }; // to the input quaternion norm times sqrt(2)
  for (int i = 0; i < 6; ++i) {
    double max = mul[i] * std::max(
      std::max(std::abs(quats[i][0]), std::abs(quats[i][1])),
      std::max(std::abs(quats[i][2]), std::abs(quats[i][3])));
    if (maxval < max) {
      maxval = max;
      imax = i;
    }
  }
  std::copy_n(scales[imax], 3, scale);
  if (maxval < 1e-6) {
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    return;
  }
  
  // Put the maximum absolute value component of quaternion as the real part
  double w = std::abs(quats[imax][0]);
  double x = std::abs(quats[imax][1]);
  double y = std::abs(quats[imax][2]);
  double z = std::abs(quats[imax][3]);
  if (w >= x && w >= y && w >= z) {
    // copy
    quat[0] = +quats[imax][0];
    quat[1] = +quats[imax][1];
    quat[2] = +quats[imax][2];
    quat[3] = +quats[imax][3];
  } else if (x >= y && x >= z) {
    // rotate -pi radians about x
    quat[0] = +quats[imax][1];
    quat[1] = -quats[imax][0];
    quat[2] = -quats[imax][3];
    quat[3] = +quats[imax][2];
  } else if (y >= z) {
    // rotate -pi radians about y
    quat[0] = +quats[imax][2];
    quat[1] = +quats[imax][3];
    quat[2] = -quats[imax][0];
    quat[3] = -quats[imax][1];
  } else {
    // rotate -pi radians about z
    quat[0] = +quats[imax][3];
    quat[1] = -quats[imax][2];
    quat[2] = +quats[imax][1];
    quat[3] = -quats[imax][0];
  }
}

static void normalizeQuat(double quat[4], bool isRealRotCoded, bool isEulerRot)
{
  // 1. Handle special rotation encoding case first
  if (!isEulerRot && isRealRotCoded) {
    // If it's a real rotation coded format, use quat[0] directly as norm_sq
    double norm_sq = quat[0];
    
    if (norm_sq > 0.0) {
      if (quat[0] < 0.0 || std::abs(norm_sq - 1.0) > 1e-12) {
        double inv_norm = 1.0 / std::sqrt(norm_sq);
        if (quat[0] < 0.0) {
          inv_norm = -inv_norm;
        }
        // Mathematically verified: setting quat[0] to 1.0 restores proper scaling
        quat[0] = 1.0;
        quat[1] *= inv_norm;
        quat[2] *= inv_norm;
        quat[3] *= inv_norm;
      }
    }
    return; // Early exit to skip standard normalization
  }

  // 2. Standard Quaternion Normalization
  double norm_sq = quat[0] * quat[0] 
                 + quat[1] * quat[1] 
                 + quat[2] * quat[2] 
                 + quat[3] * quat[3];
  
  if (norm_sq > 0.0) {
    if (quat[0] < 0.0 || std::abs(norm_sq - 1.0) > 1e-12) {
      double inv_norm = 1.0 / std::sqrt(norm_sq);
      if (quat[0] < 0.0) {
        inv_norm = -inv_norm;
      }
      quat[0] *= inv_norm;
      quat[1] *= inv_norm;
      quat[2] *= inv_norm;
      quat[3] *= inv_norm;
    }
  }
}

static void quatToEuler(double quat[4])
{
  // 1. Calculate minimum required intermediate values (7 ops)
  double q0q2 = 2.0 * quat[0] * quat[2]; // 2 * w * y
  double q0q3 = 2.0 * quat[0] * quat[3]; // 2 * w * z

  double q1q1 = 2.0 * quat[1] * quat[1]; // 2 * x * x
  double q1q2 = 2.0 * quat[1] * quat[2]; // 2 * x * y
  double q1q3 = 2.0 * quat[1] * quat[3]; // 2 * x * z

  double q2q2 = 2.0 * quat[2] * quat[2]; // 2 * y * y
  double q3q3 = 2.0 * quat[3] * quat[3]; // 2 * z * z

  // 2. Extract matrix elements strictly necessary for Pitch (quat[2]) and Yaw (quat[1])
  double r20 = q1q3 - q0q2;
  double r10 = q1q2 + q0q3;
  double r00 = 1.0 - q2q2 - q3q3;

  // 3. Compute Pitch
  double pitch = std::asin(PCCClip(-r20, -1.0, 1.0));

  if (std::abs(std::cos(pitch)) > 1e-6) {
    // Regular case (No Gimbal Lock)
    double q0q1 = 2.0 * quat[0] * quat[1]; // 2 * w * x
    double q2q3 = 2.0 * quat[2] * quat[3]; // 2 * y * z
    
    double r21 = q2q3 + q0q1;
    double r22 = 1.0 - q1q1 - q2q2;

    quat[1] = std::atan2(r10, r00);
    quat[3] = std::atan2(r21, r22);
  } else {
    // Gimbal Lock case
    double r01 = q1q2 - q0q3;
    double r11 = 1.0 - q1q1 - q3q3;

    quat[1] = std::atan2(-r01, r11);
    quat[3] = 0.0;
  }
  quat[2] = pitch;
}

const double radToDeg = 180.0 / M_PI;
const double degToRad = M_PI / 180.0;
static void eulerToQuat(double rot[4])
{
  // 1. Extract Euler angles before modifying the array
  // Assuming the sequence matches the previous quatToEuler (Z-Y-X / Yaw-Pitch-Roll)
  double yaw   = rot[1]; // Z-axis rotation
  double pitch = rot[2]; // Y-axis rotation
  double roll  = rot[3]; // X-axis rotation

  // 2. Precompute half angles to eliminate repetitive divisions
  double cy = std::cos(yaw * 0.5);
  double sy = std::sin(yaw * 0.5);
  double cp = std::cos(pitch * 0.5);
  double sp = std::sin(pitch * 0.5);
  double cr = std::cos(roll * 0.5);
  double sr = std::sin(roll * 0.5);

  // 3. Compute Quaternion elements with proper combination (Z-Y-X sequence)
  // Store results directly into the array without cross-contamination
  rot[0] = cy * cp * cr + sy * sp * sr; // w
  rot[1] = cy * cp * sr - sy * sp * cr; // x
  rot[2] = cy * sp * cr + sy * cp * sr; // y
  rot[3] = sy * cp * cr - cy * sp * sr; // z
}

const AttributeDescription& getAttributeDescription(
  const std::vector<AttributeDescription>& attrDescs,
  KnownAttributeLabel attributeLabel, int index)
{
  for(auto& desc: attrDescs) {
    if (desc.attributeLabel == attributeLabel)
    {
      if(desc.attributeLabel == KnownAttributeLabel::kSH_coefficients) {
        if (int(desc.params.sh_coefficient_label) == index) return desc;
      }
      else if (desc.attributeLabel == KnownAttributeLabel::kRotation) {
        if (index == 0 && desc.attr_num_dimensions_minus1 == 2) return desc;
        else if (index == 1 && desc.attr_num_dimensions_minus1 == 0) return desc;
      }
      else {
        return desc;
      }
    }
  }
  assert(0);
  return attrDescs[0];
}

AttributeDescription& getAttributeDescription(
  std::vector<AttributeDescription>& attrDescs,
  KnownAttributeLabel attributeLabel, int index)
{
  for(auto& desc: attrDescs) {
    if (desc.attributeLabel == attributeLabel)
    {
      if(desc.attributeLabel == KnownAttributeLabel::kSH_coefficients) {
        if (int(desc.params.sh_coefficient_label) == index) return desc;
      }
      else if (desc.attributeLabel == KnownAttributeLabel::kRotation) {
        if (index == 0 && desc.attr_num_dimensions_minus1 == 2) return desc;
        else if (index == 1 && desc.attr_num_dimensions_minus1 == 0) return desc;
      }
      else {
        return desc;
      }
    }
  }
  assert(0);
  return attrDescs[0];
}

inline void dequantShChunk(
  size_t chunkStart, size_t chunkCount, std::vector<GSPoint>& gsChunk,
  AttributeParameters& param_sh,
  const std::vector<AttributeDescription>& attrDescs, int numSHCoeffs,
  bool useShNQ,
  std::array<ShQuantRegion, AttributeParameters::kShDrNumSlots>& shRegions,
  const PCCPointSet3& cloud, bool convertColourspace)
{
  for (size_t i = 0; i < chunkCount; ++i) {
    size_t pt = chunkStart + i;
    GSPoint& g = gsChunk[i];
    for (int m = 0; m < numSHCoeffs; m++) {
      g.fxpScale[m] = SH_SCALE;
      auto& desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, m);
      if (param_sh.enableShDimensionReduction_
          && !desc.params.attr_default_value.empty()) {
        g.shFxp[m][0] = 0;
        g.shFxp[m][1] = 0;
        g.shFxp[m][2] = 0;
        continue;
      }
      auto offset = desc.params.attr_offset;
      auto attr_frac_bits = desc.params.attr_frac_bits;
      auto attr_scale = desc.params.attr_scale_minus1 + 1;
      for (int c = 0; c < 3; c++) {
        int32_t attr = int32_t(cloud.getSHCoefficient(m, pt)[c]);
        if (useShNQ && desc.params.shNumRegionsMinus1_[c] > 0) {
          int idx = m * 3 + c;
          int M = desc.params.shNumRegionsMinus1_[c];
          const auto& reg = shRegions[idx];
          int region_idx = 0;
          while (region_idx < M && attr > reg.tp_int[region_idx]) { region_idx++; }
          int64_t prev_idx_fxp = (region_idx > 0) ? (reg.tp_int[region_idx - 1] << SH_SCALE) : reg.min_offset_fxp;
          int64_t step_vec_fxp = desc.params.shStepSizesFxp_[c][region_idx];
          int64_t base_real_fxp = reg.region_start_real_fxp[region_idx] << SS_SCALE;
          int64_t attr_fxp = attr << SH_SCALE;
          int64_t val_fxp = base_real_fxp + (attr_fxp - prev_idx_fxp) * step_vec_fxp;
          g.shFxp[m][c] = val_fxp >> SS_SCALE;
        } else {
          int64_t val = attr * attr_scale + (offset << attr_frac_bits);
          g.shFxp[m][c] = adjust_frac_bits(val, attr_frac_bits, SH_SCALE);
        }
      }
    }
  }
  sphericalHarmonicsInvTransform(true, gsChunk, chunkCount, param_sh, convertColourspace);
}

bool
ply::write(
  Stopwatch* clock,
  const PCCPointSet3& cloud,
  const PropertyNameMap& attributeNames,
  double positionScale,
  Vec3<double> positionOffset,
  const std::string& fileName,
  bool asAscii,
  const std::vector<AttributeDescription>& attrDescs,
  bool asConformance,
  bool convertColourspace,
  bool externalPrePost=true)
{
  clock->stop();
  std::ofstream fout(fileName, std::ofstream::out);
  if (!fout.is_open()) {
    return false;
  }

  const size_t pointCount = cloud.getPointCount();
  fout << "ply" << std::endl;

  if (asAscii) {
    fout << "format ascii 1.0" << std::endl;
  } else {
    PCCEndianness endianess = PCCSystemEndianness();
    if (endianess == PCC_BIG_ENDIAN) {
      fout << "format binary_big_endian 1.0" << std::endl;
    } else {
      fout << "format binary_little_endian 1.0" << std::endl;
    }
  }
  fout << "element vertex " << pointCount << std::endl;
  if (asAscii) {
    fout << "property float " << attributeNames.position[0] << std::endl;
    fout << "property float " << attributeNames.position[1] << std::endl;
    fout << "property float " << attributeNames.position[2] << std::endl;
  } else {
    if (!asConformance) {
      fout << "property float " << attributeNames.position[0] << std::endl;
      fout << "property float " << attributeNames.position[1] << std::endl;
      fout << "property float " << attributeNames.position[2] << std::endl;
    }
    else {
      fout << "property float64 " << attributeNames.position[0] << std::endl;
      fout << "property float64 " << attributeNames.position[1] << std::endl;
      fout << "property float64 " << attributeNames.position[2] << std::endl;
    }
  }

  if (cloud.hasColors()) {
    fout << "property uchar green" << std::endl;
    fout << "property uchar blue" << std::endl;
    fout << "property uchar red" << std::endl;
  }
  if (cloud.hasReflectances()) {
    fout << "property uint16 refc" << std::endl;
  }
  if (cloud.hasFrameIndex()) {
    fout << "property uint8 frameindex" << std::endl;
  }
  if (cloud.hasSHCoefficients()) {
	  fout << "property float f_dc_0" << std::endl;
	  fout << "property float f_dc_1" << std::endl;
	  fout << "property float f_dc_2" << std::endl;
    const int numAcCoeffs = cloud.getNumSHCoeffs() - 1;
	  for (int i = 0; i < numAcCoeffs; i++)
		  fout << "property float f_rest_" << i << std::endl;

	  for (int i = 0; i < numAcCoeffs; i++)
		  fout << "property float f_rest_" << numAcCoeffs + i << std::endl;

	  for (int i = 0; i < numAcCoeffs; i++)
		  fout << "property float f_rest_" << 2 * numAcCoeffs + i << std::endl;
  }
  if (cloud.hasOpacities()) {
	  fout << "property float opacity" << std::endl;
  }
  if (cloud.hasScales()) {
	  fout << "property float scale_0" << std::endl;
	  fout << "property float scale_1" << std::endl;
	  fout << "property float scale_2" << std::endl;
  }
  if (cloud.hasScalarRotations() || cloud.hasRotations()) {
	  fout << "property float rot_0" << std::endl;
	  fout << "property float rot_1" << std::endl;
	  fout << "property float rot_2" << std::endl;
	  fout << "property float rot_3" << std::endl;
  }
  if (asConformance) {
    fout << "element face 0" << std::endl;
    fout << "property list uint8 int32 vertex_index" << std::endl;
  }
  fout << "end_header" << std::endl;
  clock->start();

  if (asAscii) {
    clock->stop();
    //      fout << std::setprecision(std::numeric_limits<double>::max_digits10);
    fout << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < pointCount; ++i) {
      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      fout << position.x() << " " << position.y() << " " << position.z();
      if (cloud.hasColors()) {
        const Vec3<attr_t>& color = cloud.getColor(i);
        fout << " " << static_cast<int>(color[0]) << " "
             << static_cast<int>(color[1]) << " "
             << static_cast<int>(color[2]);
      }
      if (cloud.hasReflectances()) {
        fout << " " << static_cast<int>(cloud.getReflectance(i));
      }
      if (cloud.hasFrameIndex()) {
        fout << " " << static_cast<int>(cloud.getFrameIndex(i));
      }
			if (cloud.hasSHCoefficients()) {
				const Vec3<attr_t>& SHcoefficient_dc = cloud.getSHCoefficient(0, i);
				fout << " " << SHcoefficient_dc[0] 
					<< " "  << SHcoefficient_dc[1] 
					<< " "  << SHcoefficient_dc[2];

				for (int k = 0; k < 3; k++) {
					for (int m = 1; m < cloud.getNumSHCoeffs(); m++) {
						const Vec3<int32_t>& SHcoefficient_ac = cloud.getSHCoefficient(m, i);
						fout << " " << SHcoefficient_ac[k];
					}
				}
			}
			if (cloud.hasOpacities()) {
				fout << " " << cloud.getOpacity(i);
			}
			if (cloud.hasScales()) {
				const Vec3<attr_t>& scale = cloud.getScale(i);
				fout << " " << scale[0]
					<< " " << scale[1]
					<< " " << scale[2];
			}
			if (cloud.hasScalarRotations() || cloud.hasRotations()) {
				const pcc::Vec3<attr_t>& rotation_xyz = cloud.getRotationVec(i);
        fout << " " << cloud.getRotationScalar(i)
					<< " " << rotation_xyz[0]
					<< " " << rotation_xyz[1]
					<< " " << rotation_xyz[2];
			}
      fout << std::endl;
      clock->start();
    }
  } else {
    fout.clear();
    fout.close();
    fout.open(fileName, std::ofstream::binary | std::ofstream::app);

    bool isRealRotCoded = false;
    bool isEulerRot = false;
    double rotScaleValue = 1.0;
    for (const auto& desc: attrDescs) {
      if (desc.attributeLabel == KnownAttributeLabel::kRotation) {
        if (desc.attr_num_dimensions_minus1 == 1) {
          isRealRotCoded = true;
        } else if (desc.params.rotation_label == RotationLabel::kEuler) {
          isEulerRot = true;
          if (desc.params.euler_angle_unit_type == EulerAngleUnitType::kDegree) {
            rotScaleValue = degToRad;
          }
        }
      }
    }
    

    // PCA/NQ: inverse transform (streaming via chunked GSPoint buffer)
    bool useGsdTransform2d = false;
    bool useShNQ = false;
    constexpr size_t kShChunkSize = 10000;
    int numSHCoeffs = 0;
    std::vector<GSPoint> gsChunk;
    AttributeParameters param_sh{};
    std::array<ShQuantRegion, AttributeParameters::kShDrNumSlots> shRegions;
    if (cloud.hasSHCoefficients()) {
      const AttributeDescription& desc_sh = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, 0);
      param_sh = desc_sh.params;
      useGsdTransform2d = param_sh.attr_sh_grp_transform_enabled;
      useShNQ = param_sh.attr_sh_non_uniform_quant_enabled;

      for (auto& desc : attrDescs) {
        if (desc.attributeLabel == KnownAttributeLabel::kSH_coefficients)
          numSHCoeffs++;
      }
      const bool shAcScalingEnabled = desc_sh.params.attr_ac_scaling_enabled_flag;
      if (useShNQ) {
        for (int l = 0; l < numSHCoeffs; ++l) {
          auto& params = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, l).params;
          const uint32_t bitDepth_sh =
            (shAcScalingEnabled && l > 0) ? params.attr_frac_bits : desc_sh.params.attr_frac_bits;
          int64_t minLim = params.attr_offset;
          int64_t maxLim = params.attr_offset + params.attr_scale_minus1 + 1;
          for (int c = 0; c < 3; ++c) {
            int idx = l * 3 + c;
            int M = params.shNumRegionsMinus1_[c];

            if (M > 0) {
              int64_t d_code_fxp = ((maxLim - minLim) * (1 << bitDepth_sh) + ((1ull << params.attr_frac_bits) >> 1)) / (1ull << params.attr_frac_bits);
              shRegions[idx].tp_int.resize(M);
              shRegions[idx].region_start_real_fxp.resize(M + 1);
              shRegions[idx].region_start_real_fxp[0] = minLim << SH_SCALE;
              for (int m = 0; m < M; ++m) {
                int tp_val = params.shTransitionPoints_[c][m];
                shRegions[idx].tp_int[m] = tp_val;
                int64_t buf = (tp_val * d_code_fxp) + (minLim << bitDepth_sh);
                shRegions[idx].region_start_real_fxp[m + 1] =
                    adjust_frac_bits(buf, bitDepth_sh, SH_SCALE);
              }
              shRegions[idx].min_offset_fxp = params.min_offset_fxp[c];
            }
          }
        }
      }
      if (useGsdTransform2d || useShNQ) {
        gsChunk.resize(kShChunkSize);
      }
    }

    constexpr size_t kScalePowerWriteChunkSize = 65536;
    const AttributeDescription* scaleDesc = nullptr;
    std::vector<Vec3<float>> scalePowerWriteValues;
    if (cloud.hasScales()) {
      scaleDesc = &getAttributeDescription(attrDescs, KnownAttributeLabel::kScale, 0);
      if (scaleDesc->params.scaleParameterPresent)
        scalePowerWriteValues.resize(
          std::min(kScalePowerWriteChunkSize, pointCount));
    }

    for (size_t i = 0; i < pointCount; ++i) {
      if ((useGsdTransform2d || useShNQ) && (i % kShChunkSize == 0)) {
        size_t chunkCount = std::min(kShChunkSize, pointCount - i);
        dequantShChunk(i, chunkCount, gsChunk, param_sh, attrDescs,
                       numSHCoeffs, useShNQ, shRegions, cloud,
                       convertColourspace);
      }

      if (scaleDesc && scaleDesc->params.scaleParameterPresent
          && i % kScalePowerWriteChunkSize == 0) {
        size_t chunkCount =
          std::min(kScalePowerWriteChunkSize, pointCount - i);
        inverseScalePowerTransform(
          *scaleDesc, cloud.getScale().data() + i, scalePowerWriteValues.data(),
          int(chunkCount));
      }

      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      if (!asConformance) {
        Vec3<float> positionFloat;
        positionFloat[0] = static_cast<float>(position[0]);
        positionFloat[1] = static_cast<float>(position[1]);
        positionFloat[2] = static_cast<float>(position[2]);
        clock->stop();
        fout.write(
          reinterpret_cast<const char* const>(&positionFloat), sizeof(float) * 3);
        clock->start();
      } else {
        clock->stop();
        fout.write(
          reinterpret_cast<const char* const>(&position), sizeof(double) * 3);
        clock->start();
      }

      if (cloud.hasColors()) {
        const Vec3<attr_t>& c = cloud.getColor(i);
        Vec3<uint8_t> val8b{uint8_t(c[0]), uint8_t(c[1]), uint8_t(c[2])};
        clock->stop();
        fout.write(reinterpret_cast<const char*>(&val8b), sizeof(uint8_t) * 3);
        clock->start();
      }
      if (cloud.hasReflectances()) {
        const attr_t& reflectance = cloud.getReflectance(i);
        clock->stop();
        fout.write(
          reinterpret_cast<const char*>(&reflectance), sizeof(uint16_t));
        clock->start();
      }
      if (cloud.hasFrameIndex()) {
        const uint8_t findex = uint8_t(cloud.getFrameIndex(i));
        clock->stop();
        fout.write(reinterpret_cast<const char*>(&findex), sizeof(uint8_t));
        clock->start();
      }

      //for GS descs
			if (cloud.hasSHCoefficients()) {
				if (useGsdTransform2d || useShNQ) {
					const size_t ci = i % kShChunkSize;
					GSPoint& g = gsChunk[ci];
					const float shScale = float(1 << SH_SCALE) * GS_PROPERTY_SCALE;
					Vec3<float> val32b_dc{ float(g.shFxp[0][0])/shScale, float(g.shFxp[0][1])/shScale, float(g.shFxp[0][2])/shScale };
          clock->stop();
					fout.write(reinterpret_cast<char*>(&val32b_dc), sizeof(float) * 3);
          clock->start();
					for (int k = 0; k < 3; k++) {
						for (int m = 1; m < cloud.getNumSHCoeffs(); m++) {
							float val32b_ac = float(g.shFxp[m][k])/shScale;
              clock->stop();
							fout.write(reinterpret_cast<char*>(&val32b_ac), sizeof(float));
              clock->start();
						}
					}
				} else {
          const int numShCoeffs = cloud.getNumSHCoeffs();
          double sh_buffer[numShCoeffs][3];
          for (int i = 0; i < numShCoeffs; i++)
            for (int j = 0; j < 3; j++)
              sh_buffer[i][j] = 0.0;
              
          auto& shAttrIds = cloud.getAttrIds(KnownAttributeLabel::kSH_coefficients);
          // Post-processing
          for (auto& attrId : shAttrIds)
          {
            auto& desc = attrDescs[attrId];
            int offset = desc.params.attr_offset;
            int attr_frac_bits = desc.params.attr_frac_bits;
            int attr_scale = desc.params.attr_scale_minus1 + 1;
            int numCh = desc.attr_num_dimensions_minus1 + 1;

            int shCoeffIdx = int(desc.params.sh_coefficient_label);
            const Vec3<attr_t>& attr = cloud.getSHCoefficient(shCoeffIdx, i);
            for (int c = 0 ; c < numCh ; c++ ) {
              sh_buffer[shCoeffIdx][c] = double(attr[c]) * attr_scale / (1ll<<attr_frac_bits) + offset;
            }
            if (!externalPrePost && convertColourspace){
              convertToRgb(
                desc.params.cicp_matrix_coefficients_idx, sh_buffer[shCoeffIdx]);
            }
          }

          // write SHDC
          Vec3<float> val32b_dc{ float(sh_buffer[0][0])/GS_PROPERTY_SCALE, float(sh_buffer[0][1])/GS_PROPERTY_SCALE, float(sh_buffer[0][2])/GS_PROPERTY_SCALE };
          clock->stop();
          fout.write(reinterpret_cast<char*>(&val32b_dc), sizeof(float) * 3);
          clock->start();

          // write SHAC 
          for (int k = 0; k < 3; k++) {
						for (int m = 1; m < numShCoeffs; m++) {
							float val32b_ac = float(sh_buffer[m][k])/GS_PROPERTY_SCALE;
              clock->stop();
							fout.write(reinterpret_cast<char*>(&val32b_ac), sizeof(float));
              clock->start();
						}
					}
				}
			}
			if (cloud.hasOpacities()) {
        auto desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kOpacity, 0);
        auto offset = desc.params.attr_offset;
        auto attr_frac_bits = desc.params.attr_frac_bits;
        auto attr_scale     = desc.params.attr_scale_minus1 + 1;
				const int32_t& opacity = cloud.getOpacity(i);
        double _opacity = double(opacity) * attr_scale / (1<<attr_frac_bits) + offset;
        // If the encoder applied the logit->linear transform (convertOpacity),
        // invert it here; otherwise the dequantized value is already in the
        // raw (logit) domain expected by 3DGS PLY files.
        double _opacityOut = desc.params.opacityDomainTransformEnabled
          ? opacityLinearToLogit(_opacity)
          : _opacity;
				float buf = (float)_opacityOut/GS_PROPERTY_SCALE;
        clock->stop();
				fout.write(reinterpret_cast<char*>(&buf), sizeof(float));
        clock->start();
			}
      if (cloud.hasScales()) {
        float scale[3];
        const auto& desc = *scaleDesc;
        auto offset = desc.params.attr_offset;
        auto attr_frac_bits = desc.params.attr_frac_bits;
        auto attr_scale     = desc.params.attr_scale_minus1 + 1;
        if (desc.params.scaleParameterPresent) {
          const auto& value = scalePowerWriteValues[i % kScalePowerWriteChunkSize];
          scale[0] = value[0];
          scale[1] = value[1];
          scale[2] = value[2];
        } else {
          const auto& attr = cloud.getScale()[i];
          for (int k = 0; k < 3; k++) {
            double buf = double(int32_t(attr[k])) * attr_scale
              / (1 << attr_frac_bits) + offset;
            scale[k] = buf;
          }
        }
				Vec3<float> val32b{ float(scale[0])/GS_PROPERTY_SCALE, float(scale[1])/GS_PROPERTY_SCALE, float(scale[2])/GS_PROPERTY_SCALE };
        clock->stop();
				fout.write(reinterpret_cast<char*>(&val32b), sizeof(float) * 3);
        clock->start();
			}
			if (cloud.hasScalarRotations() || cloud.hasRotations()) {
        double rot[4] = {1.0, 0.0, 0.0, 0.0};
        {
          auto desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kRotation, 0);
          auto offset = desc.params.attr_offset;
          auto attr_frac_bits = desc.params.attr_frac_bits;
          auto attr_scale     = desc.params.attr_scale_minus1 + 1;
          auto attr = cloud.getRotationVec(i);
          for (int k = 0; k < 3; k++) {
            double buf = double(attr[k]) * attr_scale / (1<<attr_frac_bits) + offset;
            rot[k+1]= buf;
          }
        }
        
        if (externalPrePost || isRealRotCoded)
        {
          auto desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kRotation, 1);
          auto offset = desc.params.attr_offset;
          auto attr_frac_bits = desc.params.attr_frac_bits;
          auto attr_scale     = desc.params.attr_scale_minus1 + 1;
          rot[0] = double(int32_t(cloud.getRotationScalar(i))) * attr_scale / (1<<attr_frac_bits) + offset;
        }
        else if (!isEulerRot)// Estimate the real part of the quaternion
        { 
          rot[0] = std::sqrt(std::min(std::max(
            1.0 - rot[1] * rot[1] - rot[2] * rot[2] - rot[3] * rot[3],
            0.0), 1.0));
        }
        
        if (isEulerRot) {
          rot[1] *= rotScaleValue;
          rot[2] *= rotScaleValue;
          rot[3] *= rotScaleValue;
          eulerToQuat(rot);
        }

				const auto rotation_w = rot[0];
				float _buf = (float)rotation_w/GS_PROPERTY_SCALE;
        clock->stop();
				fout.write(reinterpret_cast<char*>(&_buf), sizeof(float));
        clock->start();

				Vec3<float> val32b{ float(rot[1])/GS_PROPERTY_SCALE, float(rot[2])/GS_PROPERTY_SCALE, float(rot[3])/GS_PROPERTY_SCALE };
        clock->stop();
				fout.write(reinterpret_cast<char*>(&val32b), sizeof(float) * 3);
        clock->start();
			}
    }
  }
  fout.close();
  return true;
}

//============================================================================

bool ply::stringToAttributeLabelAndChannel(const std::string& propName, KnownAttributeLabel& kAttrLabel, int& channelIdx, int& maxChannelIdx) 
{
  maxChannelIdx = 3;
  if (propName == "x")
    channelIdx = 0;
  else if (propName == "y")
    channelIdx = 1;
  else if (propName == "z")
    channelIdx = 2;
  else if (propName == "red") {
    kAttrLabel = KnownAttributeLabel::kColour;
    channelIdx = 0;
  } 
  else if (propName == "green") {
    kAttrLabel = KnownAttributeLabel::kColour;
    channelIdx = 1;
  } 
  else if (propName == "blue") {
    kAttrLabel = KnownAttributeLabel::kColour;
    channelIdx = 2;
  } 
  else if (propName == "reflectance") {
    kAttrLabel = KnownAttributeLabel::kReflectance;
    channelIdx = 0;
    maxChannelIdx = 1;
  } 
  else if (propName == "opacity") {
    kAttrLabel = KnownAttributeLabel::kOpacity;
    channelIdx = 0;
    maxChannelIdx = 1;
  } 
  else if (propName == "frame_index" || propName == "frameindex") {
    kAttrLabel = KnownAttributeLabel::kFrameIndex;
    channelIdx = 0;
    maxChannelIdx = 1;
  } 
  else if (propName == "frame_number") {
    kAttrLabel = KnownAttributeLabel::kFrameNumber;
    channelIdx = 0;
    maxChannelIdx = 1;
  } 
  else if (propName == "material_id") {
    kAttrLabel = KnownAttributeLabel::kMaterialId;
    channelIdx = 0;
    maxChannelIdx = 1;
  } 
  else if (propName == "nx") {
    kAttrLabel = KnownAttributeLabel::kNormal;
    channelIdx = 0;
  } 
  else if (propName == "ny") {
    kAttrLabel = KnownAttributeLabel::kNormal;
    channelIdx = 1;
  } 
  else if (propName == "nz") {
    kAttrLabel = KnownAttributeLabel::kNormal;
    channelIdx = 2;
  } 
  else if (propName.find("scale_") == 0) {
    kAttrLabel = KnownAttributeLabel::kScale;
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1));
  }
  else if (propName.find("rot_") == 0) {
    kAttrLabel = KnownAttributeLabel::kRotation;
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1));
  }
  else if (propName.rfind("f_dc_", 0) == 0) {
    kAttrLabel = KnownAttributeLabel::kSH_coefficients;
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1));
  } else if (propName.rfind("f_rest_", 0) == 0) {
    kAttrLabel = KnownAttributeLabel::kSH_coefficients; 
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1)) + 3;
  } 
  else {
    kAttrLabel = KnownAttributeLabel::kOid;
    channelIdx = 0;
    maxChannelIdx = 1;
  }

  return true;
}

bool ply::parseHeader(
  std::ifstream& ifs, PCCPointSet3& cloud, 
  Format& format, std::vector<PropertyInfo>& properties)
{
  std::string line, token;
  properties.clear();
  bool inVertexElement = false;
  int vertexCount = 0;
  while (std::getline(ifs, line)) {
    std::istringstream iss(line);
    iss >> token;

    if (token == "format") {
      std::string formatType;
      iss >> formatType;
      format = (formatType == "ascii") ? Format::ASCII : Format::BINARY;
    } 
    else if (token == "element") {
      std::string type;
      iss >> type;
      inVertexElement = (type == "vertex");
      if (inVertexElement) {
        iss >> vertexCount;
      }
    }
    else if (token == "property" && inVertexElement) {
      std::string typeStr, nameStr;
      iss >> typeStr >> nameStr;
      
      KnownAttributeLabel attrLabel;
      int channel, maxChannel;
      if (!stringToAttributeLabelAndChannel(nameStr, attrLabel, channel, maxChannel)) {
        return false;
      }
      properties.push_back({nameStr, attrLabel, channel, maxChannel, stringToPlyType(typeStr)});
    }
    else if (token == "end_header") {
      cloud.resize(vertexCount);
      break;
      
    }
  }

  // convert SH index from exteral to internal
  if (cloud.hasSHCoefficients())
  {
    int maxShIdx = -1;
    for (const auto& p : properties) {
      if (p.attrLabel == KnownAttributeLabel::kSH_coefficients)
        maxShIdx = std::max(maxShIdx, p.channel);
    }

    if (maxShIdx >= 0) {
      int numBands = 1;
      while (3 * numBands * numBands - 1 < maxShIdx)
        ++numBands;

      const int maxNumShCoeffs = numBands * numBands;
      const int numAcCoeffs = maxNumShCoeffs - 1;

      for (auto& p : properties) {
        if (p.attrLabel != KnownAttributeLabel::kSH_coefficients)
          continue;

        p.maxChannel = maxNumShCoeffs;

        if (p.channel < 3) {
          p.channel *= maxNumShCoeffs;
        } else {
          const int acIdx = p.channel - 3;
          p.channel =
            (acIdx / numAcCoeffs) * maxNumShCoeffs +
            (acIdx % numAcCoeffs) + 1;
        }
      }

      if (maxNumShCoeffs < cloud.getNumSHCoeffs()) {
        std::cerr << "Error: SH Coefficients are not configured." << std::endl;
        return false;
      }
    }
  }

  cloud.resize(vertexCount);
  return vertexCount > 0 && !properties.empty();
}

bool ply::isRequestedAttribute(PropertyInfo& prop, const std::vector<AttributeDescription>& attrDescs) {
  for (auto &desc: attrDescs) {
    if (desc.attributeLabel.known_attribute_label == prop.attrLabel)
    {
      if (desc.attributeLabel.known_attribute_label == KnownAttributeLabel::kSH_coefficients)
      {
        if (int(desc.params.sh_coefficient_label) == (prop.channel % prop.maxChannel))
          return true;
      } else {
        return true;
      }
    }
  }
  return false;
}

bool ply::readVertices(
  std::ifstream& ifs, 
  PCCPointSet3& cloud, const std::vector<AttributeDescription>& attrDescs, double positionScale,
  Format& format, std::vector<PropertyInfo>& properties, std::vector<double>& inputAttrBuf)
{
  int pointCount = cloud.getPointCount();
  int attrIdx = 0;
  for (int i = 0; i < properties.size(); ++i) {
    PropertyInfo& prop = properties[i];
    if (isRequestedAttribute(prop, attrDescs)) {
      prop.attr_idx = attrIdx++;
    } else {
      prop.attr_idx = -1;
    }
  }
  
  inputAttrBuf.resize(attrIdx * pointCount);
  for (size_t i = 0; i < cloud.getPointCount(); ++i) {
    for (size_t j = 0; j < properties.size(); ++j) {
      PropertyInfo& prop = properties[j];
      float value;
      if (format == Format::ASCII) {
        value = readTypedDataAscii(ifs, prop.type);
      } else {
        value = readTypedDataBinary(ifs, prop.type);
      }
      if (prop.name == "x" || prop.name == "y" || prop.name == "z") {
        cloud[i][prop.channel] = value * positionScale;
      } else if (prop.attr_idx > -1) {
        inputAttrBuf[prop.attr_idx * pointCount + i] = value;
      }
    }
  }

  return true;
}

std::vector<int> getAttrIndices(const AttributeDescription& desc, const std::vector<ply::PropertyInfo>& properties)
{
  std::vector<int> indices;
  for (const auto& prop : properties)
  {
    if (desc.attributeLabel != prop.attrLabel) {
      continue;
    }

    int target_idx = -1;

    if (desc.attributeLabel == KnownAttributeLabel::kSH_coefficients)
    {
      if (prop.attr_idx < 0) continue;

      if (int(desc.params.sh_coefficient_label) == (prop.channel % prop.maxChannel)) {
        target_idx = prop.channel / prop.maxChannel;
      }
    } 
    else if (desc.attributeLabel == KnownAttributeLabel::kRotation)
    {
      if (desc.attr_num_dimensions_minus1 == 2) {
        if (prop.channel > 0) target_idx = prop.channel - 1;
      } else {
        if (prop.channel == 0) target_idx = prop.channel;
      }
    } 
    else 
    {
        target_idx = prop.channel;
    }

    if (target_idx >= 0) {
      if (target_idx >= indices.size()) {
        indices.resize(target_idx + 1, -1);
      }
      indices[target_idx] = prop.attr_idx;
    }
  }
  return indices;
}

std::vector<int> getAttrIndices(
  AttributeLabel attributeLabel, 
  const std::vector<ply::PropertyInfo>& properties)
{
    std::vector<int> indices;
    for (const auto& prop : properties)
    {
        if (attributeLabel != prop.attrLabel) {
          continue;
        }

        if (prop.channel >= indices.size()) {
            indices.resize(prop.channel + 1, -1); 
        }
        
        indices[prop.channel] = prop.attr_idx;
    }
    return indices;
}

bool
ply::convertColor (
  int pointCount,
  std::vector<double>& inputAttrBuf,
  std::vector<AttributeDescription>& attrDescs,
  std::vector<PropertyInfo>& properties)
{
  for (auto& desc : attrDescs)
  {
    if (desc.attr_num_dimensions_minus1 == 2
        && desc.attributeLabel == KnownAttributeLabel::kSH_coefficients) {

      // Skip SH color conversion when 2D transform is enabled — it handles its own column transform
      if (desc.params.attr_sh_grp_transform_enabled)
        continue;

      double* buf[3];
      std::vector<int> indices = getAttrIndices(desc, properties);

      buf[0] = &inputAttrBuf[indices[0] * pointCount];
      buf[1] = &inputAttrBuf[indices[1] * pointCount];
      buf[2] = &inputAttrBuf[indices[2] * pointCount];

      for (int i = 0; i < pointCount; i++)
      {
        double buf_tmp[3];

        buf_tmp[0] = buf[0][i];
        buf_tmp[1] = buf[1][i];
        buf_tmp[2] = buf[2][i];

        convertFromRgb(
          desc.params.cicp_matrix_coefficients_idx, buf_tmp);

        buf[0][i] = buf_tmp[0];
        buf[1][i] = buf_tmp[1];
        buf[2][i] = buf_tmp[2];
      }
    }
  }
  return true;
}

bool
ply::convertRotationAndScale (
  int pointCount,
  std::vector<double>& inputAttrBuf,
  std::vector<AttributeDescription>& attrDescs, 
  std::vector<PropertyInfo>& properties) 
{
  double* bufRot[4];
  std::vector<int> indicesRot = getAttrIndices(KnownAttributeLabel::kRotation, properties);
  if (indicesRot.empty()) return false;

  bufRot[0] = &inputAttrBuf[indicesRot[0] * pointCount];
  bufRot[1] = &inputAttrBuf[indicesRot[1] * pointCount];
  bufRot[2] = &inputAttrBuf[indicesRot[2] * pointCount];
  bufRot[3] = &inputAttrBuf[indicesRot[3] * pointCount];

  double* bufScale[3];
  std::vector<int> indicesScale = getAttrIndices(KnownAttributeLabel::kScale, properties);
  if (indicesScale.empty()) return false;

  bufScale[0] = &inputAttrBuf[indicesScale[0] * pointCount];
  bufScale[1] = &inputAttrBuf[indicesScale[1] * pointCount];
  bufScale[2] = &inputAttrBuf[indicesScale[2] * pointCount];

  bool isRealRotCoded = false;
  bool isEulerRot = false;
  double rotScaleValue = 1.0;
  
  for (auto desc: attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kRotation) {
      if (desc.attr_num_dimensions_minus1 == 1) {
        isRealRotCoded = true;
      } else if (desc.params.rotation_label == RotationLabel::kEuler) {
        isEulerRot = true;
        if (desc.params.euler_angle_unit_type == EulerAngleUnitType::kDegree) {
          rotScaleValue = radToDeg;
        }
      }
    }
  }

  double quat[4];
  float scale[3];
  for (int i = 0; i < pointCount; i++)
  {
    // input
    quat[0] = bufRot[0][i];
    quat[1] = bufRot[1][i];
    quat[2] = bufRot[2][i];
    quat[3] = bufRot[3][i];
    
    scale[0] = bufScale[0][i];
    scale[1] = bufScale[1][i];
    scale[2] = bufScale[2][i];

    // processing
    canonicalize(quat, scale);
    
    // normalization
    normalizeQuat(quat, isRealRotCoded, isEulerRot);
    
    // convert Quat to Euler
    if (isEulerRot) {
      quatToEuler(quat);
      quat[1] *= rotScaleValue;
      quat[2] *= rotScaleValue;
      quat[3] *= rotScaleValue;
    }

    // output
    bufRot[0][i] = quat[0];
    bufRot[1][i] = quat[1];
    bufRot[2][i] = quat[2];
    bufRot[3][i] = quat[3];
    
    bufScale[0][i] = scale[0];
    bufScale[1][i] = scale[1];
    bufScale[2][i] = scale[2];
  }

  return true;
}

bool
ply::convertOpacity (
  int pointCount,
  std::vector<double>& inputAttrBuf,
  std::vector<AttributeDescription>& attrDescs,
  std::vector<PropertyInfo>& properties,
  bool convertLogitToLinear,
  bool enableDeadZone)
{
  auto it = std::find_if(
    attrDescs.begin(), attrDescs.end(),
    [](const AttributeDescription& d) {
      return d.attributeLabel == KnownAttributeLabel::kOpacity;
    });
  if (it == attrDescs.end())
    return false;

  std::vector<int> indices = getAttrIndices(*it, properties);
  if (indices.empty() || indices[0] < 0)
    return false;

  double* buf = &inputAttrBuf[indices[0] * pointCount];

  int numModifiedValues = 0;
  if (convertLogitToLinear) {
    for (int i = 0; i < pointCount; ++i)
      buf[i] = opacityLogitToLinear(buf[i]);

    std::vector<int> scaleIndices =
      getAttrIndices(KnownAttributeLabel::kScale, properties);
    const bool hasScale = scaleIndices.size() >= 3
      && scaleIndices[0] >= 0 && scaleIndices[1] >= 0
      && scaleIndices[2] >= 0;
    const double* scale0 = hasScale
      ? &inputAttrBuf[scaleIndices[0] * pointCount] : nullptr;
    const double* scale1 = hasScale
      ? &inputAttrBuf[scaleIndices[1] * pointCount] : nullptr;
    const double* scale2 = hasScale
      ? &inputAttrBuf[scaleIndices[2] * pointCount] : nullptr;

    const int maxCode = (1 << it->bitdepth) - 1;
    const double baseStep = maxCode > 0
      ? (1.0) / maxCode : 0.0;
    constexpr double kOpacitySplit = 0.5;
    constexpr double kSmallGaussianScale = -8.0;

    if (enableDeadZone && baseStep > 0.0) {
      for (int i = 0; i < pointCount; ++i) {
        const bool transparent = buf[i] < kOpacitySplit;
        const bool small = hasScale && std::max(
          scale0[i], std::max(scale1[i], scale2[i])) < kSmallGaussianScale;
        const int deadZoneMultiplier = transparent
          ? (small ? 128 : 32)
          : (small ? 8 : 1);
        if (buf[i] <= deadZoneMultiplier * baseStep) {
          buf[i] = 0.0;
          ++numModifiedValues;
        }
      }
    }
  }

  Stats opacityStat;
  for (int i = 0; i < pointCount; i++)
    opacityStat.update(buf[i]);

  it->params.opacityDomainTransformEnabled = convertLogitToLinear;

  return true;
}

bool
ply::analyzeMinMax (
  int pointCount,
  std::vector<double>& inputAttrBuf,
  std::vector<AttributeDescription>& attrDescs, 
  std::vector<PropertyInfo>& properties) 
{
  for (int d = 0; d < attrDescs.size(); d++)
  {
    AttributeDescription& desc = attrDescs[d];

    if (desc.attributeLabel == KnownAttributeLabel::kFrameIndex) {
      // Keep frame_index in raw integer domain: disable linear scale/offset.
      desc.params.attr_offset = 0;
      desc.params.attr_scale_minus1 = 0;
      desc.params.attr_frac_bits = 0;
      continue;
    }
    Stats stat;
    bool hasScaleParam = desc.params.scaleParameterPresent;
    
    std::vector<int> indices = getAttrIndices(desc, properties);

    for (int c = 0; c < indices.size(); c++) {
      double* buf = &inputAttrBuf[indices[c] * pointCount];
      for (int i = 0; i < pointCount; i++) {
        double value = buf[i];
        stat.update(value);
      }
    }

    int fmin = static_cast<int>(std::floor(stat.min));
    int cmax = static_cast<int>(std::ceil(stat.max));
    int attr_offset = fmin;

    int64_t attr_frac_bits_scale = 1ll << desc.params.attr_frac_bits;
    int64_t bit_depth_scale = (1ll << desc.bitdepth) - 1ll;
    int64_t d_code = int64_t(cmax) - int64_t(fmin);
    int64_t attr_scale = (d_code * attr_frac_bits_scale + (bit_depth_scale>>1)) / bit_depth_scale;

    if (hasScaleParam) {
      d_code = std::max<int64_t>(1, d_code);
      attr_scale = std::max<int64_t>(
        1,
        (d_code * attr_frac_bits_scale + (bit_depth_scale >> 1))
          / bit_depth_scale);
    }

    desc.params.attr_offset = attr_offset;
    desc.params.attr_scale_minus1 = attr_scale - 1;

    if (hasScaleParam)
      setScalePowerTransformParams(
        desc, inputAttrBuf, indices, pointCount, d_code);
  }
  return true;
}

bool
ply::quantizationLinear (
  PCCPointSet3& cloud,
  std::vector<double>& inputAttrBuf,
  std::vector<AttributeDescription>& attrDescs, 
  std::vector<PropertyInfo>& properties) 
{
  int pointCount = cloud.getPointCount();  
  for (auto& desc : attrDescs)
  {
    bool hasScaleParam = desc.params.scaleParameterPresent;
    if (hasScaleParam)
      ensureScalePowerTransformParams(desc);

    const int64_t offset = desc.params.attr_offset;
    const int64_t attr_scale = desc.params.attr_scale_minus1 + 1;
    const int64_t attr_frac_bits_scale = 1ll << desc.params.attr_frac_bits;
    const int64_t clip_max = (1ll << desc.bitdepth) - 1ll;
    
    if (desc.attr_num_dimensions_minus1 == 0) {
      double* src;
      std::vector<int> indices = getAttrIndices(desc, properties);
      src = &inputAttrBuf[indices[0]*pointCount];

      auto& dst = cloud.getAttrScalar(desc);
      for (int i = 0; i < pointCount; i++) {
        dst[i] = attr_t(PCCClip(
                std::llround((src[i] - offset) * attr_frac_bits_scale / attr_scale),
                0, clip_max));
      }
    } else {
      double* src[3];
      std::vector<int> indices = getAttrIndices(desc, properties);
      src[0] = &inputAttrBuf[indices[0]*pointCount];
      src[1] = &inputAttrBuf[indices[1]*pointCount];
      src[2] = &inputAttrBuf[indices[2]*pointCount];

      auto& dst = cloud.getAttrVec(desc);

      if (hasScaleParam) {
        for (int c = 0; c < 3; c++)
          quantizeScalePowerTransform(desc, src[c], c, clip_max, dst, pointCount);
      } else {
        for (int c = 0; c < 3; c++) {
          for (int i = 0; i < pointCount; i++) {
            dst[i][c] = attr_t(PCCClip(
                  std::llround((src[c][i] - offset) * attr_frac_bits_scale / attr_scale),
                  0, clip_max));
          }
        }
      }
    }
    
    if (!desc.params.input_attr_default_value.empty()) {
      for (int c = 0; c <= desc.attr_num_dimensions_minus1; c++) {
        if (hasScaleParam) {
          desc.params.attr_default_value[c] = attr_t(quantizeScalePowerTransform(
            desc, desc.params.input_attr_default_value[c], c, clip_max));
        } else {
          desc.params.attr_default_value[c] =
            attr_t(PCCClip(
              std::llround((desc.params.input_attr_default_value[c] - offset) * attr_frac_bits_scale / attr_scale),
              0, clip_max));
        }
      }
    }
  }
  return true;
}

bool
ply::read(
  const std::string& fileName,
  const PropertyNameMap& attributeNames,
  double positionScale,
  PCCPointSet3& cloud,
  std::vector<AttributeDescription>& attrDescs,
  bool convertColourspace,
  bool externalPrePost,
  Stopwatch* clock,
  bool convertOpacityLogitToLinear,
  bool enableOpacityDeadZone)
{
  std::ifstream ifs(fileName,std::ifstream::in | std::ios::binary);
  if (!ifs.is_open()) {
      std::cerr << "Error: can not open file (" << fileName << ")" << std::endl;
      return false;
  }

  Format format;
  std::vector<PropertyInfo> properties;

  clock->stop();
  if (!parseHeader(ifs, cloud, format, properties)) {
    return false;
  }

  std::vector<double> inputAttrBuf;
  readVertices(ifs, cloud, attrDescs, positionScale, format, properties, inputAttrBuf);
  clock->start();

  // Check if any SH desc has 2D transform enabled
  bool gsdTransform2dPresentFlag = false;
  bool enableShNonUniformQuantization = false;
  for (auto& desc : attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kSH_coefficients) {
      if (desc.params.attr_sh_grp_transform_enabled) gsdTransform2dPresentFlag = true;
      if (desc.params.attr_sh_non_uniform_quant_enabled) enableShNonUniformQuantization = true;
    }
  }

  if (!externalPrePost)
  {
    // Color conversion (skips SH when gsdTransform2dPresentFlag)
    if (convertColourspace) {
      convertColor(cloud.getPointCount(), inputAttrBuf, attrDescs, properties);
    }

    // Rotation and Scale conversion
    convertRotationAndScale(cloud.getPointCount(), inputAttrBuf, attrDescs, properties);

    // Opacity: logit (raw 3DGS PLY domain) -> linear, prior to quantization.
    // The inverse (linear -> logit) is applied in ply::write.
    convertOpacity(
      cloud.getPointCount(), inputAttrBuf, attrDescs, properties,
      convertOpacityLogitToLinear, enableOpacityDeadZone);

    // Analyze Min-Max
    analyzeMinMax(cloud.getPointCount(), inputAttrBuf, attrDescs, properties);
  }

  if (gsdTransform2dPresentFlag || enableShNonUniformQuantization) {
    int pointCount = cloud.getPointCount();

    int numSHCoeffs = 0;
    for (auto& desc : attrDescs) {
      if (desc.attributeLabel == KnownAttributeLabel::kSH_coefficients)
        numSHCoeffs++;
    }

    bool witSHCoefficients = (numSHCoeffs > 0);
    auto& shParams = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, 0).params;

    std::vector<double*> shBuffers(numSHCoeffs * 3);
    std::vector<std::vector<int>> shIndices(numSHCoeffs);
    for (int l = 0; l < numSHCoeffs; ++l) {
      auto& desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, l);
      shIndices[l] = getAttrIndices(desc, properties);
      for (int c = 0; c < 3; ++c) {
        shBuffers[l * 3 + c] = &inputAttrBuf[shIndices[l][c] * pointCount];
      }
    }

    SHBufferView shView(shBuffers.data(), numSHCoeffs, pointCount);

    if (gsdTransform2dPresentFlag) {
      bool useSharedAc;
      if (shParams.shNqRangeMode_ == 1) {
        useSharedAc = true;
      } else if (shParams.shNqRangeMode_ == 2) {
        useSharedAc = false;
      } else {
        double top10corr = computeShAcTop10Correlation(shView, pointCount);
        useSharedAc = (top10corr > 0.5);
      }

      initializeShTransformParams(shParams);
      sphericalHarmonicsTransform(witSHCoefficients, shView, pointCount, shParams, convertColourspace);
      computeShDimensionReductionParams(shParams);
      shParams.shAcRangeShared_ = useSharedAc;
      applyShDimensionReduction(shView, pointCount, shParams);

      if (useSharedAc) {
        // Shared AC + lround: one DC range (band 0) + one AC range (bands 1-15)
        Stats dcStat, acStat;
        for (auto& desc : attrDescs) {
          if (desc.attributeLabel != KnownAttributeLabel::kSH_coefficients)
            continue;
          Stats* target = (int(desc.params.sh_coefficient_label) == 0) ? &dcStat : &acStat;
          std::vector<int> indices = getAttrIndices(desc, properties);
          for (int c = 0; c < indices.size(); c++) {
            double* buf = &inputAttrBuf[indices[c] * pointCount];
            for (int i = 0; i < pointCount; i++) {
              target->update(buf[i]);
            }
          }
        }
        auto applyStats = [](Stats& stat, AttributeParameters& params, int bitdepth) {
          int rounded_min = static_cast<int>(std::lround(stat.min));
          int rounded_max = static_cast<int>(std::lround(stat.max));
          if (rounded_min == 0) rounded_min = -1;
          if (rounded_max == 0) rounded_max = 1;
          int64_t attr_frac_bits_scale = 1ll << params.attr_frac_bits;
          int64_t bit_depth_scale = (1ll << bitdepth) - 1ll;
          int64_t d_code = int64_t(rounded_max) - int64_t(rounded_min);
          int64_t attr_scale = (d_code * attr_frac_bits_scale + (bit_depth_scale>>1)) / bit_depth_scale;
          if (attr_scale <= 0) attr_scale = 1;
          params.attr_offset = rounded_min;
          params.attr_scale_minus1 = attr_scale - 1;
        };
        for (auto& desc : attrDescs) {
          if (desc.attributeLabel != KnownAttributeLabel::kSH_coefficients)
            continue;
          if (int(desc.params.sh_coefficient_label) == 0) {
            applyStats(dcStat, desc.params, desc.bitdepth);
          } else {
            applyStats(acStat, desc.params, desc.bitdepth);
          }
        }
      } else {
        // Per-band + floor/ceil: each band gets its own range
        for (auto& desc : attrDescs) {
          if (desc.attributeLabel != KnownAttributeLabel::kSH_coefficients)
            continue;
          Stats stat;
          std::vector<int> indices = getAttrIndices(desc, properties);
          for (int c = 0; c < indices.size(); c++) {
            double* buf = &inputAttrBuf[indices[c] * pointCount];
            for (int i = 0; i < pointCount; i++) stat.update(buf[i]);
          }
          int rounded_min = static_cast<int>(std::floor(stat.min));
          int rounded_max = static_cast<int>(std::ceil(stat.max));
          int64_t attr_frac_bits_scale = 1ll << desc.params.attr_frac_bits;
          int64_t bit_depth_scale = (1ll << desc.bitdepth) - 1ll;
          int64_t d_code = int64_t(rounded_max) - int64_t(rounded_min);
          int64_t attr_scale = (d_code * attr_frac_bits_scale + (bit_depth_scale>>1)) / bit_depth_scale;
          if (attr_scale <= 0) attr_scale = 1;
          desc.params.attr_offset = rounded_min;
          desc.params.attr_scale_minus1 = attr_scale - 1;
        }
      }
    }

    if (enableShNonUniformQuantization) {
      extractSHQuantRegions(witSHCoefficients, shView, pointCount, attrDescs);

      std::array<ShQuantRegion, AttributeParameters::kShDrNumSlots> shRegions;
      const bool shAcScalingEnabled = shParams.attr_ac_scaling_enabled_flag;
      for (int l = 0; l < numSHCoeffs; ++l) {
        auto& desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, l);
        auto& params = desc.params;
        const uint32_t bitDepth_sh = (shAcScalingEnabled && l > 0) ? params.attr_frac_bits : shParams.attr_frac_bits;
        const uint32_t max_index_sh = 1u << bitDepth_sh;
        for (int c = 0; c < 3; ++c) {
          int idx = l * 3 + c;
          int M = params.shNumRegionsMinus1_[c];
          double minLim = params.attr_offset;
          double maxLim = params.attr_offset + params.attr_scale_minus1 + 1;
          double Range = maxLim - minLim;
          double d_code = std::round(Range * (1 << bitDepth_sh) / (double)(1ull << params.attr_frac_bits));

          if (M > 0) {
            shRegions[idx].tp_int.resize(M);
            shRegions[idx].region_start_real.resize(M + 1);
            shRegions[idx].region_start_real[0] = minLim;

            for (int m = 0; m < M; ++m) {
              int tp_val = params.shTransitionPoints_[c][m];
              shRegions[idx].tp_int[m] = tp_val;
              shRegions[idx].region_start_real[m + 1] = tp_val * d_code / static_cast<double>(max_index_sh) + minLim;
            }
            shRegions[idx].min_offset =
                shRegions[idx].tp_int[0] - (shRegions[idx].region_start_real[1] - minLim) / params.shStepSizes_[c][0];

            const double sh_scale = (1 << SH_SCALE);
            params.min_offset_fxp[c] = (int64_t)round(shRegions[idx].min_offset * sh_scale);
            params.min_offset[c] = shRegions[idx].min_offset;
          }
        }
      }

      const int64_t clip_max = (1ll << 16) - 1;
      for (int l = 0; l < numSHCoeffs; ++l) {
        auto& desc = getAttributeDescription(attrDescs, KnownAttributeLabel::kSH_coefficients, l);
        auto& params = desc.params;
        auto attr_frac_bits_scale = 1ll << params.attr_frac_bits;
        int64_t offset = params.attr_offset;
        int64_t attr_scale = params.attr_scale_minus1 + 1;

        for (int c = 0; c < 3; ++c) {
          double* src = &inputAttrBuf[getAttrIndices(desc, properties)[c] * pointCount];
          auto& dst = cloud.getAttrVec(desc);

          for (int i = 0; i < pointCount; i++) {
            double attr = src[i];
            int idx = l * 3 + c;
            int M = params.shNumRegionsMinus1_[c];
            if (M == 0) {
              dst[i][c] = attr_t(PCCClip(
                  std::llround((attr - offset) * attr_frac_bits_scale / attr_scale),
                  0, clip_max));
            } else {
              const auto& reg = shRegions[idx];
              int region_idx = 0;
              while (region_idx < M && attr > reg.region_start_real[region_idx + 1]) { region_idx++; }
              double prev_idx = (region_idx > 0) ? reg.tp_int[region_idx - 1] : reg.min_offset;
              double step_vec = params.shStepSizes_[c][region_idx];
              double base_real = reg.region_start_real[region_idx];
              double q_vals = std::round((attr - base_real) / step_vec + prev_idx);
              dst[i][c] = attr_t(PCCClip(static_cast<int64_t>(q_vals), int64_t(0), clip_max));
            }
          }
        }
      }
    }
  }

  if (!enableShNonUniformQuantization) {
    quantizationLinear(cloud, inputAttrBuf, attrDescs, properties);
  } else {
    int pointCount = cloud.getPointCount();
    for (auto& desc : attrDescs)
    {
      if (desc.attributeLabel == KnownAttributeLabel::kSH_coefficients)
        continue;  // SH already quantized by NQ above

      bool hasScaleParam = desc.params.scaleParameterPresent;
      if (hasScaleParam)
        ensureScalePowerTransformParams(desc);

      const int64_t offset = desc.params.attr_offset;
      const int64_t attr_scale = desc.params.attr_scale_minus1 + 1;
      const int64_t attr_frac_bits_scale = 1ll << desc.params.attr_frac_bits;
      const int64_t clip_max = (1ll << desc.bitdepth) - 1ll;

      if (desc.attr_num_dimensions_minus1 == 0) {
        double* src;
        std::vector<int> indices = getAttrIndices(desc, properties);
        src = &inputAttrBuf[indices[0]*pointCount];

        auto& dst = cloud.getAttrScalar(desc);
        for (int i = 0; i < pointCount; i++) {
          dst[i] = attr_t(PCCClip(
                  std::llround((src[i] - offset) * attr_frac_bits_scale / attr_scale),
                  0, clip_max));
        }
      } else {
        double* src[3];
        std::vector<int> indices = getAttrIndices(desc, properties);
        src[0] = &inputAttrBuf[indices[0]*pointCount];
        src[1] = &inputAttrBuf[indices[1]*pointCount];
        src[2] = &inputAttrBuf[indices[2]*pointCount];

        auto& dst = cloud.getAttrVec(desc);

        if (hasScaleParam) {
          for (int c = 0; c < 3; c++)
            quantizeScalePowerTransform(desc, src[c], c, clip_max, dst, pointCount);
        } else {
          for (int c = 0; c < 3; c++) {
            for (int i = 0; i < pointCount; i++) {
              dst[i][c] = attr_t(PCCClip(
                    std::llround((src[c][i] - offset) * attr_frac_bits_scale / attr_scale),
                    0, clip_max));
            }
          }
        }
      }

      if (!desc.params.input_attr_default_value.empty()) {
        for (int c = 0; c <= desc.attr_num_dimensions_minus1; c++) {
          if (hasScaleParam) {
            desc.params.attr_default_value[c] = attr_t(quantizeScalePowerTransform(
              desc, desc.params.input_attr_default_value[c], c, clip_max));
          } else {
            desc.params.attr_default_value[c] =
              attr_t(PCCClip(
                std::llround((desc.params.input_attr_default_value[c] - offset) * attr_frac_bits_scale / attr_scale),
                0, clip_max));
          }
        }
      }
    }
  }

  clock->stop();

  return true;
}

//============================================================================

}  // namespace pcc
