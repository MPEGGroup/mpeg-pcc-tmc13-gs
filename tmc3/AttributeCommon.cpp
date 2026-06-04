/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2019, ISO/IEC
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

#include "AttributeCommon.h"
#include "PCCTMC3Common.h"

namespace pcc {

//============================================================================
// AttributeLods methods

void
AttributeLods::generate(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const PCCPointSet3& cloud)
{
  _aps = aps;
  _abh = abh;

  if (minGeomNodeSizeLog2 > 0)
    assert(aps.scalable_lifting_enabled_flag);

  buildPredictorsFast(
    aps, abh, cloud, minGeomNodeSizeLog2, geom_num_points_minus1, predictors,
    numPointsInLod, indexes);

  assert(predictors.size() == cloud.getPointCount());
  for (auto& predictor : predictors) {
    predictor.computeWeights();
    if (aps.attr_encoding == AttributeEncoding::kPredictingTransform)
      if (aps.pred_weight_blending_enabled_flag)
        predictor.blendWeights(cloud, indexes);
  }
}

//----------------------------------------------------------------------------

bool
AttributeLods::isReusable(
  const AttributeParameterSet& aps, const AttributeBrickHeader& abh) const
{
  // No LoDs cached => can be reused by anything
  if (numPointsInLod.empty())
    return true;

  // If the other aps doesn't use LoDs, it is compatible.
  // Otherwise, if both use LoDs, check each parameter
  if (!(_aps.lodParametersPresent() && aps.lodParametersPresent()))
    return true;

  // NB: the following comparison order needs to be the same as the i/o
  // order otherwise comparisons may involve undefined values

  if (
    _aps.num_pred_nearest_neighbours_minus1
    != aps.num_pred_nearest_neighbours_minus1)
    return false;

  if (_aps.inter_lod_search_range != aps.inter_lod_search_range)
    return false;

  if (_aps.intra_lod_search_range != aps.intra_lod_search_range)
    return false;

  if (_aps.num_detail_levels_minus1 != aps.num_detail_levels_minus1)
    return false;

  if (_aps.lodNeighBias != aps.lodNeighBias)
    return false;

  // until this feature is stable, always generate LoDs.
  if (_aps.scalable_lifting_enabled_flag || aps.scalable_lifting_enabled_flag)
    return false;

  if (_aps.lod_decimation_type != aps.lod_decimation_type)
    return false;

  if (_aps.dist2 + _abh.attr_dist2_delta != aps.dist2 + abh.attr_dist2_delta)
    return false;

  if (_aps.lodSamplingPeriod != aps.lodSamplingPeriod)
    return false;

  if (
    _aps.intra_lod_prediction_skip_layers
    != aps.intra_lod_prediction_skip_layers)
    return false;

  if (_aps.canonical_point_order_flag != aps.canonical_point_order_flag)
    return false;

  if (
    _aps.pred_weight_blending_enabled_flag
    != aps.pred_weight_blending_enabled_flag)
    return false;

  return true;
}

//============================================================================

bool
predModeEligibleAttrMultiComp(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexes,
  const PCCPredictor& predictor)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;

  if (predictor.neighborCount <= 1 || !aps.max_num_direct_predictors)
    return false;

  Vec3<int64_t> minValue = {0, 0, 0};
  Vec3<int64_t> maxValue = {0, 0, 0};

  for (int i = 0; i < predictor.neighborCount; ++i) {
    Vec3<attr_t> attributeNeighbor = getAttribute(
      desc, pointCloud, indexes[predictor.neighbors[i].predictorIndex]);
    for (int k = 0; k < attribCount; ++k) {
      if (i == 0 || attributeNeighbor[k] < minValue[k])
        minValue[k] = attributeNeighbor[k];
      if (i == 0 || attributeNeighbor[k] > maxValue[k])
        maxValue[k] = attributeNeighbor[k];
    }
  }

  int64_t maxDiff = 0;
  for (int k = 0; k < attribCount; ++k)
    maxDiff = std::max(maxDiff, maxValue[k] - minValue[k]);
  return maxDiff >= aps.adaptivePredictionThreshold(desc);
}

//----------------------------------------------------------------------------

bool
predModeEligibleAttrSingleComp(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexes,
  const PCCPredictor& predictor)
{
  if (predictor.neighborCount <= 1 || !aps.max_num_direct_predictors)
    return false;

  int64_t minValue = 0;
  int64_t maxValue = 0;
  for (int i = 0; i < predictor.neighborCount; ++i) {
    attr_t attributeNeighbor = getAttribute(
      desc, pointCloud, indexes[predictor.neighbors[i].predictorIndex])[0];
    if (i == 0 || attributeNeighbor < minValue)
      minValue = attributeNeighbor;
    if (i == 0 || attributeNeighbor > maxValue)
      maxValue = attributeNeighbor;
  }

  auto maxDiff = maxValue - minValue;
  return maxDiff >= aps.adaptivePredictionThreshold(desc);
}

//============================================================================

Vec3<attr_t>
getAttribute(
  const AttributeDescription& desc,
  const PCCPointSet3& pointCloud,
  uint32_t pointIndex)
{
  Vec3<attr_t> attribute(0);

  if (desc.attributeLabel == KnownAttributeLabel::kColour) {
    attribute = pointCloud.getColor(pointIndex);
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kScale) {
    attribute = pointCloud.getScale(pointIndex);
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kRotation) {
    if (desc.params.rotation_component_label ==
        RotationComponentLabel::kRotation_XYZ) {
      attribute = pointCloud.getRotationVec(pointIndex);
    }
    else {
      attribute[0] = pointCloud.getRotationScalar(pointIndex);
    }
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kReflectance) {
    attribute[0] = pointCloud.getReflectance(pointIndex);
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kOpacity) {
    attribute[0] = pointCloud.getOpacity(pointIndex);
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kSH_coefficients) {
    int sh = int(desc.params.sh_coefficient_label);
    attribute = pointCloud.getSHCoefficient(sh, pointIndex);
  }

  return attribute;
}

//----------------------------------------------------------------------------

void
setAttribute(
  const AttributeDescription& desc,
  PCCPointSet3& pointCloud,
  uint32_t pointIndex,
  const Vec3<attr_t>& attribute)
{
  if (desc.attributeLabel == KnownAttributeLabel::kColour) {
    pointCloud.setColor(pointIndex, attribute);
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kScale) {
    pointCloud.setScale(pointIndex, attribute);
  }
  else if (desc.attributeLabel == KnownAttributeLabel::kRotation) {
    if (desc.params.rotation_component_label ==
        RotationComponentLabel::kRotation_XYZ)
      pointCloud.setRotationVec(pointIndex, attribute);
    else
      pointCloud.setRotationScalar(pointIndex, attribute[0]);
  }
  else if (desc.attributeLabel ==
           KnownAttributeLabel::kReflectance) {
    pointCloud.setReflectance(pointIndex, attribute[0]);
  }
  else if (desc.attributeLabel ==
           KnownAttributeLabel::kOpacity) {
    pointCloud.setOpacity(pointIndex, attribute[0]);
  }
  else if (desc.attributeLabel ==
           KnownAttributeLabel::kSH_coefficients) {
    int SHindex = int(desc.params.sh_coefficient_label);
    pointCloud.setSHCoefficient(SHindex, pointIndex, attribute);
  }
}

//============================================================================

}  // namespace pcc
