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

#include "AttributeDecoder.h"
#include <set>
#include "AttributeCommon.h"
#include "DualLutCoder.h"
#include "attribute_raw.h"
#include "constants.h"
#include "entropy.h"
#include "hls.h"
#include "io_hls.h"
#include "RAHT.h"
#include "FixedPoint.h"

namespace pcc {

//============================================================================
// An encapsulation of the entropy decoding methods used in attribute coding

class PCCResidualsDecoder : protected AttributeContexts {
public:
  PCCResidualsDecoder(
    const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem);

  EntropyDecoder arithmeticDecoder;

  const AttributeContexts& getCtx() const { return *this; }

  void start(const SequenceParameterSet& sps, const char* buf, int buf_len);
  void stop();

  int decodeRunLength();
  int decodeSymbol(int k1, int k2, int k3);
  void decode(int32_t* values, int attribCount);
};

//----------------------------------------------------------------------------

PCCResidualsDecoder::PCCResidualsDecoder(
  const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem)
  : AttributeContexts(ctxtMem)
{}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::start(
  const SequenceParameterSet& sps, const char* buf, int buf_len)
{
  arithmeticDecoder.setBuffer(buf_len, buf);
  arithmeticDecoder.enableBypassStream(sps.cabac_bypass_stream_enabled_flag);
  arithmeticDecoder.start();
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::stop()
{
  arithmeticDecoder.stop();
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeRunLength()
{
  int runLength = 0;
  auto* ctx = ctxRunLen;
  for (; runLength < 3; runLength++, ctx++) {
    int bin = arithmeticDecoder.decode(*ctx);
    if (!bin)
      return runLength;
  }

  for (int i = 0; i < 4; i++) {
    int bin = arithmeticDecoder.decode(*ctx);
    if (!bin) {
      runLength += arithmeticDecoder.decode();
      return runLength;
    }
    runLength += 2;
  }

  runLength += arithmeticDecoder.decodeExpGolomb(2, *++ctx);
  return runLength;
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeSymbol(int k1, int k2, int k3)
{
  if (!arithmeticDecoder.decode(ctxCoeffGtN[0][k1]))
    return 0;

  if (!arithmeticDecoder.decode(ctxCoeffGtN[1][k2]))
    return 1;

  int coeff_abs_minus2 = arithmeticDecoder.decodeExpGolomb(
    1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);

  return coeff_abs_minus2 + 2;
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::decode(int32_t* values, int attribCount)
{
  if (attribCount == 1) {
    int32_t mag = decodeSymbol(0, 0, 0) + 1;
    int sign = arithmeticDecoder.decode();
    values[0] = sign ? -mag : mag;
    return;
  }

  // Decode component 1 (non-last)
  values[1] = decodeSymbol(0, 0, 1);
  int b0 = (values[1] == 0);
  int b1 = (values[1] <= 1);
  if (values[1] && arithmeticDecoder.decode()) {
    values[1] = -values[1];
  }

  // Decode component 2 (only for 3 components)
  int b2 = 0, b3 = 0;
  if (attribCount == 3) {
    values[2] = decodeSymbol(1 + b0, 1 + b1, 1);
    b2 = (values[2] == 0);
    b3 = (values[2] <= 1);
    if (values[2] && arithmeticDecoder.decode()) {
      values[2] = -values[2];
    }
  }

  // Decode component 0 (last), conditioned on previous components
  values[0] = decodeSymbol(3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 0);

  if (attribCount == 2) {
    values[0] += b0;
  } else { // attribCount == 3
    values[0] += (b0 && b2);
  }

  if (values[0] && arithmeticDecoder.decode()) {
    values[0] = -values[0];
  }
}

//============================================================================
// AttributeDecoderIntf

AttributeDecoderIntf::~AttributeDecoderIntf() = default;

//============================================================================
// AttributeDecoder factory

std::unique_ptr<AttributeDecoderIntf>
makeAttributeDecoder()
{
  return std::unique_ptr<AttributeDecoder>(new AttributeDecoder());
}

//============================================================================
// AttributeDecoder Members

void
AttributeDecoder::decode(
  const SequenceParameterSet& sps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const char* payload,
  size_t payloadLen,
  AttributeContexts& ctxtMem,
  PCCPointSet3& pointCloud,
  int rahtMaxAttributeGroupNum,
  std::unordered_map<int, std::vector<int>>* bufferRahtCoefficients)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    AttrRawDecoder::decode(
      attr_desc, attr_aps, abh, payload, payloadLen, pointCloud);
    return;
  }

  QpSet qpSet = deriveQpSet(attr_desc, attr_aps, abh);

  PCCResidualsDecoder decoder(abh, ctxtMem);
  decoder.start(sps, payload, payloadLen);

  // generate LoDs if necessary
  if (attr_aps.lodParametersPresent() && _lods.empty())
    _lods.generate(
      attr_aps, abh, geom_num_points_minus1, minGeomNodeSizeLog2, pointCloud);

  switch (attr_aps.attr_encoding) {
  case AttributeEncoding::kRAHTransform:
    decodeAttributeRaht(attr_desc, attr_aps, qpSet, decoder, pointCloud, rahtMaxAttributeGroupNum, bufferRahtCoefficients);
    break;

  case AttributeEncoding::kPredictingTransform:
    decodeAttributePred(attr_desc, attr_aps, abh, qpSet, decoder, pointCloud);
    break;

  case AttributeEncoding::kLiftingTransform:
    decodeAttributeLift(
      attr_desc, attr_aps, abh, qpSet, geom_num_points_minus1,
      minGeomNodeSizeLog2, decoder, pointCloud);
    break;

  case AttributeEncoding::kRaw:
    // Already handled
    break;
  }

  decoder.stop();

  // save the context state for re-use by a future slice if required
  ctxtMem = decoder.getCtx();
}

//----------------------------------------------------------------------------

bool
AttributeDecoder::isReusable(
  const AttributeParameterSet& aps, const AttributeBrickHeader& abh) const
{
  return _lods.isReusable(aps, abh);
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodePredModeAttrSingleComp(
  const AttributeParameterSet& aps, int32_t& coeff, PCCPredictor& predictor)
{
  int coeffAbs = abs(coeff);
  int coeffSign = coeff < 0 ? -1 : 1;
  int mode;

  int maxcand =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;
  switch (maxcand) {
  case 4:
    mode = coeffAbs & 3;
    coeff = coeffSign * (coeffAbs >> 2);
    break;

  case 3:
    mode = coeffAbs & 1;
    coeffAbs >>= 1;
    if (mode > 0) {
      mode += coeffAbs & 1;
      coeffAbs >>= 1;
    }
    coeff = coeffSign * coeffAbs;
    break;

  case 2:
    mode = coeffAbs & 1;
    coeff = coeffSign * (coeffAbs >> 1);
    break;

  default: mode = 0;
  }

  predictor.predMode = mode + aps.direct_avg_predictor_disabled_flag;
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodePredModeAttrMultiComp(
  const AttributeParameterSet& aps,
  Vec3<int32_t>& coeff,
  PCCPredictor& predictor,
  int attribCount)
{
  int signk1 = coeff[1] < 0 ? -1 : 1;
  int signk2 = coeff[2] < 0 ? -1 : 1;
  int coeffAbsk1 = abs(coeff[1]);
  int coeffAbsk2 = abs(coeff[2]);

  int mode;
  int maxcand =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;

  if (attribCount == 2) {
    mode = coeffAbsk1 % maxcand;
    coeff[1] = signk1 * (coeffAbsk1 / maxcand);
    predictor.predMode = mode + aps.direct_avg_predictor_disabled_flag;
    return;
  }

  switch (maxcand) {
    int parityk1, parityk2;
  case 4:
    parityk1 = coeffAbsk1 & 1;
    parityk2 = coeffAbsk2 & 1;
    coeff[1] = signk1 * (coeffAbsk1 >> 1);
    coeff[2] = signk2 * (coeffAbsk2 >> 1);

    mode = (parityk1 << 1) + parityk2;
    break;

  case 3:
    parityk1 = coeffAbsk1 & 1;
    coeff[1] = signk1 * (coeffAbsk1 >> 1);
    mode = parityk1;
    if (parityk1) {
      parityk2 = coeffAbsk2 & 1;
      coeff[2] = signk2 * (coeffAbsk2 >> 1);
      mode += parityk2;
    }
    break;

  case 2:
    parityk1 = coeffAbsk1 & 1;
    coeff[1] = signk1 * (coeffAbsk1 >> 1);
    mode = parityk1;
    break;

  default: assert(maxcand >= 2); mode = 0;
  }

  predictor.predMode = mode + aps.direct_avg_predictor_disabled_flag;
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeAttributePred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;

  int64_t clipMax = (1 << desc.bitdepth) - 1;
  int32_t values[3] = {0, 0, 0};

  bool icpPresent = abh.icpPresent(desc, aps) && attribCount >= 2;

  int zeroRunRem = 0;
  const size_t lodCount = _lods.numPointsInLod.size();
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    int quantLayer = std::min(int(qpSet.layers.size()) - 1, int(lod));
    auto icpCoeff = icpPresent ? abh.icpCoeffs[lod] : 0;

    auto coeffRange = LodCoeffIdxRange(_lods.numPointsInLod, lod);
    for (auto predictorIndex : coeffRange) {
      const uint32_t pointIndex = _lods.indexes[predictorIndex];
      auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
      auto& predictor = _lods.predictors[predictorIndex];
      predictor.predMode = 0;

      if (--zeroRunRem < 0)
        zeroRunRem = decoder.decodeRunLength();

      if (zeroRunRem)
        values[0] = values[1] = values[2] = 0;
      else
        decoder.decode(values, attribCount);

      if (attribCount == 1) {
        if (predModeEligibleAttrSingleComp(
              desc, aps, pointCloud, _lods.indexes, predictor))
          decodePredModeAttrSingleComp(aps, values[0], predictor);
      } else {
        if (predModeEligibleAttrMultiComp(
              desc, aps, pointCloud, _lods.indexes, predictor))
          decodePredModeAttrMultiComp(aps, *reinterpret_cast<Vec3<int32_t>*>(values), predictor, attribCount);
      }

      Vec3<attr_t> predictedAttribute;
      if (attribCount == 1) {
        predictedAttribute = predictor.predictAttrSingleComp(pointCloud, _lods.indexes, desc);
      } else {
        predictedAttribute =
          predictor.predictAttrMultiComp(pointCloud, _lods.indexes, desc);
      }

      Vec3<attr_t> attribute(0);
      int64_t residual0 = 0;
      for (int k = 0; k < attribCount; ++k) {
        const auto& q = quant[std::min(k, 1)];
        int64_t residual =
          divExp2RoundHalfUp(q.scale(values[k]), kFixedPointAttributeShift);

        // ICP applies to all non-base components (k >= 1)
        // supports 2 and 3 component attributes
        if (icpPresent && k >= 1)
          residual = residual + ((int64_t(icpCoeff[k]) * residual0 + 2) >> 2);

        int64_t recon =
          int64_t(predictedAttribute[k]) + residual;
        attribute[k] = attr_t(PCCClip(recon, int64_t(0), clipMax));
        if (k == 0 && aps.inter_component_prediction_enabled_flag)
          residual0 = residual;
      }
      setAttribute(desc, pointCloud, pointIndex, attribute);
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeAttributeRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  int rahtMaxAttributeGroupNum,
  std::unordered_map<int, std::vector<int>>* bufferRahtCoefficients)
{
  const int voxelCount = int(pointCloud.getPointCount());
  if (_sortedMortonCode.size() == 0) {
    std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
    for (int n = 0; n < voxelCount; n++) {
      packedVoxel[n].mortonCode = mortonAddr(pointCloud[n]);
      packedVoxel[n].index = n;
    }
    sort(packedVoxel.begin(), packedVoxel.end());

    _sortedMortonCode.resize(voxelCount);
    _sortedIndex.resize(voxelCount);
    for (int n = 0; n < voxelCount; n++) {
      _sortedMortonCode[n] = packedVoxel[n].mortonCode;
      _sortedIndex[n] = packedVoxel[n].index;
    }
  }

  // Entropy decode
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  std::vector<int> coefficients(attribCount * voxelCount);
  std::vector<Qps> pointQpOffsets(voxelCount);

  int zeroRunRem = 0;
  for (int n = 0; n < voxelCount; ++n)
    pointQpOffsets[n] = qpSet.regionQpOffset(pointCloud[_sortedIndex[n]]);

  std::vector<int32_t> values(attribCount);
  for (int n = 0; n < voxelCount; ++n) {
    if (--zeroRunRem < 0)
      zeroRunRem = decoder.decodeRunLength();

    if (!zeroRunRem) {
      decoder.decode(values.data(), attribCount);
    } else {
      std::fill(values.begin(), values.end(), 0);
    }

    for (int d = 0; d < attribCount; ++d)
      coefficients[voxelCount * d + n] = values[d];
  }

  if (rahtMaxAttributeGroupNum > 1){
    auto aps_attr_parameter_set_id = aps.aps_attr_parameter_set_id;
    (*bufferRahtCoefficients)[aps_attr_parameter_set_id] = coefficients;
      return;
  }

  std::vector<int> attributes(attribCount * voxelCount);
  const int rahtPredThreshold[2] = {aps.raht_prediction_threshold0,
                                    aps.raht_prediction_threshold1};

  FixedPoint::kFracBits = aps.raht_fixed_point_precision_minus15 + 15;

  init_LUT();

  regionAdaptiveHierarchicalInverseTransform(
    aps.raht_prediction_enabled_flag, rahtPredThreshold, qpSet,
    pointQpOffsets.data(), _sortedMortonCode.data(), attributes.data(), attribCount,
    voxelCount, coefficients.data());

  const int64_t clipMax = (1 << desc.bitdepth) - 1;
  for (int n = 0; n < voxelCount; ++n) {
    const uint32_t sortedIndex = _sortedIndex[n];
    Vec3<attr_t> attribute(0);
    for (int d = 0; d < attribCount; ++d) {
      int value = attributes[attribCount * n + d];
      attribute[d] = attr_t(PCCClip(value, 0, clipMax));
    }

    setAttribute(desc, pointCloud, sortedIndex, attribute);
  }
}

//----------------------------------------------------------------------------
std::vector<std::vector<int>>
AttributeGroup(
  std::vector<AttributeBrickHeader>& allAbh,
  std::map<int, AttributeParameterSet> _apss,
  int Group_maxNum)
{
  std::vector<std::vector<int>> Groups;
  std::vector<int> Group;
  std::set<std::array<int, 2>> uniquePredictionThreshold;

  for (int i = 0; i < allAbh.size();i++) {
    auto abh = allAbh[i];
    int attrIdx = abh.attr_attr_parameter_set_id;
    const auto it_attr_aps = _apss.find(attrIdx);
    const auto& attr_aps = it_attr_aps->second;

    if (attr_aps.raht_prediction_enabled_flag) {
      auto predictionThreshold0 = attr_aps.raht_prediction_threshold0;
      auto predictionThreshold1 = attr_aps.raht_prediction_threshold1;
      uniquePredictionThreshold.insert(
        {predictionThreshold0, predictionThreshold1});
    } else {
      Group.push_back(i);
      if (Group.size() >= Group_maxNum) {
        Groups.push_back(Group);
        Group.clear();
      }
    }
  }
  if (Group.size() != 0)
    Groups.push_back(Group);
  Group.clear();

  for (const auto Threshold : uniquePredictionThreshold) {
    for (int i = 0; i < allAbh.size(); i++) {
      auto abh = allAbh[i];
      int attrIdx = abh.attr_attr_parameter_set_id;
      const auto it_attr_aps = _apss.find(attrIdx);
      const auto& attr_aps = it_attr_aps->second;

      if (attr_aps.raht_prediction_enabled_flag) {
        auto predictionThreshold0 = attr_aps.raht_prediction_threshold0;
        auto predictionThreshold1 = attr_aps.raht_prediction_threshold1;
        if (
          predictionThreshold0 == Threshold[0]
          && predictionThreshold1 == Threshold[1]) {
          Group.push_back(i);
          if (Group.size() >= Group_maxNum) {
            Groups.push_back(Group);
            Group.clear();
          }
        }
      }
    }
    if (Group.size() != 0)
      Groups.push_back(Group);
    Group.clear();
  }

  return Groups;
}

void
AttributeDecoder::decodeMultiAttributeRAHT(
  const SequenceParameterSet& sps,
  std::map<int, AttributeParameterSet> _apss,
  std::vector<AttributeBrickHeader>& _allAttrAbh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  std::unordered_map<int, std::vector<int>>& bufferRahtCoefficients,
  PCCPointSet3& pointCloud,
  int maxAttributeNumberPerGroup)
{
  std::vector<std::vector<int>> AttributeGroups = AttributeGroup(_allAttrAbh, _apss, maxAttributeNumberPerGroup);  //Attribute Set for Multi-Attribute decoding
  const int voxelCount = int(pointCloud.getPointCount());
  if (_sortedMortonCode.size() == 0) {
    std::vector<MortonCodeWithIndex> packedVoxel(voxelCount);
    for (int n = 0; n < voxelCount; n++) {
      packedVoxel[n].mortonCode = mortonAddr(pointCloud[n]);
      packedVoxel[n].index = n;
    }
    sort(packedVoxel.begin(), packedVoxel.end());

    _sortedMortonCode.resize(voxelCount);
    _sortedIndex.resize(voxelCount);
    for (int n = 0; n < voxelCount; n++) {
      _sortedMortonCode[n] = packedVoxel[n].mortonCode;
      _sortedIndex[n] = packedVoxel[n].index;
    }
  }

  std::vector<Qps> pointQpOffsets(voxelCount);
  for (int i = 0; i < voxelCount; i++)
    pointQpOffsets[i] = {0, 0}; // qpSet.regionQpOffset(pointCloud[index]);
     
  for (auto AttributeGroup : AttributeGroups) {
    // Transform.
    const int attribCount = 3 * AttributeGroup.size();
    std::vector<int> coefficients(attribCount * voxelCount);
    std::vector<int> attributes(attribCount * voxelCount);
    std::vector<QpSet> qpSets;

    int rahtPredThreshold[2];
    bool raht_prediction_enabled_flag;
    for (int i = 0; i < AttributeGroup.size(); i++){
      FixedPoint::kFracBits =_apss[AttributeGroup[i]].raht_fixed_point_precision_minus15 + 15;
      init_LUT();

      auto abh = _allAttrAbh[AttributeGroup[i]];
      int attr_attr_parameter_set_id = abh.attr_attr_parameter_set_id;
      auto start =bufferRahtCoefficients[attr_attr_parameter_set_id].begin();
      auto end = start + 3 * voxelCount;
      std::copy(start, end, coefficients.begin() + 3 * i * voxelCount);
      const auto it_attr_aps = _apss.find(abh.attr_attr_parameter_set_id);
      const auto& attr_aps = it_attr_aps->second;
      rahtPredThreshold[0] = attr_aps.raht_prediction_threshold0;
      rahtPredThreshold[1] = attr_aps.raht_prediction_threshold1;
      raht_prediction_enabled_flag = attr_aps.raht_prediction_enabled_flag;
      const auto& attr_sps = sps.attributeSets[abh.attr_sps_attr_idx];
      QpSet qpSet = deriveQpSet(attr_sps, attr_aps, abh);
      qpSets.push_back(qpSet);
    }
    
    regionAdaptiveHierarchicalInverseTransformforMultiAttr(
      raht_prediction_enabled_flag, rahtPredThreshold, pointQpOffsets.data(),
      qpSets,_sortedMortonCode.data(), attributes.data(),
      attribCount, voxelCount, coefficients.data());

    for (int i = 0; i < AttributeGroup.size(); i++) {
      auto abh = _allAttrAbh[AttributeGroup[i]];
      auto desc = sps.attributeSets[abh.attr_sps_attr_idx];
      int clipMax = (1 << desc.bitdepth) - 1;
      for (int n = 0; n < voxelCount; n++) {
        auto index = _sortedIndex[n];
        Vec3<attr_t> attribute(0);
        for (int d = 0; d < attribCount; ++d) {
          int value = attributes[attribCount * n + 3 * i + d];
          attribute[d] = attr_t(PCCClip(value, 0, clipMax));
        }

        setAttribute(desc, pointCloud, index, attribute);
      }
    }
  }

  bufferRahtCoefficients.clear();
}
//----------------------------------------------------------------------------

void
AttributeDecoder::decodeAttributeLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const QpSet& qpSet,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint64_t> weights;

  if (!aps.scalable_lifting_enabled_flag) {
    computeQuantizationWeights(
      _lods.numPointsInLod, _lods.predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      _lods.numPointsInLod, _lods.predictors, geom_num_points_minus1 + 1,
      minGeomNodeSizeLog2, weights);
  }

  std::vector<Vec3<int64_t>> attributes(pointCount, Vec3<int64_t>(0));

  // decompress
  int zeroRunRem = 0;
  const size_t lodCount = _lods.numPointsInLod.size();
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    int quantLayer = std::min(int(qpSet.layers.size()) - 1, int(lod));

    const bool lcpEnabled = aps.last_component_prediction_enabled_flag && attribCount == 3;
    const int8_t lastCompPredCoeff = lcpEnabled ? abh.attrLcpCoeffs[lod] : 0;

    auto coeffRange = LodCoeffIdxRange(_lods.numPointsInLod, lod);
    for (auto predictorIndex : coeffRange) {
      const uint32_t pointIndex = _lods.indexes[predictorIndex];
      auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);

      if (--zeroRunRem < 0)
        zeroRunRem = decoder.decodeRunLength();

      int32_t values[3] = {0, 0, 0};
      if (!zeroRunRem)
        decoder.decode(values, attribCount);

      const int64_t iQuantWeight = fastIrsqrt(weights[predictorIndex]);

      for (int d = 0; d < attribCount; ++d) {
        const int quantIdx = std::min(d, 1);
        int64_t scaled = quant[quantIdx].scale(values[d]);

        // NB: undo LCP component 2 only
        if (lcpEnabled && d == 2) {
          const int64_t scaled1 = quant[1].scale(values[1]);
          scaled += (int64_t(lastCompPredCoeff) * scaled1) >> 2;
        }

        attributes[predictorIndex][d] = divExp2RoundHalfInf(scaled * iQuantWeight, 36);
      }
    }
  }

  // reconstruct
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    if (lod == 0 /* coarsest */)
      continue;

    auto coeffRange = LodCoeffIdxRange(_lods.numPointsInLod, lod);
    PCCLiftUpdate (_lods.predictors, weights, coeffRange, false, attributes);
    PCCLiftPredict(_lods.predictors, coeffRange, false, attributes);
  }

  int64_t clipMax = (1 << desc.bitdepth) - 1;
  for (size_t f = 0; f < pointCount; ++f) {
    Vec3<attr_t> attribute(0);
    for (int d = 0; d < attribCount; ++d) {
      const auto attribute0 =
        divExp2RoundHalfInf(attributes[f][d], kFixedPointAttributeShift);
      attribute[d] = attr_t(PCCClip(attribute0, 0, clipMax));
    }

    setAttribute(desc, pointCloud, _lods.indexes[f], attribute);
  }
}

//============================================================================

} /* namespace pcc */
