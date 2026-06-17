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

#include "AttributeEncoder.h"
#include <set>
#include "DualLutCoder.h"
#include "attribute_raw.h"
#include "constants.h"
#include "entropy.h"
#include "io_hls.h"
#include "quantization.h"
#include "RAHT.h"
#include "FixedPoint.h"

#include <algorithm>

// todo(df): promote to per-attribute encoder parameter
static const double kAttrPredLambdaR = 0.01;
static const double kAttrPredLambdaC = 0.14;

namespace pcc {
//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

class PCCResidualsEncoder : protected AttributeContexts {
public:
  PCCResidualsEncoder(
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    const AttributeContexts& ctxtMem);

  EntropyEncoder arithmeticEncoder;

  const AttributeContexts& getCtx() const { return *this; }

  void start(const SequenceParameterSet& sps, int numPoints);
  int stop();

  void encodeRunLength(int runLength);
  void encodeSymbol(uint32_t value, int k1, int k2, int k3);
  void encode(const int32_t* values, int attribCount);

  int availPredModes;
  double bitsPtAttrMultiComp(Vec3<int32_t> value, int parity, int attribCount);
  double bitsPtAttrSingleComp(int32_t value, int parity);

  // Encoder side residual cost calculation
  const int scaleRes = 1 << 20;
  const int windowLog2 = 6;
  int probResGt0[3];  //prob of residuals larger than 0: 1 for each component
  int probResGt1[3];  //prob of residuals larger than 1: 1 for each component
  void resStatUpdateAttrMultiComp(Vec3<int32_t> values, int attribCount);
  void resStatUpdateAttrSingleComp(int32_t values);
  void resStatReset();
};

//----------------------------------------------------------------------------

PCCResidualsEncoder::PCCResidualsEncoder(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const AttributeContexts& ctxtMem)
  : AttributeContexts(ctxtMem)
{
  availPredModes =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;

  resStatReset();
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::start(const SequenceParameterSet& sps, int pointCount)
{
  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = pointCount * 3 * 2 + 1024;
  arithmeticEncoder.setBuffer(maxAcBufLen, nullptr);
  arithmeticEncoder.enableBypassStream(sps.cabac_bypass_stream_enabled_flag);
  arithmeticEncoder.start();
}

//----------------------------------------------------------------------------

int
PCCResidualsEncoder::stop()
{
  return arithmeticEncoder.stop();
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatReset()
{
  for (int k = 0; k < 3; k++)
    probResGt0[k] = probResGt1[k] = (scaleRes >> 1);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatUpdateAttrMultiComp(Vec3<int32_t> value, int attribCount)
{
  for (int k = 0; k < attribCount; k++) {
    probResGt0[k] += value[k] ? (scaleRes - probResGt0[k]) >> windowLog2
                              : -((probResGt0[k]) >> windowLog2);
    if (value[k])
      probResGt1[k] += abs(value[k]) > 1
        ? (scaleRes - probResGt1[k]) >> windowLog2
        : -((probResGt1[k]) >> windowLog2);
  }
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::resStatUpdateAttrSingleComp(int32_t value)
{
  probResGt0[0] += value ? (scaleRes - probResGt0[0]) >> windowLog2
                         : -(probResGt0[0] >> windowLog2);
  if (value)
    probResGt1[0] += abs(value) > 1 ? (scaleRes - probResGt1[0]) >> windowLog2
                                    : -(probResGt1[0] >> windowLog2);
}

//----------------------------------------------------------------------------

double
PCCResidualsEncoder::bitsPtAttrMultiComp(Vec3<int32_t> value, int mode, int attribCount)
{
  if (attribCount == 2) {
    if (availPredModes > 1)
      value[1] = availPredModes * abs(value[1]) + mode;
  } else if (availPredModes == 4) {
    value[1] = 2 * abs(value[1]) + (mode >> 1);
    value[2] = 2 * abs(value[2]) + (mode & 1);
  } else if (availPredModes == 3) {
    value[1] = 2 * abs(value[1]) + (mode > 0);
    if (mode > 0)
      value[2] = 2 * abs(value[2]) + (mode - 1);
  } else if (availPredModes == 2) {
    value[1] = 2 * abs(value[1]) + (mode & 1);
  }

  int log2scaleRes = ilog2(uint32_t(scaleRes));
  double bits = 0;
  for (int k = 0; k < attribCount; k++) {
    bits += value[k] ? log2scaleRes - log2(probResGt0[k])
                     : log2scaleRes - log2(scaleRes - probResGt0[k]);  //Gt0
    int mag = abs(value[k]);
    if (mag) {
      bits += mag > 1 ? log2scaleRes - log2(probResGt1[k])
                      : log2scaleRes - log2(scaleRes - probResGt1[k]);  //Gt1
      bits += 1;  //sign bit.
      if (mag > 1)
        bits += 2.0 * log2(mag - 1.0) + 1.0;  //EG0 approximation.
    }
  }

  return bits;
}

//----------------------------------------------------------------------------

double
PCCResidualsEncoder::bitsPtAttrSingleComp(int32_t value, int mode)
{
  if (availPredModes == 4) {
    value = (abs(value) << 2) + (mode);
  } else if (availPredModes == 3) {
    if (mode > 0)
      value = (abs(value) << 1) + (mode - 1);
    value = (abs(value) << 1) + (mode > 0);
  } else if (availPredModes == 2) {
    value = (abs(value) << 1) + (mode & 1);
  }

  int log2scaleRes = ilog2((uint32_t)scaleRes);
  double bits = 0;
  bits += value ? log2scaleRes - log2(probResGt0[0])
                : log2scaleRes - log2(scaleRes - probResGt0[0]);  //Gt0
  int mag = abs(value);
  if (mag) {
    bits += mag > 1 ? log2scaleRes - log2(probResGt1[0])
                    : log2scaleRes - log2(scaleRes - probResGt1[0]);  //Gt1
    bits += 1;  //sign bit.
    if (mag > 1)
      bits += 2.0 * log2(mag - 1.0) + 1.0;  //EG0 approximation.
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeRunLength(int runLength)
{
  auto* ctx = ctxRunLen;
  for (int i = 0; i < std::min(3, runLength); i++, ctx++)
    arithmeticEncoder.encode(1, *ctx);

  if (runLength < 3) {
    arithmeticEncoder.encode(0, *ctx);
    return;
  }
  runLength -= 3;

  auto prefix = runLength >> 1;
  for (int i = 0; i < std::min(4, prefix); i++)
    arithmeticEncoder.encode(1, *ctx);

  if (runLength < 8) {
    arithmeticEncoder.encode(0, *ctx);
    arithmeticEncoder.encode(runLength & 1);
    return;
  }
  runLength -= 8;

  arithmeticEncoder.encodeExpGolomb(runLength, 2, *++ctx);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeSymbol(uint32_t value, int k1, int k2, int k3)
{
  arithmeticEncoder.encode(value > 0, ctxCoeffGtN[0][k1]);
  if (!value)
    return;

  arithmeticEncoder.encode(--value > 0, ctxCoeffGtN[1][k2]);
  if (!value)
    return;

  arithmeticEncoder.encodeExpGolomb(
    --value, 1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(const int32_t* value, int attribCount)
{
  int mag0 = abs(value[0]);

  int bZeroPrev = 0;
  int bLe1Prev  = 0;

  // Packed flags representing zero / <=1 conditions of preceding components
  int zeroBits = 0;
  int le1Bits  = 0;

  // True if all preceding components are zero
  int allZero = 1;

  if (attribCount >= 2) {
    // Encode component 1 (first non-last)
    int mag1 = abs(value[1]);

    int b0 = (mag1 == 0);
    int b1 = (mag1 <= 1);

    bZeroPrev = b0;
    bLe1Prev  = b1;

    int w0 = 1 << (attribCount - 2);
    zeroBits += b0 * w0;
    le1Bits  += b1 * w0;
    allZero  &= b0;

    encodeSymbol(mag1, 0, 0, 1);
    if (mag1)
      arithmeticEncoder.encode(value[1] < 0);

    // Encode remaining non-last components
    for (int k = 2; k < attribCount; ++k) {
      int magK = abs(value[k]);

      int bZeroK = (magK == 0);
      int bLe1K  = (magK <= 1);

      encodeSymbol(magK, 1 + bZeroPrev, 1 + bLe1Prev, 1);
      if (magK)
        arithmeticEncoder.encode(value[k] < 0);

      int w = 1 << (attribCount - 1 - k);
      zeroBits += bZeroK * w;
      le1Bits  += bLe1K  * w;
      allZero  &= bZeroK;

      bZeroPrev = bZeroK;
      bLe1Prev  = bLe1K;
    }
  }

  // Encode component 0 (last)
  // Adjustment matches decoder:
  //  - single component: mag0 - 1
  //  - otherwise: subtract 1 if all previous components are zero
  int mag0minusX =
    (attribCount == 1) ? (mag0 - 1) : (mag0 - allZero);

  assert(mag0minusX >= 0);

  // Contexts derived from accumulated flags
  int ctxZ   = 3 + zeroBits;
  int ctxLe1 = 3 + le1Bits;

  encodeSymbol(mag0minusX, ctxZ, ctxLe1, 0);

  if (mag0)
    arithmeticEncoder.encode(value[0] < 0);
}

//============================================================================
// An encapsulation of the entropy coding methods used in attribute coding

struct PCCResidualsEntropyEstimator {
  size_t freq0[kAttributeResidualAlphabetSize + 1];
  size_t freq1[kAttributeResidualAlphabetSize + 1];
  size_t symbolCount0;
  size_t symbolCount1;
  size_t isZero0Count;
  size_t isZero1Count;
  PCCResidualsEntropyEstimator() { init(); }
  void init();
  double bitsDetail(
    const uint32_t detail,
    const size_t symbolCount,
    const size_t* const freq) const;
  double bits(const uint32_t value0) const;
  void update(const uint32_t value0);
  double bits(
    const uint32_t value0, const uint32_t value1, const uint32_t value2) const;
  void
  update(const uint32_t value0, const uint32_t value1, const uint32_t value2);
};

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::init()
{
  for (size_t i = 0; i <= kAttributeResidualAlphabetSize; ++i) {
    freq0[i] = 1;
    freq1[i] = 1;
  }
  symbolCount0 = kAttributeResidualAlphabetSize + 1;
  symbolCount1 = kAttributeResidualAlphabetSize + 1;
  isZero1Count = isZero0Count = symbolCount0 / 2;
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bitsDetail(
  const uint32_t detail,
  const size_t symbolCount,
  const size_t* const freq) const
{
  const uint32_t detailClipped =
    std::min(detail, uint32_t(kAttributeResidualAlphabetSize));
  const double pDetail =
    PCCClip(double(freq[detailClipped]) / symbolCount, 0.001, 0.999);
  double bits = -log2(pDetail);
  if (detail >= kAttributeResidualAlphabetSize) {
    const double x = double(detail) - double(kAttributeResidualAlphabetSize);
    bits += 2.0 * std::floor(log2(x + 1.0)) + 1.0;
  }
  return bits;
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bits(const uint32_t value0) const
{
  const bool isZero0 = value0 == 0;
  const double pIsZero0 = isZero0
    ? double(isZero0Count) / symbolCount0
    : double(symbolCount0 - isZero0Count) / symbolCount0;
  double bits = -log2(PCCClip(pIsZero0, 0.001, 0.999));
  if (!isZero0) {
    bits += bitsDetail(value0 - 1, symbolCount0, freq0);
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::update(const uint32_t value0)
{
  const bool isZero0 = value0 == 0;
  ++symbolCount0;
  if (!isZero0) {
    ++freq0[std::min(value0 - 1, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero0Count;
  }
}

//----------------------------------------------------------------------------

double
PCCResidualsEntropyEstimator::bits(
  const uint32_t value0, const uint32_t value1, const uint32_t value2) const
{
  const bool isZero0 = value0 == 0;
  const double pIsZero0 = isZero0
    ? double(isZero0Count) / symbolCount0
    : double(symbolCount0 - isZero0Count) / symbolCount0;
  double bits = -log2(PCCClip(pIsZero0, 0.001, 0.999));
  if (!isZero0) {
    bits += bitsDetail(value0 - 1, symbolCount0, freq0);
  }

  const bool isZero1 = value1 == 0 && value2 == 0;
  const double pIsZero1 = isZero1
    ? double(isZero1Count) / symbolCount0
    : double(symbolCount0 - isZero1Count) / symbolCount0;
  bits -= log2(PCCClip(pIsZero1, 0.001, 0.999));
  if (!isZero1) {
    bits += bitsDetail(value1, symbolCount1, freq1);
    bits += bitsDetail(value2, symbolCount1, freq1);
  }
  return bits;
}

//----------------------------------------------------------------------------

void
PCCResidualsEntropyEstimator::update(
  const uint32_t value0, const uint32_t value1, const uint32_t value2)
{
  const bool isZero0 = value0 == 0;
  ++symbolCount0;
  if (!isZero0) {
    ++freq0[std::min(value0 - 1, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero0Count;
  }

  const bool isZero1 = value1 == 0 && value2 == 0;
  symbolCount1 += 2;
  if (!isZero1) {
    ++freq1[std::min(value1, uint32_t(kAttributeResidualAlphabetSize))];
    ++freq1[std::min(value2, uint32_t(kAttributeResidualAlphabetSize))];
  } else {
    ++isZero1Count;
  }
}

//============================================================================
// AttributeEncoderIntf

AttributeEncoderIntf::~AttributeEncoderIntf() = default;

//============================================================================
// AttributeEncoder factory

std::unique_ptr<AttributeEncoderIntf>
makeAttributeEncoder()
{
  return std::unique_ptr<AttributeEncoder>(new AttributeEncoder());
}

//============================================================================
// AttributeEncoder Members

void
AttributeEncoder::encode(
  const SequenceParameterSet& sps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  AttributeBrickHeader& abh,
  AttributeContexts& ctxtMem,
  PCCPointSet3& pointCloud,
  PayloadBuffer* payload,
  std::unordered_map<int, std::vector<int>> *bufferRahtCoefficients)
{
  if (attr_aps.attr_encoding == AttributeEncoding::kRaw) {
    AttrRawEncoder::encode(sps, desc, attr_aps, abh, pointCloud, payload);
    return;
  }

  // Encoders are able to modify the slice header:
  _abh = &abh;

  QpSet qpSet = deriveQpSet(desc, attr_aps, abh);

  // generate LoDs if necessary
  if (attr_aps.lodParametersPresent() && _lods.empty())
    _lods.generate(
      attr_aps, abh, pointCloud.getPointCount() - 1, 0, pointCloud);

  PCCResidualsEncoder encoder(attr_aps, abh, ctxtMem);
  encoder.start(sps, int(pointCloud.getPointCount()));

  switch (attr_aps.attr_encoding) {
  case AttributeEncoding::kRAHTransform:
    encodeAttributeTransformRaht(
      desc, attr_aps, qpSet, pointCloud, encoder, bufferRahtCoefficients);
    break;

  case AttributeEncoding::kPredictingTransform:
    encodeAttributePred(desc, attr_aps, qpSet, pointCloud, encoder);
    break;

  case AttributeEncoding::kLiftingTransform:
    encodeAttributeLift(desc, attr_aps, qpSet, pointCloud, encoder);
    break;

  case AttributeEncoding::kRaw:
    // Already handled
    break;
  }

  uint32_t acDataLen = encoder.stop();

  // write abh
  write(sps, attr_aps, abh, payload);
  _abh = nullptr;

  std::copy_n(
    encoder.arithmeticEncoder.buffer(), acDataLen,
    std::back_inserter(*payload));

  // save the context state for re-use by a future slice if required
  ctxtMem = encoder.getCtx();
}

//----------------------------------------------------------------------------

bool
AttributeEncoder::isReusable(
  const AttributeParameterSet& aps, const AttributeBrickHeader& abh) const
{
  return _lods.isReusable(aps, abh);
}

//----------------------------------------------------------------------------

int64_t
AttributeEncoder::computeAttrSingleCompResidual(
  const uint64_t attribute,
  const uint64_t predictedAttribute,
  const Quantizer& quant)
{
  const int64_t quantAttValue = attribute;
  const int64_t quantPredAttValue = predictedAttribute;
  const int64_t delta = quant.quantize(
    (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
  return delta;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::decidePredModeAttrSingleComp(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexesLOD,
  const uint32_t predictorIndex,
  PCCPredictor& predictor,
  PCCResidualsEncoder& encoder,
  PCCResidualsEntropyEstimator& context,
  const Quantizer& quant)
{
  uint64_t attrValue = uint64_t(getAttribute(desc, pointCloud, indexesLOD[predictorIndex])[0]);

  // base case: start with the first neighbour
  // NB: skip evaluation of mode 0 (weighted average of n neighbours)
  int startPredIndex = aps.direct_avg_predictor_disabled_flag;
  predictor.predMode = startPredIndex;
  uint64_t attrPred = predictor.predictAttrSingleComp(pointCloud, indexesLOD, desc);
  int64_t attrResidualQuant =
    computeAttrSingleCompResidual(attrValue, attrPred, quant);

  // NB: idxBits is not included in the score
  int mode = predictor.predMode - aps.direct_avg_predictor_disabled_flag;
  int64_t best_score = encoder.bitsPtAttrSingleComp(attrResidualQuant, mode);

  for (int i = startPredIndex; i < predictor.neighborCount; i++) {
    if (i == aps.max_num_direct_predictors)
      break;
    const uint32_t neighPointIndex =
      indexesLOD[predictor.neighbors[i].predictorIndex];
    attrPred = uint64_t(getAttribute(desc, pointCloud, neighPointIndex)[0]);
    attrResidualQuant = computeAttrSingleCompResidual(attrValue, attrPred, quant);

    mode = i + !aps.direct_avg_predictor_disabled_flag;
    const int64_t score = encoder.bitsPtAttrSingleComp(attrResidualQuant, mode);
    if (score < best_score) {
      best_score = score;
      predictor.predMode = i + 1;
      // NB: setting predictor.neighborCount = 1 will cause issues
      // with reconstruction.
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodePredModeAttrSingleComp(
  const AttributeParameterSet& aps, int predMode, int32_t& coeff)
{
  int coeffSign = coeff < 0 ? -1 : 1;
  int coeffAbs = abs(coeff);

  int mode = predMode - aps.direct_avg_predictor_disabled_flag;
  int maxcand =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;
  switch (maxcand) {
  case 4: coeff = coeffSign * ((coeffAbs << 2) + mode); break;

  case 3:
    if (mode > 0)
      coeffAbs = ((coeffAbs << 1) + (mode - 1));
    coeffAbs = ((coeffAbs << 1) + (mode > 0));
    coeff = coeffSign * coeffAbs;
    break;

  case 2: coeff = coeffSign * ((coeffAbs << 1) + mode); break;

  default: assert(mode == 0);
  }
}

//----------------------------------------------------------------------------

Vec3<int64_t>
AttributeEncoder::computeAttrMultiCompResiduals(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const Vec3<attr_t> attribute,
  const Vec3<attr_t> predictedAttribute,
  const Vec3<int8_t> icpCoeff,
  const Quantizers& quant)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  Vec3<int64_t> residuals;
  const int64_t quantAttValue = attribute[0];
  const int64_t quantPredAttValue = predictedAttribute[0];
  const int64_t delta = quant[0].quantize(
    (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
  residuals[0] = delta;
  const int64_t residual0 =
    divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
  for (size_t k = 1; k < attribCount; ++k) {
    const int64_t quantAttValue = attribute[k];
    const int64_t quantPredAttValue = predictedAttribute[k];

    if (aps.inter_component_prediction_enabled_flag) {
      auto err = quantAttValue - quantPredAttValue
        - ((icpCoeff[k] * residual0 + 2) >> 2);

      auto delta = quant[1].quantize(err << kFixedPointAttributeShift);
      residuals[k] = delta;
    } else {
      const int64_t delta = quant[1].quantize(
        (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift);
      residuals[k] = delta;
    }
  }
  return residuals;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::decidePredModeAttrMultiComp(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const PCCPointSet3& pointCloud,
  const std::vector<uint32_t>& indexesLOD,
  const uint32_t predictorIndex,
  PCCPredictor& predictor,
  PCCResidualsEncoder& encoder,
  PCCResidualsEntropyEstimator& context,
  const Vec3<int8_t>& icpCoeff,
  const Quantizers& quant)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;

  const uint32_t pointIndex = indexesLOD[predictorIndex];
  Vec3<attr_t> attrValue = getAttribute(desc, pointCloud, pointIndex);

  // base case: weighted average of n neighbours
  int startpredIndex = aps.direct_avg_predictor_disabled_flag;
  predictor.predMode = startpredIndex;
  Vec3<attr_t> attrPred = predictor.predictAttrMultiComp(pointCloud, indexesLOD, desc);

  Vec3<int64_t> attrResidualQuant =
    computeAttrMultiCompResiduals(desc, aps, attrValue, attrPred, icpCoeff, quant);
  auto attrDistortion =
    computeAttrMultiCompDistortions(desc, attrValue, attrPred, quant);

  double rate = encoder.bitsPtAttrMultiComp(attrResidualQuant, 0, attribCount);
  double best_score = attrDistortion
    + rate * kAttrPredLambdaC
      * (quant[0].stepSize() >> kFixedPointAttributeShift);

  for (int i = startpredIndex; i < predictor.neighborCount; i++) {
    if (i == aps.max_num_direct_predictors)
      break;

    const uint32_t neighPointIndex = indexesLOD[predictor.neighbors[i].predictorIndex];
    attrPred = getAttribute(desc, pointCloud, neighPointIndex);

    attrResidualQuant =
      computeAttrMultiCompResiduals(desc, aps, attrValue, attrPred, icpCoeff, quant);
    attrDistortion = computeAttrMultiCompDistortions(desc, attrValue, attrPred, quant);

    int sigIdx = i + !aps.direct_avg_predictor_disabled_flag;
    double rate = encoder.bitsPtAttrMultiComp(attrResidualQuant, sigIdx, attribCount);
    double score = attrDistortion
      + rate * kAttrPredLambdaC
        * (quant[0].stepSize() >> kFixedPointAttributeShift);

    if (score < best_score) {
      best_score = score;
      predictor.predMode = i + 1;
      // NB: setting predictor.neighborCount = 1 will cause issues
      // with reconstruction.
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodePredModeAttrMultiComp(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps, int predMode, Vec3<int32_t>& coeff)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;

  int signk1 = coeff[1] < 0 ? -1 : 1;
  int signk2 = coeff[2] < 0 ? -1 : 1;
  int coeffAbsk1 = abs(coeff[1]);
  int coeffAbsk2 = abs(coeff[2]);

  int mode = predMode - aps.direct_avg_predictor_disabled_flag;
  int maxcand =
    aps.max_num_direct_predictors + !aps.direct_avg_predictor_disabled_flag;

  if (attribCount == 2) {
    int signk1 = coeff[1] < 0 ? -1 : 1;
    int coeffAbsk1 = abs(coeff[1]);
    coeff[1] = signk1 * (coeffAbsk1 * maxcand + mode);
    return;
  }

  switch (maxcand) {
    int parityk1, parityk2;
  case 4:
    parityk1 = mode >> 1;
    parityk2 = mode & 1;
    coeff[1] = signk1 * ((coeffAbsk1 << 1) + parityk1);
    coeff[2] = signk2 * ((coeffAbsk2 << 1) + parityk2);
    break;

  case 3:
    parityk1 = mode ? 1 : 0;
    coeff[1] = signk1 * ((coeffAbsk1 << 1) + parityk1);
    if (parityk1) {
      parityk2 = mode - parityk1;
      coeff[2] = signk2 * ((coeffAbsk2 << 1) + parityk2);
    }
    break;

  case 2:
    parityk1 = mode;
    coeff[1] = signk1 * ((coeffAbsk1 << 1) + parityk1);
    break;

  default: assert(mode == 0);
  }
}

//----------------------------------------------------------------------------

std::vector<Vec3<int8_t>>
AttributeEncoder::computeInterComponentPredictionCoeffs(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps, const PCCPointSet3& pointCloud)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  int maxNumDetailLevels = aps.maxNumDetailLevels();
  assert(_lods.numPointsInLod.size() <= maxNumDetailLevels);

  // Two secondary colour components (positive sign set)
  // NB: k=0 is never used
  std::vector<Vec3<int8_t>> signs(maxNumDetailLevels, {0, 1, 1});

  // Estimate residual using original neighbour as predictor
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<Vec3<int32_t>> residual(pointCount);

  for (size_t predIdx = 0; predIdx < pointCount; ++predIdx) {
    const auto pointIdx = _lods.indexes[predIdx];
    auto& predictor = _lods.predictors[predIdx];

    // taking first neighbor for simplicity
    predictor.predMode = 1;
    auto predAttr = predictor.predictAttrMultiComp(pointCloud, _lods.indexes, desc);

    Vec3<attr_t> srcAttr = getAttribute(desc, pointCloud, pointIdx);

    residual[predIdx] = Vec3<int>(srcAttr) - Vec3<int>(predAttr);

    // reset is needed as RD would be done later.
    predictor.predMode = 0;
  }

  const int nWeights = 8;
  const int nShift = 2;  // from log2(nWeights >> 1)
  std::vector<Vec3<int64_t>> sumPredCoeff(nWeights, 0);
  Vec3<int64_t> sumOrigCoeff = 0;

  const size_t lodCount = _lods.numPointsInLod.size();
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    for (auto predIdx : LodCoeffIdxRange(_lods.numPointsInLod, lod)) {
      Vec3<int32_t> resid = residual[predIdx];

      for (int w = 0; w < nWeights; w++)
        for (int k = 1; k < attribCount; k++) {
          auto pred = signs[lod][k] * (((w + 1) * resid[0] + 2) >> nShift);
          sumPredCoeff[w][k] += abs(resid[k] - pred);
        }

      for (int k = 1; k < attribCount; k++)
        sumOrigCoeff[k] += abs(resid[k]);
    }

    // find the best weight
    for (int k = 1; k < attribCount; k++) {
      auto best = std::min_element(
        sumPredCoeff.begin(), sumPredCoeff.end(),
        [=](Vec3<int64_t>& a, Vec3<int64_t>& b) { return a[k] < b[k]; });

      int coeff = 1 + std::distance(sumPredCoeff.begin(), best);
      signs[lod][k] *= coeff;

      assert(signs[lod][k] < nWeights + 1 && signs[lod][k] > -(nWeights + 1));

      if ((*best)[k] > sumOrigCoeff[k])
        signs[lod][k] = 0;
    }

    for (int w = 0; w < nWeights; w++)
      sumPredCoeff[w] = 0;
    sumOrigCoeff = 0;
  }

  // NB: there may be more coefficients than actual detail levels
  // Set any unused detail level coefficients to 0
  for (auto lod = lodCount; lod < maxNumDetailLevels; lod++)
    signs[lod] = 0;

  return signs;
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeAttributePred(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  const bool isSingleComp = (attribCount == 1);

  int64_t clipMax = (int64_t(1) << desc.bitdepth) - 1;
  PCCResidualsEntropyEstimator context;
  int zeroRunAcc = 0;

  bool icpPresent = false;
  if (aps.inter_component_prediction_enabled_flag && !isSingleComp) {
    icpPresent = _abh->icpPresent(desc, aps);
    if (icpPresent)
      _abh->icpCoeffs = computeInterComponentPredictionCoeffs(desc, aps, pointCloud);
  }

  const size_t lodCount = _lods.numPointsInLod.size();
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    int quantLayer = std::min<int>(qpSet.layers.size() - 1, int(lod));
    Vec3<int8_t> icpCoeff = (icpPresent && !isSingleComp) ? _abh->icpCoeffs[lod] : Vec3<int8_t>(0);

    auto coeffRange = LodCoeffIdxRange(_lods.numPointsInLod, lod);
    for (auto predictorIndex : coeffRange) {
      const auto pointIndex = _lods.indexes[predictorIndex];
      auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);
      auto& predictor = _lods.predictors[predictorIndex];
      predictor.predMode = 0;

      if (isSingleComp) {
        bool predModeEligible =
          predModeEligibleAttrSingleComp(desc, aps, pointCloud, _lods.indexes, predictor);
        if (predModeEligible)
          decidePredModeAttrSingleComp(
            desc, aps, pointCloud, _lods.indexes, predictorIndex, predictor,
            encoder, context, quant[0]);

        Vec3<attr_t> attribute = getAttribute(desc, pointCloud, pointIndex);
        const attr_t predictedAttribute =
          predictor.predictAttrSingleComp(pointCloud, _lods.indexes, desc);
        const int64_t quantAttValue = attribute[0];
        const int64_t quantPredAttValue = predictedAttribute;
        const int32_t delta = int32_t(quant[0].quantize(
            (quantAttValue - quantPredAttValue) << kFixedPointAttributeShift));
        int32_t attValue0 = delta;

        if (predModeEligible)
          encodePredModeAttrSingleComp(aps, predictor.predMode, attValue0);

        const int64_t reconstructedDelta =
          divExp2RoundHalfUp(quant[0].scale(delta), kFixedPointAttributeShift);
        const int64_t reconstructedQuantAttValue =
          quantPredAttValue + reconstructedDelta;
        Vec3<attr_t> reconstructedAttribute(0);
          reconstructedAttribute[0] = attr_t(PCCClip(reconstructedQuantAttValue, int64_t(0), clipMax));

        if (!attValue0) {
          ++zeroRunAcc;
        } else {
          encoder.encodeRunLength(zeroRunAcc);
          encoder.encode(&attValue0, 1);
          zeroRunAcc = 0;
        }
        setAttribute(desc, pointCloud, pointIndex, reconstructedAttribute);
        encoder.resStatUpdateAttrSingleComp(attValue0);
      } else {
        bool predModeEligible =
          predModeEligibleAttrMultiComp(desc, aps, pointCloud, _lods.indexes, predictor);
        if (predModeEligible)
          decidePredModeAttrMultiComp(
            desc, aps, pointCloud, _lods.indexes, predictorIndex, predictor,
            encoder, context, icpCoeff, quant);
        Vec3<attr_t> attribute = getAttribute(desc, pointCloud, pointIndex);

        const Vec3<attr_t> predictedAttribute =
          predictor.predictAttrMultiComp(pointCloud, _lods.indexes, desc);

        Vec3<int32_t> values = 0;
        Vec3<attr_t> reconstructedAttribute = 0;
        int64_t residual0 = 0;
        for (int k = 0; k < attribCount; ++k) {
          const auto& q = quant[std::min(k, 1)];
          int64_t residual = int64_t(attribute[k]) - int64_t(predictedAttribute[k]);

          int64_t residualQ = q.quantize(residual << kFixedPointAttributeShift);
          int64_t residualR =
            divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);

          if (aps.inter_component_prediction_enabled_flag && k > 0) {
            residual = residual - ((icpCoeff[k] * residual0 + 2) >> 2);
            residualQ = q.quantize(residual << kFixedPointAttributeShift);
            residualR = ((icpCoeff[k] * residual0 + 2) >> 2)
              + divExp2RoundHalfUp(
                q.scale(residualQ), kFixedPointAttributeShift);
          }

          if (k == 0)
            residual0 = residualR;

          values[k] = int32_t(residualQ);

          int64_t recon = int64_t(predictedAttribute[k]) + residualR;
          reconstructedAttribute[k] = attr_t(PCCClip(recon, int64_t(0), clipMax));
        }

        if (predModeEligible)
          encodePredModeAttrMultiComp(desc, aps, predictor.predMode, values);

        setAttribute(desc, pointCloud, pointIndex, reconstructedAttribute);

        encoder.resStatUpdateAttrMultiComp(values, attribCount);

        bool allZero = true;
        for (int k = 0; k < attribCount; ++k)
          allZero &= (values[k] == 0);

        if (allZero) {
          ++zeroRunAcc;
        } else {
          encoder.encodeRunLength(zeroRunAcc);
          int32_t residualCoeffs[3] = { values[0], values[1], values[2] };
          encoder.encode(residualCoeffs, attribCount);
          zeroRunAcc = 0;
        }
      }
    }
  }
  if (zeroRunAcc)
    encoder.encodeRunLength(zeroRunAcc);
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeAttributeTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  std::unordered_map<int, std::vector<int>>* bufferRahtCoefficients)
{
  const int voxelCount = int(pointCloud.getPointCount());
  if (_sortedMortonCode.size() == 0)
  {
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

  // Allocate arrays
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  auto aps_attr_parameter_set_id = aps.aps_attr_parameter_set_id;
  auto it = (*bufferRahtCoefficients).find(aps_attr_parameter_set_id);
  if (it != bufferRahtCoefficients->end()) {
    // Entropy encode.
    std::vector<int> values(attribCount);
    int zeroRun = 0;
    for (int n = 0; n < voxelCount; ++n) {
      bool allZero = true;
      for (int d = 0; d < attribCount; ++d) {
        values[d] = it->second[voxelCount * d + n];
        if (values[d] != 0) allZero = false;
      }
      if (allZero)
        ++zeroRun;
      else {
        encoder.encodeRunLength(zeroRun);
        encoder.encode(values.data(), attribCount);
        zeroRun = 0;
      }
    }
    if (zeroRun)
      encoder.encodeRunLength(zeroRun);
      
    return;
  } 

  std::vector<int> attributes(attribCount * voxelCount);
  std::vector<int> coefficients(attribCount * voxelCount);
  std::vector<Qps> pointQpOffsets(voxelCount);

  // Populate input arrays.
  for (int n = 0; n < voxelCount; n++) {
    auto& index = _sortedIndex[n];
    Vec3<attr_t> attribute = getAttribute(desc, pointCloud, index);
    for (int d = 0; d < attribCount; ++d)
      attributes[attribCount * n + d] = int(attribute[d]);
    pointQpOffsets[n] = qpSet.regionQpOffset(pointCloud[index]);
  }

  const int rahtPredThreshold[2] = {aps.raht_prediction_threshold0,
                                    aps.raht_prediction_threshold1};

  FixedPoint::kFracBits = aps.raht_fixed_point_precision_minus15 + 15;
  init_LUT();
  // Transform.
  regionAdaptiveHierarchicalTransform(
    aps.raht_prediction_enabled_flag, rahtPredThreshold, qpSet,
    pointQpOffsets.data(), _sortedMortonCode.data(), attributes.data(), attribCount,
    voxelCount, coefficients.data());

  // Entropy encode
  std::vector<int> values(attribCount);
  int zeroRun = 0;
  for (int n = 0; n < voxelCount; ++n) {
    bool allZero = true;
    for (int d = 0; d < attribCount; ++d) {
      values[d] = coefficients[voxelCount * d + n];
      if (values[d] != 0) allZero = false;
    }
    if (allZero)
      ++zeroRun;
    else {
      encoder.encodeRunLength(zeroRun);
      encoder.encode(values.data(), attribCount);
      zeroRun = 0;
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

  int64_t clipMax = (1ll << desc.bitdepth) - 1;
  for (int n = 0; n < voxelCount; ++n) {
    auto index = _sortedIndex[n];
    Vec3<attr_t> attribute(0);
    for (int d = 0; d < attribCount; ++d) {
      int value = attributes[attribCount * n + d];
      attribute[d] = attr_t(PCCClip(value, 0, clipMax));
    }

    setAttribute(desc, pointCloud, index, attribute);
  }
}

//----------------------------------------------------------------------------

std::vector<std::vector<int>>
AttributeGroup(
  std::vector<std::pair<std::string, int>> attributeIdxMap,
  std::vector<AttributeParameterSet*> attr_aps,
  int Group_maxNum)
{
  std::vector<std::vector<int>> Groups;
  std::vector<int> Group;
  std::set<std::array<int, 2>> uniquePredictionThreshold;
  for (const auto& it : attributeIdxMap) {
    int attrIdx = it.second;
    if (attr_aps[attrIdx]->raht_prediction_enabled_flag) {
      auto predictionThreshold0 =
        attr_aps[attrIdx]->raht_prediction_threshold0;
      auto predictionThreshold1 =
        attr_aps[attrIdx]->raht_prediction_threshold1;
      uniquePredictionThreshold.insert(
        {predictionThreshold0, predictionThreshold1});
    } else {
      Group.push_back(attrIdx);
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
    for (const auto& it : attributeIdxMap) {
      int attrIdx = it.second;
      if (attr_aps[attrIdx]->raht_prediction_enabled_flag) {
        auto predictionThreshold0 =
          attr_aps[attrIdx]->raht_prediction_threshold0;
        auto predictionThreshold1 =
          attr_aps[attrIdx]->raht_prediction_threshold1;
        if (
          predictionThreshold0 == Threshold[0]
          && predictionThreshold1 == Threshold[1]) {
          Group.push_back(attrIdx);
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
AttributeEncoder::encodeMultiAttributeRAHT(
  const SequenceParameterSet& _sps,
  std::vector<AttributeParameterSet*> attr_aps,
  PCCPointSet3& pointCloud,
  std::vector<std::pair<std::string, int>> attributeIdxMap,
  int maxAttributeNumberPerGroup,
  std::unordered_map<int, std::vector<int>>* bufferRahtCoefficients)
{

  int rahtAttrNum = 0;
  for (auto it = attributeIdxMap.begin(); it != attributeIdxMap.end();) {
    int attrIdx = it->second;
    const auto& sps = _sps.attributeSets[attrIdx];
    const auto& aps = attr_aps[attrIdx];

    if (
      aps->attr_encoding == AttributeEncoding::kRAHTransform
      && sps.attr_num_dimensions_minus1 == 2) {
      rahtAttrNum++;
      ++it;
    } else {
      it = attributeIdxMap.erase(it);
    }
  }
  if (rahtAttrNum <= 1)
    return;
  
  std::vector<std::vector<int>> AttributeGroups = AttributeGroup(
    attributeIdxMap, attr_aps, maxAttributeNumberPerGroup);  //Attribute Set for Multi-Attribute Encoding
 
  const int voxelCount = int(pointCloud.getPointCount());
  if (_sortedMortonCode.size() == 0)
  {
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

  // Allocate arrays.
  for (auto AttributeGroup : AttributeGroups) {
    const int attribCount = 3 * AttributeGroup.size();

    std::vector<int> attributes(attribCount * voxelCount);
    std::vector<int> coefficients(attribCount * voxelCount);
    std::vector<Qps> pointQpOffsets(voxelCount);

    // Populate input arrays.
    for (int AttributeGroupIdx = 0; AttributeGroupIdx < AttributeGroup.size();
         AttributeGroupIdx++) {
      auto attrIdx = AttributeGroup[AttributeGroupIdx];
      auto desc = _sps.attributeSets[attrIdx];

      for (int n = 0; n < voxelCount; n++) {
        auto& index = _sortedIndex[n];
        Vec3<attr_t> attribute = getAttribute(desc, pointCloud, index);
        for (int d = 0; d < attribCount; ++d)
          attributes[attribCount * n + 3 * AttributeGroupIdx + d] = attribute[d];
        pointQpOffsets[n] = {0, 0}; // qpSet.regionQpOffset(pointCloud[index]);
      }
    }
    
    int rahtPredThreshold[2];
    bool raht_prediction_enabled_flag;
    // Transform.
    FixedPoint::kFracBits = attr_aps[AttributeGroup[0]]->raht_fixed_point_precision_minus15 + 15;

    init_LUT();

    std::vector<QpSet> qpSets;

    for (const auto& attrIdx : AttributeGroup) {
      auto attr_sps = _sps.attributeSets[attrIdx];
      auto aps = attr_aps[attrIdx];
      rahtPredThreshold[0] = aps->raht_prediction_threshold0;
      rahtPredThreshold[1] = aps->raht_prediction_threshold1;
      raht_prediction_enabled_flag = aps->raht_prediction_enabled_flag;
      AttributeBrickHeader abh;
      abh.attr_sps_attr_idx = attrIdx;
      abh.attr_qp_delta_luma = 0;
      abh.attr_qp_delta_chroma = 0;

      QpSet qpSet = deriveQpSet(attr_sps, *aps, abh);
      qpSets.push_back(qpSet);
    }

    regionAdaptiveHierarchicalTransformforMultiAttr(
      raht_prediction_enabled_flag, rahtPredThreshold, pointQpOffsets.data(),qpSets,
      _sortedMortonCode.data(), attributes.data(), attribCount, voxelCount,
      coefficients.data());

    for (int AttributeGroupIdx = 0; AttributeGroupIdx < AttributeGroup.size();
         AttributeGroupIdx++) {
      auto attrIdx = AttributeGroup[AttributeGroupIdx];
      auto desc = _sps.attributeSets[attrIdx];

      std::vector<int> bufferRahtCoefficient;
      bufferRahtCoefficient.reserve(3 * voxelCount);

      auto start = coefficients.begin() + 3 * AttributeGroupIdx * voxelCount;
      auto end = start + 3 * voxelCount;
      bufferRahtCoefficient.assign(start, end); 
      auto aps_attr_parameter_set_id =attr_aps[attrIdx]->aps_attr_parameter_set_id;
      (*bufferRahtCoefficients)[aps_attr_parameter_set_id] =bufferRahtCoefficient;

      int clipMax = (1 << desc.bitdepth) - 1;
      for (int n = 0; n < voxelCount; n++) {
        auto index = _sortedIndex[n];
        Vec3<attr_t> attribute(0);
        for (int d = 0; d < attribCount; ++d) {
          int value = attributes[attribCount * n + 3 * AttributeGroupIdx + d];
          attribute[d] = attr_t(PCCClip(value, 0, clipMax));
        }

        setAttribute(desc, pointCloud, index, attribute);
      }
    }
  }
}

//----------------------------------------------------------------------------
void
AttributeEncoder::encodeAttributeLift(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder)
{
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  const size_t pointCount = pointCloud.getPointCount();
  std::vector<uint64_t> weights;

  if (!aps.scalable_lifting_enabled_flag) {
    computeQuantizationWeights(
      _lods.numPointsInLod, _lods.predictors, weights);
  } else {
    computeQuantizationWeightsScalable(
      _lods.numPointsInLod, _lods.predictors, pointCount, 0, weights);
  }

  std::vector<int64_t> iQuantWeights(pointCount);
  std::vector<int64_t> quantWeights(pointCount);
  for (size_t i = 0; i < pointCount; ++i) {
    iQuantWeights[i] = fastIrsqrt(weights[i]);
    quantWeights[i] = divExp2RoundHalfUp(weights[i] * iQuantWeights[i], 28);
  }

  std::vector<Vec3<int64_t>> attributes(pointCount, Vec3<int64_t>(0));

  for (size_t index = 0; index < pointCount; ++index) {
    Vec3<attr_t> attribute = getAttribute(desc, pointCloud, _lods.indexes[index]);

    for (int d = 0; d < attribCount; ++d) {
      attributes[index][d] =
        int64_t(attribute[d]) << kFixedPointAttributeShift;
    }
  }

  const size_t lodCount = _lods.numPointsInLod.size();
  for (auto lod : LodOrder::FinestToCoarsest(lodCount)) {
    if (lod == 0 /* coarsest */)
      break;

    auto coeffRange = LodCoeffIdxRange(_lods.numPointsInLod, lod);
    PCCLiftPredict(_lods.predictors, coeffRange, true, attributes);
    PCCLiftUpdate (_lods.predictors, weights, coeffRange, true, attributes);
  }

  const bool lcpEnabled = aps.last_component_prediction_enabled_flag && attribCount == 3;

  // Per level-of-detail coefficients for last component prediction
  if (lcpEnabled)
    _abh->attrLcpCoeffs = computeLastComponentPredictionCoeff(aps, attributes);

  int zeroRun = 0;
  std::vector<int32_t> values(attribCount);
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    int quantLayer = std::min(int(qpSet.layers.size()) - 1, int(lod));

    int8_t lastCompPredCoeff =
      lcpEnabled ? _abh->attrLcpCoeffs[lod] : 0;

    auto coeffRange = LodCoeffIdxRange(_lods.numPointsInLod, lod);
    for (auto predictorIndex : coeffRange) {
      const auto pointIndex = _lods.indexes[predictorIndex];
      auto quant = qpSet.quantizers(pointCloud[pointIndex], quantLayer);

      const int64_t iQuantWeight = iQuantWeights[predictorIndex];
      const int64_t quantWeight = quantWeights[predictorIndex];

      int64_t reconstructedDelta = 0;

      for (int d = 0; d < attribCount; ++d) {
        int quantIdx = std::min(d, 1);

        int64_t attribute = attributes[predictorIndex][d];

        // LCP subtraction before quant
        if (lcpEnabled && d == 2)
          attribute -= (int64_t(lastCompPredCoeff) * attributes[predictorIndex][1]) >> 2;

        values[d] = quant[quantIdx].quantize((attribute * quantWeight) >> 16);
        int64_t scaled = quant[quantIdx].scale(values[d]);

        // LCP add-back after inverse quant
        if (lcpEnabled && d == 2)
          scaled += (int64_t(lastCompPredCoeff) * reconstructedDelta) >> 2;

        attributes[predictorIndex][d] = divExp2RoundHalfInf(scaled * iQuantWeight, 36);

        if (d == 1)
          reconstructedDelta = scaled;
      }

      bool allZero = true;
      for (int d = 0; d < attribCount; ++d)
        allZero &= (values[d] == 0);

      if (allZero)
        ++zeroRun;
      else {
        encoder.encodeRunLength(zeroRun);
        encoder.encode(values.data(), attribCount);
        zeroRun = 0;
      }
    }
  }
  if (zeroRun)
    encoder.encodeRunLength(zeroRun);

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

//----------------------------------------------------------------------------

std::vector<int8_t>
AttributeEncoder::computeLastComponentPredictionCoeff(
  const AttributeParameterSet& aps, const std::vector<Vec3<int64_t>>& coeffs)
{
  int maxNumDetailLevels = aps.maxNumDetailLevels();
  assert(_lods.numPointsInLod.size() <= maxNumDetailLevels);
  std::vector<int8_t> signs(maxNumDetailLevels, 0);

  const size_t lodCount = _lods.numPointsInLod.size();
  for (auto lod : LodOrder::CoarsestToFinest(lodCount)) {
    int64_t sumk1k2 = 0;
    int64_t sumk1k1 = 0;

    for (auto coeffIdx : LodCoeffIdxRange(_lods.numPointsInLod, lod)) {
      auto& attr = coeffs[coeffIdx];
      int mult = attr[1] * attr[2];
      int mult2 = attr[1] * attr[1];
      sumk1k2 += mult;
      sumk1k1 += mult2;
    }

    int scale = 0;
    if (sumk1k2 && sumk1k1) {
      // sign(sumk1k2) * sign(sumk1k1)
      int sign = (sumk1k2 < 0) ^ (sumk1k1 < 0) ? -1 : 1;
      scale = ((sumk1k2 << 2) + sign * (sumk1k1 >> 1)) / sumk1k1;
    }

    // NB: coding range is limited to +-8
    signs[lod] = PCCClip(scale, -8, 8);
  }

  // NB: there may be more coefficients than actual detail levels
  // Propagate the last value to all unused levels to minimise useless cost
  for (auto lod = lodCount; lod < maxNumDetailLevels; lod++)
    signs[lod] = signs[lod - 1];

  return signs;
}

//============================================================================

int
AttributeEncoder::computeAttrMultiCompDistortions(
  const AttributeDescription& desc,
  const Vec3<attr_t> attribute,
  const Vec3<attr_t> predictedAttribute,
  const Quantizers& quant)
{
  int64_t clipMax = (1 << desc.bitdepth) - 1;

  Vec3<attr_t> reconstructedAttribute;
  const int attribCount = desc.attr_num_dimensions_minus1 + 1;
  for (int k = 0; k < attribCount; ++k) {
    const auto& q = quant[std::min(k, 1)];
    int64_t residual = attribute[k] - predictedAttribute[k];

    int64_t residualQ = q.quantize(residual << kFixedPointAttributeShift);
    int64_t residualR =
      divExp2RoundHalfUp(q.scale(residualQ), kFixedPointAttributeShift);

    int64_t recon = predictedAttribute[k] + residualR;
    reconstructedAttribute[k] = attr_t(PCCClip(recon, int64_t(0), clipMax));
  }

  int distortion = 0;
  for (int k = 0; k < attribCount; ++k)
    distortion += std::abs(attribute[k] - reconstructedAttribute[k]);

  return distortion;
}

//============================================================================
// estimation of dist2

int
estimateDist2(
  const PCCPointSet3& cloud,
  int32_t samplingPeriod,
  int32_t searchRange,
  float percentileEstimate)
{
  int32_t pointCount = cloud.getPointCount();
  if (pointCount < 2)
    return 0;

  std::vector<int64_t> dists;
  dists.reserve(pointCount / samplingPeriod + 1);

  for (int32_t index = 0; index < pointCount; index += samplingPeriod) {
    auto k0 = std::max(0, index - searchRange);
    auto k1 = std::min(pointCount - 1, index + searchRange);
    auto d2 = std::numeric_limits<int64_t>::max();
    for (auto k = k0; k <= k1; ++k) {
      if (k == index)
        continue;

      d2 = std::min(d2, (cloud[index] - cloud[k]).getNorm2<int64_t>());
    }
    dists.push_back(d2);
  }

  int p = int(std::floor(dists.size() * percentileEstimate));

  std::nth_element(dists.begin(), dists.begin() + p, dists.end());
  int64_t dist2 = dists[p];
  int shiftBits = 0;
  while ((int64_t(3) << (shiftBits << 1)) < dist2 && shiftBits < 20)
    ++shiftBits;

  return shiftBits;
}

//============================================================================

} /* namespace pcc */
