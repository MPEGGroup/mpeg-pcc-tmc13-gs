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

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace pcc {

//============================================================================

static double
scalePowerUniformStep(const AttributeParameters& params)
{
  const int64_t fracScale = int64_t(1) << std::max(params.attr_frac_bits, 0);
  return double(params.attr_scale_minus1 + 1) / double(fracScale);
}

static int64_t
scalePowerCodeDomain(const AttributeParameters& params)
{
  int bits = std::max(0, std::min(params.attr_frac_bits, 30));
  return int64_t(1) << bits;
}

static double
scalePowerType3Range(const AttributeDescription& desc)
{
  const double codeMax = double((int64_t(1) << desc.bitdepth) - 1);
  return std::max(scalePowerUniformStep(desc.params) * codeMax, 1e-12);
}

static double
scalePowerClampUnit(double value)
{
  return std::max(0.0, std::min(1.0, value));
}

static double
scalePowerFxpToDouble(int64_t value, int scale)
{
  return double(value) / double(int64_t(1) << scale);
}

static int64_t
scalePowerDoubleToFxp(double value, int scale)
{
  return int64_t(std::llround(value * double(int64_t(1) << scale)));
}

static double
scalePowerFactor(const AttributeParameters& params, int c)
{
  return std::max(
    scalePowerFxpToDouble(
      params.scalePowerParams.powerFactorFxp[c],
      GS_SCALE_POWER_FACTOR_SCALE),
    1e-12);
}

static double
scalePowerScalingFactor(const AttributeParameters& params, int c)
{
  return std::max(
    scalePowerFxpToDouble(
      params.scalePowerParams.scalingFactorFxp[c],
      GS_SCALE_SCALING_FACTOR_SCALE),
    1e-12);
}

static int64_t
scalePowerEffectiveCodeDomain(const AttributeParameters& params, int c)
{
  const int64_t fullDomain = scalePowerCodeDomain(params);
  if (fullDomain <= 1)
    return fullDomain;

  const double scalingFactor =
    std::min(scalePowerScalingFactor(params, c), 1.0);
  const int64_t codeMax = PCCClip(
    int64_t(std::llround(double(fullDomain - 1) * scalingFactor)),
    int64_t(1), fullDomain - 1);
  return codeMax + 1;
}

struct ScalePowerMetric {
  double psnr = 0.0;
  double rate = 0.0;
};

static constexpr size_t kScalePowerSearchSegmentLimit = 50000;

struct ScalePowerSearchContext {
  int numBins = 0;
  double qMax = 0.0;
  double valueMin = 0.0;
  double valueRange = 0.0;
  double sumWeights = 0.0;
  double log10MaxValue = 0.0;
  std::vector<float> xSorted;
};

static ScalePowerSearchContext
scalePowerMakeSearchContext(
  std::vector<float>&& normalizedValuesSorted,
  int bitdepth,
  double valueMin,
  double valueMax)
{
  ScalePowerSearchContext ctxt;
  ctxt.numBins = 1 << bitdepth;
  ctxt.qMax = double(ctxt.numBins - 1);
  ctxt.valueMin = valueMin;
  ctxt.valueRange = valueMax - valueMin;
  ctxt.xSorted = std::move(normalizedValuesSorted);

  double maxValue = ctxt.xSorted.empty()
    ? 0.0
    : std::exp(
        double(ctxt.xSorted.back()) * ctxt.valueRange + ctxt.valueMin);
  ctxt.log10MaxValue = std::log10(std::max(maxValue, 1e-300));

  return ctxt;
}

static std::vector<ScalePowerMetric>
scalePowerMetricsForGridStreaming(
  ScalePowerSearchContext& ctxt, const std::vector<double>& powerGrid)
{
  const int numPowers = int(powerGrid.size());
  std::vector<ScalePowerMetric> metrics(numPowers);
  if (ctxt.xSorted.empty())
    return metrics;

  const int n = int(ctxt.xSorted.size());

  std::vector<double> invPower(numPowers);
  std::vector<std::vector<float>> edges(numPowers);
  for (int p = 0; p < numPowers; p++) {
    invPower[p] = 1.0 / powerGrid[p];
    edges[p].resize(ctxt.numBins - 1);
    for (int q = 0; q + 1 < ctxt.numBins; q++)
      edges[p][q] =
        float(std::pow((double(q) + 0.5) / ctxt.qMax, invPower[p]));
  }

  std::vector<int> binIdx(numPowers, 0);
  std::vector<size_t> binStart(numPowers, 0);
  std::vector<double> startW1(numPowers, 0.0);
  std::vector<double> startW2(numPowers, 0.0);
  std::vector<double> startW3(numPowers, 0.0);
  std::vector<double> weightedError(numPowers, 0.0);

  double W1 = 0.0;
  double W2 = 0.0;
  double W3 = 0.0;

  auto closeBin = [&](int p, size_t end) {
    int q = binIdx[p];
    int count = int(end - binStart[p]);
    if (count > 0) {
      double qNorm =
        scalePowerClampUnit(std::max(double(q) / ctxt.qMax, 1e-7));
      double rec = std::exp(
        std::pow(qNorm, invPower[p]) * ctxt.valueRange + ctxt.valueMin);
      double s1 = W1 - startW1[p];
      double s2 = W2 - startW2[p];
      double s3 = W3 - startW3[p];
      weightedError[p] += s3 - 2.0 * rec * s2 + rec * rec * s1;
    }

    binStart[p] = end;
    startW1[p] = W1;
    startW2[p] = W2;
    startW3[p] = W3;
    binIdx[p]++;
  };

  for (size_t i = 0; i < ctxt.xSorted.size(); i++) {
    float x = ctxt.xSorted[i];
    for (int p = 0; p < numPowers; p++) {
      while (binIdx[p] + 1 < ctxt.numBins && x >= edges[p][binIdx[p]])
        closeBin(p, i);
    }

    double weight = std::exp(double(x) * ctxt.valueRange + ctxt.valueMin);
    W1 += weight;
    W2 += weight * weight;
    W3 += weight * weight * weight;
  }

  for (int p = 0; p < numPowers; p++)
    closeBin(p, ctxt.xSorted.size());

  ctxt.sumWeights = W1;
  if (ctxt.sumWeights <= 0.0)
    return std::vector<ScalePowerMetric>(numPowers);

  for (int p = 0; p < numPowers; p++) {
    double weightedMse =
      std::max(weightedError[p] / ctxt.sumWeights, 1e-12);
    metrics[p].psnr =
      20.0 * (ctxt.log10MaxValue - 0.5 * std::log10(weightedMse));
  }

  return metrics;
}

static std::vector<ScalePowerMetric>
scalePowerMetricsForGrid(
  ScalePowerSearchContext& ctxt, const std::vector<double>& powerGrid)
{
  const int numPowers = int(powerGrid.size());
  std::vector<ScalePowerMetric> metrics(numPowers);
  if (ctxt.xSorted.empty())
    return metrics;

  const size_t sampleCount = ctxt.xSorted.size();
  std::vector<std::vector<size_t>> binEndSampleIdx(numPowers);
  std::vector<size_t> cuts;
  cuts.reserve(size_t(numPowers) * size_t(std::max(ctxt.numBins - 1, 0)) + 2);
  cuts.push_back(0);
  cuts.push_back(sampleCount);

  for (int p = 0; p < numPowers; p++) {
    double invPower = 1.0 / powerGrid[p];
    binEndSampleIdx[p].resize(ctxt.numBins - 1);
    for (int q = 0; q + 1 < ctxt.numBins; q++) {
      float edge =
        float(std::pow((double(q) + 0.5) / ctxt.qMax, invPower));
      size_t end = std::lower_bound(
                     ctxt.xSorted.begin(), ctxt.xSorted.end(), edge)
        - ctxt.xSorted.begin();
      binEndSampleIdx[p][q] = end;
      cuts.push_back(end);
    }
  }

  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

  const size_t segmentCount = cuts.size() - 1;
  if (segmentCount > kScalePowerSearchSegmentLimit)
    return scalePowerMetricsForGridStreaming(ctxt, powerGrid);

  std::vector<std::vector<int>> binEndSegmentIdx(numPowers);
  for (int p = 0; p < numPowers; p++) {
    binEndSegmentIdx[p].resize(ctxt.numBins - 1);
    for (int q = 0; q + 1 < ctxt.numBins; q++) {
      binEndSegmentIdx[p][q] = int(
        std::lower_bound(cuts.begin(), cuts.end(), binEndSampleIdx[p][q])
        - cuts.begin());
    }
  }

  std::vector<double> prefixW1(cuts.size(), 0.0);
  std::vector<double> prefixW2(cuts.size(), 0.0);
  std::vector<double> prefixW3(cuts.size(), 0.0);

  double W1 = 0.0;
  double W2 = 0.0;
  double W3 = 0.0;
  size_t cutIdx = 1;
  for (size_t i = 0; i < sampleCount; i++) {
    while (cutIdx < cuts.size() && cuts[cutIdx] == i) {
      prefixW1[cutIdx] = W1;
      prefixW2[cutIdx] = W2;
      prefixW3[cutIdx] = W3;
      cutIdx++;
    }

    double weight =
      std::exp(double(ctxt.xSorted[i]) * ctxt.valueRange + ctxt.valueMin);
    W1 += weight;
    W2 += weight * weight;
    W3 += weight * weight * weight;
  }

  while (cutIdx < cuts.size()) {
    prefixW1[cutIdx] = W1;
    prefixW2[cutIdx] = W2;
    prefixW3[cutIdx] = W3;
    cutIdx++;
  }

  ctxt.sumWeights = W1;
  if (ctxt.sumWeights <= 0.0)
    return std::vector<ScalePowerMetric>(numPowers);

  const int n = int(sampleCount);
  for (int p = 0; p < numPowers; p++) {
    const double invPower = 1.0 / powerGrid[p];
    double weightedError = 0.0;
    size_t start = 0;

    for (int q = 0; q < ctxt.numBins; q++) {
      size_t end = q + 1 < ctxt.numBins
        ? size_t(binEndSegmentIdx[p][q])
        : segmentCount;
      int count = int(cuts[end] - cuts[start]);
      if (count > 0) {
        double qNorm =
          scalePowerClampUnit(std::max(double(q) / ctxt.qMax, 1e-7));
        double rec = std::exp(
          std::pow(qNorm, invPower) * ctxt.valueRange + ctxt.valueMin);
        double s1 = prefixW1[end] - prefixW1[start];
        double s2 = prefixW2[end] - prefixW2[start];
        double s3 = prefixW3[end] - prefixW3[start];
        weightedError += s3 - 2.0 * rec * s2 + rec * rec * s1;
      }

      start = end;
    }

    double weightedMse = std::max(weightedError / ctxt.sumWeights, 1e-12);
    metrics[p].psnr =
      20.0 * (ctxt.log10MaxValue - 0.5 * std::log10(weightedMse));
  }

  return metrics;
}

static bool
scalePowerSelectBestPower(
  const std::vector<double>& powerGrid,
  const std::vector<ScalePowerMetric>& metrics,
  const ScalePowerMetric& baseline,
  double* bestPower,
  double* rateMin,
  double* ratePtp,
  double* psnrMin,
  double* psnrPtp,
  const double* fixedRateMin = nullptr,
  const double* fixedRatePtp = nullptr,
  const double* fixedPsnrMin = nullptr,
  const double* fixedPsnrPtp = nullptr)
{
  struct Candidate {
    double power;
    double rateGain;
    double psnrDrop;
  };
  std::vector<Candidate> candidates;

  for (size_t i = 0; i < powerGrid.size(); i++) {
    const ScalePowerMetric& metric = metrics[i];
    double rateGain = baseline.rate - metric.rate;
    if (rateGain <= 0.0)
      continue;
    candidates.push_back(
      {powerGrid[i], rateGain, baseline.psnr - metric.psnr});
  }

  if (candidates.empty())
    return false;

  double localRateMin = candidates[0].rateGain;
  double localRateMax = candidates[0].rateGain;
  double localPsnrMin = candidates[0].psnrDrop;
  double localPsnrMax = candidates[0].psnrDrop;
  for (const auto& candidate : candidates) {
    localRateMin = std::min(localRateMin, candidate.rateGain);
    localRateMax = std::max(localRateMax, candidate.rateGain);
    localPsnrMin = std::min(localPsnrMin, candidate.psnrDrop);
    localPsnrMax = std::max(localPsnrMax, candidate.psnrDrop);
  }

  double rMin = fixedRateMin ? *fixedRateMin : localRateMin;
  double rPtp = fixedRatePtp ? *fixedRatePtp : localRateMax - localRateMin;
  double pMin = fixedPsnrMin ? *fixedPsnrMin : localPsnrMin;
  double pPtp = fixedPsnrPtp ? *fixedPsnrPtp : localPsnrMax - localPsnrMin;
  rPtp += 1e-12;
  pPtp += 1e-12;

  const Candidate* best = &candidates[0];
  double bestScore = -std::numeric_limits<double>::infinity();
  for (const auto& candidate : candidates) {
    double score =
      (candidate.rateGain - rMin) / rPtp
      - (candidate.psnrDrop - pMin) / pPtp;
    if (score > bestScore) {
      bestScore = score;
      best = &candidate;
    }
  }

  *bestPower = best->power;
  if (rateMin)
    *rateMin = rMin;
  if (ratePtp)
    *ratePtp = rPtp;
  if (psnrMin)
    *psnrMin = pMin;
  if (psnrPtp)
    *psnrPtp = pPtp;
  return true;
}

static std::vector<double>
scalePowerDiffEntropyForGrid(
  const std::vector<float>& normalizedValuesPointOrder,
  double qMax,
  const std::vector<double>& powerGrid)
{
  const int numPowers = int(powerGrid.size());
  std::vector<double> rate(numPowers, 0.0);
  const size_t n = normalizedValuesPointOrder.size();
  if (n < 2)
    return rate;

  const int qMaxInt = int(std::llround(qMax));
  if (qMaxInt <= 0)
    return rate;

  std::vector<int64_t> diffCounts(size_t(2 * qMaxInt + 1), 0);
  std::vector<int> touchedBins;
  touchedBins.reserve(diffCounts.size());

  auto quantizedCode = [&](size_t i, double invPower) {
    double value = scalePowerClampUnit(double(normalizedValuesPointOrder[i]));
    if (value <= 0.0)
      return 0;
    if (value >= 1.0)
      return qMaxInt;
    return int(std::llround(std::pow(value, invPower) * qMax));
  };

  auto uniformQuantizedCode = [&](size_t i) {
    double value = scalePowerClampUnit(double(normalizedValuesPointOrder[i]));
    if (value <= 0.0)
      return 0;
    if (value >= 1.0)
      return qMaxInt;
    return int(std::llround(value * qMax));
  };

  auto accumulateDiff = [&](int code, int& prevCode, int64_t& total) {
    int diff = code - prevCode;
    prevCode = code;

    int bin = diff + qMaxInt;
    if (!diffCounts[bin])
      touchedBins.push_back(bin);
    diffCounts[bin]++;
    total++;
  };

  for (int p = 0; p < numPowers; p++) {
    int prevCode = powerGrid[p] == 1.0
      ? uniformQuantizedCode(0)
      : quantizedCode(0, 1.0 / powerGrid[p]);
    int64_t total = 0;

    if (powerGrid[p] == 1.0) {
      for (size_t i = 1; i < n; i++)
        accumulateDiff(uniformQuantizedCode(i), prevCode, total);
    } else {
      double invPower = 1.0 / powerGrid[p];
      for (size_t i = 1; i < n; i++)
        accumulateDiff(quantizedCode(i, invPower), prevCode, total);
    }

    double entropy = 0.0;
    for (int bin : touchedBins) {
      double prob = double(diffCounts[bin]) / double(total);
      entropy -= prob * std::log2(prob);
      diffCounts[bin] = 0;
    }
    touchedBins.clear();
    rate[p] = entropy;
  }
  return rate;
}

static double
scalePowerFindOptimalPower(
  std::vector<float>&& normalizedValuesSorted,
  const std::vector<float>& normalizedValuesPointOrder,
  int bitdepth,
  double valueMin,
  double valueMax)
{
  if (normalizedValuesSorted.empty() || bitdepth <= 0 || bitdepth > 16)
    return 1.0;

  const int searchBitdepth = std::min(bitdepth, 10);
  ScalePowerSearchContext ctxt = scalePowerMakeSearchContext(
    std::move(normalizedValuesSorted), searchBitdepth, valueMin, valueMax);

  const int coarseSearchCount = 20;
  std::vector<double> coarseGrid;
  coarseGrid.push_back(1.0);
  for (int i = 0; i < coarseSearchCount; i++)
    coarseGrid.push_back(
      3.5 + (8.0 - 3.5) * double(i) / double(coarseSearchCount - 1));

  std::vector<ScalePowerMetric> coarseMetrics =
    scalePowerMetricsForGrid(ctxt, coarseGrid);
  std::vector<double> coarseDiffRate = scalePowerDiffEntropyForGrid(
    normalizedValuesPointOrder, ctxt.qMax, coarseGrid);
  for (size_t i = 0; i < coarseMetrics.size(); i++)
    coarseMetrics[i].rate = coarseDiffRate[i];

  ScalePowerMetric baseline = coarseMetrics.front();
  coarseGrid.erase(coarseGrid.begin());
  coarseMetrics.erase(coarseMetrics.begin());

  double coarsePower = 3.5;
  double rateMin = 0.0;
  double ratePtp = 1.0;
  double psnrMin = 0.0;
  double psnrPtp = 1.0;
  if (!scalePowerSelectBestPower(
        coarseGrid, coarseMetrics, baseline, &coarsePower, &rateMin,
        &ratePtp, &psnrMin, &psnrPtp)) {
    return 1.0;
  }

  const int fineSearchCount = 10;
  std::vector<double> fineGrid;
  double lo = std::max(3.5, coarsePower - 0.2);
  double hi = std::min(8.0, coarsePower + 0.2);
  for (int i = 0; i < fineSearchCount; i++)
    fineGrid.push_back(
      lo + (hi - lo) * double(i) / double(fineSearchCount - 1));

  std::vector<ScalePowerMetric> fineMetrics =
    scalePowerMetricsForGrid(ctxt, fineGrid);
  std::vector<double> fineDiffRate = scalePowerDiffEntropyForGrid(
    normalizedValuesPointOrder, ctxt.qMax, fineGrid);
  for (size_t i = 0; i < fineMetrics.size(); i++)
    fineMetrics[i].rate = fineDiffRate[i];

  double finePower = coarsePower;
  if (!scalePowerSelectBestPower(
        fineGrid, fineMetrics, baseline, &finePower, nullptr, nullptr,
        nullptr, nullptr, &rateMin, &ratePtp, &psnrMin, &psnrPtp)) {
    finePower = coarsePower;
  }

  finePower = std::round(finePower * 100.0) / 100.0;

  std::vector<ScalePowerMetric> winnerMetric =
    scalePowerMetricsForGrid(ctxt, std::vector<double>{finePower});
  if (baseline.psnr - winnerMetric[0].psnr > 1.0)
    return 1.0;

  return finePower;
}

//============================================================================

void
setScalePowerTransformParams(
  AttributeDescription& desc, const std::vector<double>& inputAttrBuf,
  const std::vector<int>& indices, int pointCount, int64_t dCode)
{
  if (desc.attributeLabel != KnownAttributeLabel::kScale)
    return;

  dCode = std::max<int64_t>(1, dCode);

  double valueMin = double(desc.params.attr_offset);
  double valueMax = valueMin + double(dCode);
  double valueRange = std::max(valueMax - valueMin, 1e-12);

  constexpr size_t kMaxPointOrderSamples = size_t(1) << 14;
  const size_t numComp = std::max<size_t>(1, indices.size());
  const size_t pointOrderRunLen =
    std::max<size_t>(1, kMaxPointOrderSamples / numComp);

  std::vector<float> normalizedSamples;
  std::vector<float> scaleSamplesPointOrder;
  normalizedSamples.reserve(size_t(pointCount) * indices.size());
  scaleSamplesPointOrder.reserve(pointOrderRunLen * indices.size());
  for (int c = 0; c < int(indices.size()); c++) {
    const double* buf = &inputAttrBuf[indices[c] * pointCount];
    size_t runRemaining = pointOrderRunLen;
    for (int i = 0; i < pointCount; i++) {
      double value = buf[i];
      if (!std::isfinite(value))
        continue;

      float raw = float(value);
      float norm =
        float(scalePowerClampUnit((double(raw) - valueMin) / valueRange));
      normalizedSamples.push_back(norm);
      if (runRemaining) {
        scaleSamplesPointOrder.push_back(norm);
        runRemaining--;
      }
    }
  }

  std::sort(normalizedSamples.begin(), normalizedSamples.end());

  if (!normalizedSamples.empty()) {
    size_t trim = size_t(double(normalizedSamples.size()) * 0.0005);
    size_t first = trim;
    size_t last = normalizedSamples.size() - trim;
    if (first >= last) {
      first = 0;
      last = normalizedSamples.size();
    }

    size_t out = 0;
    for (size_t i = first; i < last; i++)
      normalizedSamples[out++] = normalizedSamples[i];
    normalizedSamples.resize(out);
  }

  double power = scalePowerFindOptimalPower(
    std::move(normalizedSamples), scaleSamplesPointOrder, desc.bitdepth,
    valueMin, valueMax);

  auto& scaleParams = desc.params.scalePowerParams;
  for (int c = 0; c <= desc.attr_num_dimensions_minus1; c++)
    scaleParams.powerFactorFxp[c] =
      scalePowerDoubleToFxp(power, GS_SCALE_POWER_FACTOR_SCALE);

  ensureScalePowerTransformParams(desc);
}

void
ensureScalePowerTransformParams(AttributeDescription& desc)
{
  if (desc.attributeLabel != KnownAttributeLabel::kScale)
    return;

  auto& params = desc.params;
  auto& scaleParams = params.scalePowerParams;
  double uniformStep = scalePowerUniformStep(params);
  if (!std::isfinite(uniformStep) || uniformStep <= 0.0)
    uniformStep = 1.0;
  int64_t uniformStepFxp =
    scalePowerDoubleToFxp(uniformStep, GS_SCALE_STEP_SIZE_SCALE);

  for (int c = 0; c <= desc.attr_num_dimensions_minus1; c++) {
    auto& transitionPoints = scaleParams.quantRegionTransitionPoint[c];
    auto& stepSizes = scaleParams.regionStepSizeFxp[c];

    double range = scalePowerType3Range(desc);
    double scalingFactor = scalePowerScalingFactor(params, c);
    double threshold0 =
      scalePowerClampUnit((-10.0 - double(params.attr_offset)) / range);
    double threshold1 =
      scalePowerClampUnit((-8.0 - double(params.attr_offset)) / range);
    if (threshold1 < threshold0)
      std::swap(threshold0, threshold1);

    double yMax = std::min(scalingFactor, 1.0);
    int64_t codeDomain = scalePowerEffectiveCodeDomain(params, c);

    if (codeDomain < 4 || threshold1 <= 0.0) {
      scaleParams.numQuantRegionsMinus1[c] = 0;
      scaleParams.numBitsForTransitionPoint[c] = params.attr_frac_bits;
      transitionPoints.clear();
      stepSizes.assign(1, uniformStepFxp);
      continue;
    }

    double power = scalePowerFactor(params, c);
    if (power <= 1.0) {
      int transition0 = int(PCCClip(
        int64_t(std::llround(double(codeDomain) * threshold0)),
        int64_t(1), codeDomain - 2));
      int transition1 = int(PCCClip(
        int64_t(std::llround(double(codeDomain) * threshold1)),
        int64_t(transition0 + 1), codeDomain - 1));

      scaleParams.numQuantRegionsMinus1[c] = 2;
      scaleParams.numBitsForTransitionPoint[c] = params.attr_frac_bits;
      transitionPoints.assign({transition0, transition1});
      stepSizes.assign(3, uniformStepFxp);
      continue;
    }

    double y0 = scalingFactor * std::pow(threshold0, power);
    double y1 = scalingFactor * std::pow(threshold1, power);
    y0 = std::min(y0, yMax);
    y1 = std::min(y1, yMax);
    double dy0 = y0;
    double dy1 = y1 - y0;
    double dy2 = yMax - y1;

    if (dy0 <= 0.0 || dy1 <= 0.0 || dy2 <= 0.0) {
      scaleParams.numQuantRegionsMinus1[c] = 0;
      scaleParams.numBitsForTransitionPoint[c] = params.attr_frac_bits;
      transitionPoints.clear();
      stepSizes.assign(1, uniformStepFxp);
      continue;
    }

    const double ratio0 = 24.0;
    const double ratio1 = 4.0;
    const double ratio2 = 1.0;
    double w0 = dy0 / ratio0;
    double w1 = dy1 / ratio1;
    double w2 = dy2 / ratio2;
    double totalWeight = w0 + w1 + w2;
    if (totalWeight <= 0.0) {
      scaleParams.numQuantRegionsMinus1[c] = 0;
      scaleParams.numBitsForTransitionPoint[c] = params.attr_frac_bits;
      transitionPoints.clear();
      stepSizes.assign(1, uniformStepFxp);
      continue;
    }

    int transition0 = int(std::llround(codeDomain * w0 / totalWeight));
    int transition1 =
      int(std::llround(codeDomain * (w0 + w1) / totalWeight));
    transition0 = int(PCCClip(
      int64_t(transition0), int64_t(1), codeDomain - 2));
    transition1 = int(PCCClip(
      int64_t(transition1), int64_t(transition0 + 1), codeDomain - 1));

    int64_t count0 = transition0;
    int64_t count1 = transition1 - transition0;
    int64_t count2 = codeDomain - transition1;
    double step0 = dy0 * range / double(count0);
    double step1 = dy1 * range / double(count1);
    double step2 = dy2 * range / double(count2);

    if (step0 <= 0.0 || step1 <= 0.0 || step2 <= 0.0) {
      scaleParams.numQuantRegionsMinus1[c] = 0;
      scaleParams.numBitsForTransitionPoint[c] = params.attr_frac_bits;
      transitionPoints.clear();
      stepSizes.assign(1, uniformStepFxp);
      continue;
    }

    scaleParams.numQuantRegionsMinus1[c] = 2;
    scaleParams.numBitsForTransitionPoint[c] = params.attr_frac_bits;
    transitionPoints.assign({transition0, transition1});
    stepSizes.assign({
      scalePowerDoubleToFxp(step0, GS_SCALE_STEP_SIZE_SCALE),
      scalePowerDoubleToFxp(step1, GS_SCALE_STEP_SIZE_SCALE),
      scalePowerDoubleToFxp(step2, GS_SCALE_STEP_SIZE_SCALE),
    });
  }
}

struct ScalePowerRuntimeContext {
  const std::vector<int>* transitionPoints = nullptr;
  int numRegionsMinus1 = 0;
  int maxTransitionPoint = 0;
  int64_t attrScale = 0;
  int64_t fracScale = 0;
  double attrOffset = 0.0;
  double range = 1.0;
  double scalingFactor = 1.0;
  double powerFactor = 1.0;
  double invPowerFactor = 1.0;
  double uniformStep = 1.0;
  bool useRegions = false;
  bool uniformIdentityRegions = false;
  std::vector<double> offsetReal;
  std::vector<double> stepSizes;
};

static double
scalePowerPow(double base, double exponent)
{
  return exponent == 1.0 ? base : std::pow(base, exponent);
}

static ScalePowerRuntimeContext
scalePowerMakeRuntimeContext(const AttributeDescription& desc, int c)
{
  const auto& params = desc.params;
  const auto& scaleParams = params.scalePowerParams;

  ScalePowerRuntimeContext ctxt;
  ctxt.numRegionsMinus1 = scaleParams.numQuantRegionsMinus1[c];
  ctxt.transitionPoints = &scaleParams.quantRegionTransitionPoint[c];
  ctxt.attrScale = params.attr_scale_minus1 + 1;
  ctxt.fracScale = int64_t(1) << params.attr_frac_bits;
  ctxt.attrOffset = double(params.attr_offset);
  ctxt.range = scalePowerType3Range(desc);
  ctxt.scalingFactor = scalePowerScalingFactor(params, c);
  ctxt.powerFactor = scalePowerFactor(params, c);
  ctxt.invPowerFactor = 1.0 / ctxt.powerFactor;
  ctxt.uniformStep = scalePowerUniformStep(params);

  const auto& stepSizesFxp = scaleParams.regionStepSizeFxp[c];
  ctxt.useRegions =
    ctxt.numRegionsMinus1 > 0
    && int(ctxt.transitionPoints->size()) == ctxt.numRegionsMinus1
    && int(stepSizesFxp.size()) == ctxt.numRegionsMinus1 + 1;

  if (!ctxt.useRegions)
    return ctxt;

  ctxt.maxTransitionPoint = int(scalePowerEffectiveCodeDomain(params, c));
  ctxt.offsetReal.resize(ctxt.numRegionsMinus1 + 1);
  ctxt.stepSizes.resize(ctxt.numRegionsMinus1 + 1);
  ctxt.offsetReal[0] = 0.0;

  for (int m = 0; m < ctxt.numRegionsMinus1 + 1; m++) {
    ctxt.stepSizes[m] =
      scalePowerFxpToDouble(stepSizesFxp[m], GS_SCALE_STEP_SIZE_SCALE);
    if (m == 0)
      continue;

    int range = m == 1
      ? (*ctxt.transitionPoints)[m - 1]
      : (*ctxt.transitionPoints)[m - 1] - (*ctxt.transitionPoints)[m - 2];
    ctxt.offsetReal[m] =
      ctxt.offsetReal[m - 1] + range * ctxt.stepSizes[m - 1];
  }

  ctxt.uniformIdentityRegions =
    ctxt.powerFactor <= 1.0 && ctxt.scalingFactor == 1.0
    && !ctxt.stepSizes.empty();
  for (double stepSize : ctxt.stepSizes)
    ctxt.uniformIdentityRegions =
      ctxt.uniformIdentityRegions && stepSize == ctxt.stepSizes[0];

  return ctxt;
}

static double
scalePowerForwardTransform(
  const ScalePowerRuntimeContext& ctxt, double value)
{
  double normalized =
    scalePowerClampUnit((value - ctxt.attrOffset) / ctxt.range);
  double transformedNorm = scalePowerClampUnit(
    ctxt.scalingFactor * scalePowerPow(normalized, ctxt.powerFactor));
  double transformed = ctxt.range * transformedNorm;
  return std::isfinite(transformed) ? transformed : 0.0;
}

static double
scalePowerInverseTransform(
  const ScalePowerRuntimeContext& ctxt, double value)
{
  double normalized =
    scalePowerClampUnit(value / (ctxt.range * ctxt.scalingFactor));
  double transformed =
    ctxt.attrOffset + ctxt.range * scalePowerPow(normalized, ctxt.invPowerFactor);
  return std::isfinite(transformed) ? transformed : ctxt.attrOffset;
}

static double
scalePowerInverseRegionValue(const ScalePowerRuntimeContext& ctxt, int code)
{
  double value = double(code) * ctxt.uniformStep;

  if (!ctxt.useRegions)
    return value;

  if (ctxt.uniformIdentityRegions)
    return double(code) * ctxt.stepSizes[0];

  for (int m = 0; m < ctxt.numRegionsMinus1 + 1; m++) {
    int transitionPoint = m == ctxt.numRegionsMinus1
      ? ctxt.maxTransitionPoint
      : (*ctxt.transitionPoints)[m];
    if (code < transitionPoint) {
      int startInt = m == 0 ? 0 : (*ctxt.transitionPoints)[m - 1];
      return (code - startInt) * ctxt.stepSizes[m] + ctxt.offsetReal[m];
    }
  }

  return value;
}

static int
scalePowerQuantizeUniform(
  const ScalePowerRuntimeContext& ctxt, double target, int64_t clipMax)
{
  if (ctxt.attrScale <= 0)
    return 0;

  return int(PCCClip(
    std::llround(target * double(ctxt.fracScale) / double(ctxt.attrScale)),
    int64_t(0), clipMax));
}

static int
scalePowerQuantizeTransform(
  const ScalePowerRuntimeContext& ctxt, double value, int64_t clipMax)
{
  double target = scalePowerForwardTransform(ctxt, value);

  if (ctxt.useRegions) {
    if (ctxt.uniformIdentityRegions) {
      int64_t code = std::llround(target / ctxt.stepSizes[0]);
      return int(PCCClip(
        code, int64_t(0),
        std::min<int64_t>(clipMax, int64_t(ctxt.maxTransitionPoint - 1))));
    }

    for (int m = 0; m < ctxt.numRegionsMinus1 + 1; m++) {
      int startInt = m == 0 ? 0 : (*ctxt.transitionPoints)[m - 1];
      int transitionPoint = m == ctxt.numRegionsMinus1
        ? ctxt.maxTransitionPoint
        : (*ctxt.transitionPoints)[m];
      double stepSize = ctxt.stepSizes[m];
      if (stepSize <= 0.0)
        break;

      double regionEnd =
        ctxt.offsetReal[m] + (transitionPoint - startInt) * stepSize;
      if (target < regionEnd || m == ctxt.numRegionsMinus1) {
        int64_t code =
          startInt + std::llround((target - ctxt.offsetReal[m]) / stepSize);
        return int(PCCClip(
          code, int64_t(0),
          std::min<int64_t>(clipMax, int64_t(ctxt.maxTransitionPoint - 1))));
      }
    }
  }

  return scalePowerQuantizeUniform(ctxt, target, clipMax);
}

void
inverseScalePowerTransform(
  const AttributeDescription& desc, const Vec3<attr_t>* src, Vec3<float>* dst,
  int pointCount)
{
  if (!src || !dst || pointCount <= 0)
    return;

  for (int c = 0; c <= desc.attr_num_dimensions_minus1 && c < 3; c++) {
    ScalePowerRuntimeContext ctxt = scalePowerMakeRuntimeContext(desc, c);
    const bool useCodeCache = ctxt.powerFactor != 1.0 && pointCount > 32768;
    std::vector<float> codeCache;
    std::vector<uint8_t> codeSeen;
    if (useCodeCache) {
      size_t cacheSize = size_t(std::numeric_limits<attr_t>::max()) + 1;
      codeCache.resize(cacheSize);
      codeSeen.resize(cacheSize, 0);
    }

    for (int i = 0; i < pointCount; i++) {
      int code = int(src[i][c]);
      if (useCodeCache) {
        size_t cacheIdx = size_t(code);
        if (!codeSeen[cacheIdx]) {
          codeCache[cacheIdx] = float(scalePowerInverseTransform(
            ctxt, scalePowerInverseRegionValue(ctxt, code)));
          codeSeen[cacheIdx] = 1;
        }
        dst[i][c] = codeCache[cacheIdx];
      } else {
        dst[i][c] = float(scalePowerInverseTransform(
          ctxt, scalePowerInverseRegionValue(ctxt, code)));
      }
    }
  }
}

int
quantizeScalePowerTransform(
  const AttributeDescription& desc, double value, int c, int64_t clipMax)
{
  ScalePowerRuntimeContext ctxt = scalePowerMakeRuntimeContext(desc, c);
  return scalePowerQuantizeTransform(ctxt, value, clipMax);
}

void
quantizeScalePowerTransform(
  const AttributeDescription& desc, const double* src, int c, int64_t clipMax,
  std::vector<Vec3<attr_t>>& dst, int pointCount)
{
  if (!src || pointCount <= 0)
    return;

  ScalePowerRuntimeContext ctxt = scalePowerMakeRuntimeContext(desc, c);
  for (int i = 0; i < pointCount; i++)
    dst[i][c] = attr_t(scalePowerQuantizeTransform(ctxt, src[i], clipMax));
}

//============================================================================

}  // namespace pcc
