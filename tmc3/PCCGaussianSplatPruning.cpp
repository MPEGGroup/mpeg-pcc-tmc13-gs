#include "PCCGaussianSplatPruning.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace pcc {

PCCPointSet3 getPartition(const PCCPointSet3& src, const std::vector<int32_t>& indexes);

bool
pruneGaussianSplatsByScore(PCCPointSet3& cloud, double pruningPercentage)
{
  if (pruningPercentage == 0.)
    return true;

  if (pruningPercentage < 0. || pruningPercentage >= 100.) {
    std::cerr << "[ERROR] pruningPercentage must be in range [0,100), got "
              << pruningPercentage << "\n";
    return false;
  }

  if (!cloud.hasScores()) {
    std::cerr
      << "[ERROR] score-based pruning requires input 'score' attribute. "
         "Provide an input PLY containing per-point 'score'.\n";
    return false;
  }

  const size_t numInputPoints = cloud.getPointCount();
  const size_t numToPrune = static_cast<size_t>(
    std::floor(double(numInputPoints) * pruningPercentage / 100.0));

  if (numToPrune == 0)
    return true;

  const size_t numToKeep = numInputPoints - numToPrune;
  std::vector<int32_t> rankedIndexes(numInputPoints);
  std::iota(rankedIndexes.begin(), rankedIndexes.end(), 0);

  std::stable_sort(
    rankedIndexes.begin(), rankedIndexes.end(),
    [&cloud](int32_t a, int32_t b) {
      const float sa = cloud.getScore(a);
      const float sb = cloud.getScore(b);
      if (sa != sb)
        return sa > sb;
      return a < b;
    });

  rankedIndexes.resize(numToKeep);
  std::sort(rankedIndexes.begin(), rankedIndexes.end());
  cloud = getPartition(cloud, rankedIndexes);
  return true;
}

}  // namespace pcc