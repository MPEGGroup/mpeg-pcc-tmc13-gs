/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such rights
 * are granted under this licence.
 */

#pragma once

#include "PCCPointSet.h"

namespace pcc {

bool pruneGaussianSplatsByScore(PCCPointSet3& cloud, double pruningPercentage);

}  // namespace pcc