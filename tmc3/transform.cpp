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
#include "PCCSymEigen.h"

#include <numeric>
#include <string>
#include <unordered_map>

#include <fstream>
#include <string>
#include <vector>

namespace pcc {

//============================================================================
/// Initialize SH 2D transform parameters
void initializeShTransformParams(AttributeParameters& params) {
  if (!params.gsdTransform2dPresentFlag_) {
    params.shTrTypeCol_.clear();
    params.shTrGrpRow_.clear();
    params.shTrTypeRow_.clear();
    params.shKltFlag_.clear();
  }
  else if (params.dataAdaptiveShTransformEnabled_) {
    // parameter for group PCA m75790
    params.shTrTypeCol_.assign(16, 0);
    params.shTrGrpRow_ = { { 1, 5, 11 }, { 3, 7, 13 }, { 4, 10 }, { 8, 14 }, { 2, 6, 12 }, { 9, 15 } };
    params.shTrTypeRow_.assign(6, 0);
    params.shKltFlag_.assign(6, 1);
    if (params.gsdPcaDcFlag_) {
      params.shTrGrpRow_.push_back({0});
      params.shTrTypeRow_.push_back(0);
      params.shKltFlag_.push_back(1);
    }
  }
  else {
    // parameter for 2D transform m75038
    params.shTrTypeCol_.assign(16, 0);
    params.shTrTypeCol_[0] = 0;
    params.shTrGrpRow_ = { { 1, 5, 11 }, { 3, 7, 13 }, { 4, 10 }, { 8, 14 }, { 2, 6, 12 }, { 9, 15 } };
    params.shTrTypeRow_.assign(6, 2);
    params.shKltFlag_.assign(6, 0);
  }
}

bool computeShKlt(const std::vector<GSPoint>& gsBuf, size_t pointCount,
  const std::vector<int>& groupidx, PCCKltParameters& klt) {
  if ( pointCount == 0 ) return false;
  const size_t size     = groupidx.size();
  const size_t fullSize = size * 3;

  klt.eigvec.assign( fullSize, std::vector<double>( fullSize, 0.0 ) );
  klt.mean.assign( fullSize, 0.0 );
  klt.sigma.assign( fullSize, 1.0 );
  klt.eigvecFxp.assign(fullSize, std::vector<int64_t>(fullSize, 0));
  klt.meanFxp.assign(fullSize, 0);
  klt.sigmaFxp.assign(fullSize, 1 << KLT_SIGMA_SCALE);

  for ( size_t p = 0; p < pointCount; ++p ) {
    for ( int i = 0; i < size; ++i ) {
    for ( int ch = 0; ch < 3; ++ch ) { klt.mean[i * 3 + ch] += gsBuf[p].sh[groupidx[i]][ch]; }
    }
  }
  for ( int j = 0; j < fullSize; ++j ) { klt.mean[j] /= (double)pointCount; }

  std::vector<double> sigma( fullSize, 0.0 );
  std::vector<double> centered( fullSize );

  for ( size_t p = 0; p < pointCount; ++p ) {
    for ( int j = 0; j < fullSize; ++j ) {
    centered[j] = (double)gsBuf[p].sh[groupidx[j / 3]][j % 3] - klt.mean[j];
    sigma[j] += centered[j] * centered[j];
    }
  }

  for ( int j = 0; j < fullSize; ++j ) {
    sigma[j] = std::sqrt( sigma[j] / (double)pointCount );
    if ( sigma[j] < 1e-10 ) sigma[j] = 1.0;
    klt.sigma[j] = sigma[j];
  }

  std::vector<double> cov( fullSize * fullSize, 0.0 );
  std::vector<double> normalized( fullSize );
  for ( size_t p = 0; p < pointCount; ++p ) {
    for ( int j = 0; j < fullSize; ++j ) { normalized[j] = ( (double)gsBuf[p].sh[groupidx[j / 3]][j % 3] - klt.mean[j] ) / klt.sigma[j]; }
    for ( int r = 0; r < fullSize; ++r ) {
    for ( int c = 0; c < fullSize; ++c ) { cov[r + c * fullSize] += normalized[r] * normalized[c]; }
    }
  }

  double div = double( pointCount > 1 ? pointCount - 1 : 1 );
  for ( int i = 0; i < fullSize * fullSize; ++i ) { cov[i] /= div; }
  // Eigen Decomposition using built-in Jacobi method
  pcc::SymEigenResult eigResult;
  if ( !pcc::symEigenDecomposeJacobi( cov, fullSize, eigResult ) ) {
    std::cerr << "symEigenDecomposeJacobi failed for matrix size " << fullSize << std::endl;
    return false;
  }

  // Sort in descending order (Jacobi returns ascending)
  const auto& w = eigResult.eigenvalues;
  const auto& V = eigResult.eigenvectors;
  int n = (int)fullSize;
  std::vector<std::pair<double, int>> eigPairs;
  for ( int i = 0; i < n; ++i ) { eigPairs.push_back( { w[i], i } ); }
  std::sort( eigPairs.rbegin(), eigPairs.rend() );

  // Jacobi returns eigenvectors in columns of V (column-major).
  for ( int j = 0; j < n; ++j ) {
    int srcIdx = eigPairs[j].second;
    for ( int i = 0; i < n; ++i ) {
      klt.eigvec[j][i] = V[i + srcIdx * n];
    }
  }

  // Populate fixed-point representations
  for ( int j = 0; j < n; ++j ) {
    for ( int i = 0; i < n; ++i ) {
      klt.eigvecFxp[j][i] = (int64_t)round(klt.eigvec[j][i] * (1 << KLT_SCALE));
    }
  }
  for ( int j = 0; j < fullSize; ++j ) {
    klt.meanFxp[j] = (int64_t)round(klt.mean[j] * (1 << KLT_MEAN_SCALE));
    klt.sigmaFxp[j] = (int64_t)round(klt.sigma[j] * (1 << KLT_SIGMA_SCALE));
  }

  return true;
}

  

bool sphericalHarmonicsInvTransform2D(bool witSHCoefficients, std::vector<GSPoint>& gsBuf,
  size_t pointCount, AttributeParameters& params) {
  if (!witSHCoefficients) return false;

  const double itr3x3Flp[3][3][3] = {
    // BT.709
    { {1.0, 0.0, +1.5748}, {1.0, -0.1873, -0.4681}, {1.0, +1.8556, 0.0} },
    // 601
    { { 1.00000, 0.00000, 1.40200 }, { 1.00000, -0.34000, -0.71400 }, { 1.00000, 1.77200, 0.00000 } },
    // DCT
    { { 1.00000, 1.22474, 0.70711 }, { 1.00000, 0.00000, -1.41421 }, { 1.00000, -1.22474, 0.70711 } }
  };
  const double itr2x2Flp[3][2][2] = {
      { { 1.0, 0.0 }, { 0.0, 1.0 } },
      { { 1.0, 0.0 }, { 0.0, 1.0 } },
      { { 1.0, 1.0 }, { 1.0, -1.0 } }
  };
  double dct_scale = (double)(1 << DCT_SCALE);
  int64_t itr3x3[3][3][3];
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      for (auto k = 0; k < 3; ++k) {
        itr3x3[i][j][k] = (int)round(itr3x3Flp[i][j][k] * dct_scale);
      }
    }
  }
  int64_t itr2x2[3][2][2];
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 2; ++j) {
      for (auto k = 0; k < 2; ++k) {
        itr2x2[i][j][k] = (int)round(itr2x2Flp[i][j][k] * dct_scale);
      }
    }
  }

  std::vector<int> colTr(16, 1);
  for (size_t i = 0; i < params.shTrGrpRow_.size(); ++i) {
    if (params.shKltFlag_[i] == 1) {
      for (int idx : params.shTrGrpRow_[i]) { colTr[idx] = 0; }
    }
  }

  int   ki = 0;
  for (size_t i = 0; i < params.shTrGrpRow_.size(); ++i) {
    const auto& groupidx = params.shTrGrpRow_[i];
    int         size = groupidx.size();
    if (params.shKltFlag_[i] == 1) {
      auto& klt = params.shKltParams_[ki++];
      int         fullSize = size * 3;
      bool        isDcGroup = (groupidx.size() == 1 && groupidx[0] == 0);

      if (isDcGroup) {
        // DC PCA inverse: x = eigvecFxp^T * y >> KLT_SCALE + meanFxp
        const int dcSize = 3;
        std::vector<int64_t> Ydc(dcSize);
        for (size_t p = 0; p < pointCount; ++p) {
          auto& sh = gsBuf[p].shFxp;
          for (int j = 0; j < dcSize; ++j) Ydc[j] = sh[0][j];
          for (int i_ = 0; i_ < dcSize; i_++) {
            int64_t x_norm = 0;
            for (int j = 0; j < dcSize; j++)
              x_norm += klt.eigvecFxp[j][i_] * Ydc[j];
            x_norm = x_norm >> KLT_SCALE;
            sh[0][i_] = x_norm + klt.meanFxp[i_];
          }
        }
      } else {
        std::vector<int64_t> Y(fullSize);
        for (size_t p = 0; p < pointCount; ++p) {
          auto& sh = gsBuf[p].shFxp;
          for (int j = 0; j < fullSize; j++) { Y[j] = sh[groupidx[j / 3]][j % 3]; }

          for (int i_ = 0; i_ < fullSize; i_++) {
            int64_t x_norm = 0;
            for (int j = 0; j < fullSize; j++) { x_norm += klt.eigvecFxp[j][i_] * Y[j] * 35; }
            x_norm = x_norm >> (KLT_SCALE + 1);
            x_norm = (x_norm * klt.sigmaFxp[i_]) >> KLT_SIGMA_SCALE;
            sh[groupidx[i_ / 3]][i_ % 3] = x_norm + klt.meanFxp[i_];
          }
        }
      }
    }
    else {
      int type = params.shTrTypeRow_[i];
      if (groupidx.size() == 3) {
        for (size_t p = 0; p < pointCount; ++p) {
          gsBuf[p].fxpScale[groupidx[0]] += DCT_SCALE;
          gsBuf[p].fxpScale[groupidx[1]] += DCT_SCALE;
          gsBuf[p].fxpScale[groupidx[2]] += DCT_SCALE;

          auto& sh = gsBuf[p].shFxp;
          for (int c = 0; c < 3; ++c) {
            auto v0 = sh[groupidx[0]][c], v1 = sh[groupidx[1]][c], v2 = sh[groupidx[2]][c];
            sh[groupidx[0]][c] = itr3x3[type][0][0] * v0 + itr3x3[type][0][1] * v1 + itr3x3[type][0][2] * v2;
            sh[groupidx[1]][c] = itr3x3[type][1][0] * v0 + itr3x3[type][1][1] * v1 + itr3x3[type][1][2] * v2;
            sh[groupidx[2]][c] = itr3x3[type][2][0] * v0 + itr3x3[type][2][1] * v1 + itr3x3[type][2][2] * v2;
          }
        }
      }
      else if (groupidx.size() == 2) {
        for (size_t p = 0; p < pointCount; ++p) {
          gsBuf[p].fxpScale[groupidx[0]] += DCT_SCALE;
          gsBuf[p].fxpScale[groupidx[1]] += DCT_SCALE;

          auto& sh = gsBuf[p].shFxp;
          for (int c = 0; c < 3; ++c) {
            auto v0 = sh[groupidx[0]][c], v1 = sh[groupidx[1]][c];
            sh[groupidx[0]][c] = itr2x2[type][0][0] * v0 + itr2x2[type][0][1] * v1;
            sh[groupidx[1]][c] = itr2x2[type][1][0] * v0 + itr2x2[type][1][1] * v1;
          }
        }
      }
    }
  }

  // Inverse Column transform
  for (size_t p = 0; p < pointCount; ++p) {
    auto& sh = gsBuf[p].shFxp;
    auto& shFlp = gsBuf[p].sh;
    for (size_t i = 0; i < 16; ++i) {
      if (colTr[i] == 1) {
        gsBuf[p].fxpScale[i] += DCT_SCALE;
        int    type = params.shTrTypeCol_[i];
        auto y = sh[i][0], u = sh[i][1], v = sh[i][2];
        sh[i][0] = itr3x3[type][0][0] * y + itr3x3[type][0][1] * u + itr3x3[type][0][2] * v;
        sh[i][1] = itr3x3[type][1][0] * y + itr3x3[type][1][1] * u + itr3x3[type][1][2] * v;
        sh[i][2] = itr3x3[type][2][0] * y + itr3x3[type][2][1] * u + itr3x3[type][2][2] * v;
      }
      int shift = gsBuf[p].fxpScale[i] - SH_SCALE;
      sh[i][0] = sh[i][0] >> shift;
      sh[i][1] = sh[i][1] >> shift;
      sh[i][2] = sh[i][2] >> shift;
    }
  }
  return true;
}

bool sphericalHarmonicsInvTransform(bool witSHCoefficients, std::vector<GSPoint>& gsBuf,
  size_t pointCount, AttributeParameters& params, bool convertColourspace){
  if ( params.gsdTransform2dPresentFlag_ ) {
    return sphericalHarmonicsInvTransform2D(witSHCoefficients, gsBuf, pointCount, params ); }
  else{
    if (!convertColourspace) return true;

    for(size_t p = 0; p < pointCount; p++) {
      auto& sh = gsBuf[p].shFxp;
      for (size_t m = 0; m < 16; m++)
        convertToRgb(params.cicp_matrix_coefficients_idx, sh[m]);
    }
  }
  return true;
}

//============================================================================
// SHBufferView-based transform implementations (zero-copy Option 3)
//============================================================================

static bool computeShKltFromView(SHBufferView& shView, size_t pointCount,
  const std::vector<int>& groupidx, PCCKltParameters& klt) {
  if (pointCount == 0) return false;
  const size_t size = groupidx.size();
  const size_t fullSize = size * 3;

  klt.eigvec.assign(fullSize, std::vector<double>(fullSize, 0.0));
  klt.mean.assign(fullSize, 0.0);
  klt.sigma.assign(fullSize, 1.0);
  klt.eigvecFxp.assign(fullSize, std::vector<int64_t>(fullSize, 0));
  klt.meanFxp.assign(fullSize, 0);
  klt.sigmaFxp.assign(fullSize, 1 << KLT_SIGMA_SCALE);

  // Compute mean
  for (size_t p = 0; p < pointCount; ++p) {
    for (int i = 0; i < size; ++i) {
      for (int ch = 0; ch < 3; ++ch) {
        klt.mean[i * 3 + ch] += shView.sh(p, groupidx[i], ch);
      }
    }
  }
  for (int j = 0; j < fullSize; ++j) { klt.mean[j] /= (double)pointCount; }

  // Compute sigma
  std::vector<double> sigma(fullSize, 0.0);
  std::vector<double> centered(fullSize);
  for (size_t p = 0; p < pointCount; ++p) {
    for (int j = 0; j < fullSize; ++j) {
      centered[j] = shView.sh(p, groupidx[j / 3], j % 3) - klt.mean[j];
      sigma[j] += centered[j] * centered[j];
    }
  }
  for (int j = 0; j < fullSize; ++j) {
    sigma[j] = std::sqrt(sigma[j] / (double)pointCount);
    if (sigma[j] < 1e-10) sigma[j] = 1.0;
    klt.sigma[j] = sigma[j];
  }

  // Compute covariance
  std::vector<double> cov(fullSize * fullSize, 0.0);
  std::vector<double> normalized(fullSize);
  for (size_t p = 0; p < pointCount; ++p) {
    for (int j = 0; j < fullSize; ++j) {
      normalized[j] = (shView.sh(p, groupidx[j / 3], j % 3) - klt.mean[j]) / klt.sigma[j];
    }
    for (int r = 0; r < fullSize; ++r) {
      for (int c = 0; c < fullSize; ++c) {
        cov[r + c * fullSize] += normalized[r] * normalized[c];
      }
    }
  }

  double div = double(pointCount > 1 ? pointCount - 1 : 1);
  for (int i = 0; i < fullSize * fullSize; ++i) { cov[i] /= div; }

  // Eigen decomposition
  pcc::SymEigenResult eigResult;
  if (!pcc::symEigenDecomposeJacobi(cov, fullSize, eigResult)) {
    std::cerr << "symEigenDecomposeJacobi failed for matrix size " << fullSize << std::endl;
    return false;
  }

  // Sort in descending order
  const auto& w = eigResult.eigenvalues;
  const auto& V = eigResult.eigenvectors;
  int n = (int)fullSize;
  std::vector<std::pair<double, int>> eigPairs;
  for (int i = 0; i < n; ++i) { eigPairs.push_back({w[i], i}); }
  std::sort(eigPairs.rbegin(), eigPairs.rend());

  // Store eigenvectors
  for (int j = 0; j < n; ++j) {
    int srcIdx = eigPairs[j].second;
    for (int i = 0; i < n; ++i) {
      klt.eigvec[j][i] = V[i + srcIdx * n];
    }
  }

  // Fixed-point representations
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      klt.eigvecFxp[j][i] = (int64_t)round(klt.eigvec[j][i] * (1 << KLT_SCALE));
    }
  }
  for (int j = 0; j < fullSize; ++j) {
    klt.meanFxp[j] = (int64_t)round(klt.mean[j] * (1 << KLT_MEAN_SCALE));
    klt.sigmaFxp[j] = (int64_t)round(klt.sigma[j] * (1 << KLT_SIGMA_SCALE));
  }

  return true;
}

bool sphericalHarmonicsTransform2D(bool witSHCoefficients, SHBufferView& shView,
  size_t pointCount, AttributeParameters& params) {
  if (!witSHCoefficients) return false;

  const double tr3x3[3][3][3] = {
    { {+0.2126, 0.7152, 0.0722}, { -0.1146, -0.3854, 0.5 }, { +0.5, -0.4542, -0.0458 } },
    { { 0.29900, 0.58700, 0.11400 }, { -0.16900, -0.33100, 0.50000 }, { 0.50000, -0.41900, -0.08100 } },
    { { 0.33333, 0.33333, 0.33333 }, { 0.40825, 0.00000, -0.40825 }, { 0.23570, -0.47140, 0.23570 } }
  };
  const double tr2x2[3][2][2] = {
    { { 1.0, 0.0 }, { 0.0, 1.0 } },
    { { 1.0, 0.0 }, { 0.0, 1.0 } },
    { { 0.5, 0.5 }, { 0.5, -0.5 } }
  };

  std::vector<int> colTr(16, 1);
  for (size_t gi = 0; gi < params.shTrGrpRow_.size(); ++gi) {
    if (params.shKltFlag_[gi] == 1) {
      for (int idx : params.shTrGrpRow_[gi]) { colTr[idx] = 0; }
    }
  }

  // Column transform
  for (size_t p = 0; p < pointCount; ++p) {
    for (size_t i = 0; i < 16; ++i) {
      if (colTr[i] == 1) {
        int type = params.shTrTypeCol_[i];
        double r = shView.sh(p, i, 0), g = shView.sh(p, i, 1), b = shView.sh(p, i, 2);
        shView.sh(p, i, 0) = tr3x3[type][0][0] * r + tr3x3[type][0][1] * g + tr3x3[type][0][2] * b;
        shView.sh(p, i, 1) = tr3x3[type][1][0] * r + tr3x3[type][1][1] * g + tr3x3[type][1][2] * b;
        shView.sh(p, i, 2) = tr3x3[type][2][0] * r + tr3x3[type][2][1] * g + tr3x3[type][2][2] * b;
      }
    }
  }

  // Row transform
  params.shKltParams_.clear();
  for (size_t gi = 0; gi < params.shTrGrpRow_.size(); ++gi) {
    const auto& groupidx = params.shTrGrpRow_[gi];
    int size = groupidx.size();
    if (params.shKltFlag_[gi] == 1) {
      bool isDcGroup = (groupidx.size() == 1 && groupidx[0] == 0);

      if (isDcGroup) {
        // DC PCA forward: 3x3 PCA on RGB of SH coeff 0, no sigma
        const int dcSize = 3;
        PCCKltParameters klt;
        klt.eigvec.assign(dcSize, std::vector<double>(dcSize, 0.0));
        klt.mean.assign(dcSize, 0.0);
        klt.sigma.assign(dcSize, 1.0);
        klt.eigvecFxp.assign(dcSize, std::vector<int64_t>(dcSize, 0));
        klt.meanFxp.assign(dcSize, 0);
        klt.sigmaFxp.assign(dcSize, 1 << KLT_SIGMA_SCALE);

        for (size_t p = 0; p < pointCount; ++p) {
          for (int ch = 0; ch < dcSize; ++ch)
            klt.mean[ch] += shView.sh(p, 0, ch);
        }
        for (int j = 0; j < dcSize; ++j) klt.mean[j] /= (double)pointCount;

        std::vector<double> cov(dcSize * dcSize, 0.0);
        std::vector<double> centeredDc(dcSize);
        for (size_t p = 0; p < pointCount; ++p) {
          for (int j = 0; j < dcSize; ++j) centeredDc[j] = shView.sh(p, 0, j) - klt.mean[j];
          for (int r = 0; r < dcSize; ++r)
            for (int c = 0; c < dcSize; ++c)
              cov[r + c * dcSize] += centeredDc[r] * centeredDc[c];
        }
        double div = double(pointCount > 1 ? pointCount - 1 : 1);
        for (int i = 0; i < dcSize * dcSize; ++i) cov[i] /= div;

        pcc::SymEigenResult eigResult;
        if (!pcc::symEigenDecomposeJacobi(cov, dcSize, eigResult)) {
          std::cerr << "DC PCA: symEigenDecomposeJacobi failed" << std::endl;
          continue;
        }
        std::vector<std::pair<double, int>> eigPairs;
        for (int i = 0; i < dcSize; ++i)
          eigPairs.push_back({eigResult.eigenvalues[i], i});
        std::sort(eigPairs.rbegin(), eigPairs.rend());
        for (int j = 0; j < dcSize; ++j) {
          int srcIdx = eigPairs[j].second;
          for (int i = 0; i < dcSize; ++i)
            klt.eigvec[j][i] = eigResult.eigenvectors[i + srcIdx * dcSize];
        }

        for (int j = 0; j < dcSize; ++j) {
          for (int i = 0; i < dcSize; ++i)
            klt.eigvecFxp[j][i] = (int64_t)round(klt.eigvec[j][i] * (1 << KLT_SCALE));
          klt.meanFxp[j] = (int64_t)round(klt.mean[j] * (1 << KLT_MEAN_SCALE));
          klt.sigmaFxp[j] = (int64_t)round(klt.sigma[j] * (1 << KLT_SIGMA_SCALE));
        }

        // Apply DC KLT: y = eigvec * (x - mean), no sigma
        std::vector<double> X(dcSize);
        for (size_t p = 0; p < pointCount; ++p) {
          for (int j = 0; j < dcSize; ++j) X[j] = (double)shView.sh(p, 0, j) - klt.mean[j];
          for (int j = 0; j < dcSize; ++j) {
            double y = 0.0;
            for (int i = 0; i < dcSize; ++i) y += klt.eigvec[j][i] * X[i];
            shView.sh(p, 0, j) = static_cast<float>(y);
          }
        }

        params.shKltParams_.push_back(klt);
        for (int idx : groupidx) { params.shTrTypeCol_[idx] = -1; }
      } else {
        PCCKltParameters klt;
        if (!computeShKltFromView(shView, pointCount, groupidx, klt)) continue;

        int fullSize = size * 3;
        std::vector<double> X(fullSize);
        for (size_t p = 0; p < pointCount; ++p) {
          for (int j = 0; j < fullSize; j++) {
            X[j] = (shView.sh(p, groupidx[j / 3], j % 3) - klt.mean[j]) / klt.sigma[j];
          }
          for (int j = 0; j < fullSize; j++) {
            double y = 0.0;
            for (int i = 0; i < fullSize; i++) { y += klt.eigvec[j][i] * X[i]; }
            shView.sh(p, groupidx[j / 3], j % 3) = static_cast<float>(y / 17.5);
          }
        }
        params.shKltParams_.push_back(klt);
        for (int idx : groupidx) { params.shTrTypeCol_[idx] = -1; }
      }
    }
    else {
      int type = params.shTrTypeRow_[gi];
      if (groupidx.size() == 3) {
        for (size_t p = 0; p < pointCount; ++p) {
          for (int c = 0; c < 3; ++c) {
            double v0 = shView.sh(p, groupidx[0], c);
            double v1 = shView.sh(p, groupidx[1], c);
            double v2 = shView.sh(p, groupidx[2], c);
            shView.sh(p, groupidx[0], c) = tr3x3[type][0][0] * v0 + tr3x3[type][0][1] * v1 + tr3x3[type][0][2] * v2;
            shView.sh(p, groupidx[1], c) = tr3x3[type][1][0] * v0 + tr3x3[type][1][1] * v1 + tr3x3[type][1][2] * v2;
            shView.sh(p, groupidx[2], c) = tr3x3[type][2][0] * v0 + tr3x3[type][2][1] * v1 + tr3x3[type][2][2] * v2;
          }
        }
      }
      else if (groupidx.size() == 2) {
        for (size_t p = 0; p < pointCount; ++p) {
          for (int c = 0; c < 3; ++c) {
            double v0 = shView.sh(p, groupidx[0], c);
            double v1 = shView.sh(p, groupidx[1], c);
            shView.sh(p, groupidx[0], c) = tr2x2[type][0][0] * v0 + tr2x2[type][0][1] * v1;
            shView.sh(p, groupidx[1], c) = tr2x2[type][1][0] * v0 + tr2x2[type][1][1] * v1;
          }
        }
      }
    }
  }
  return true;
}

bool sphericalHarmonicsTransform(bool witSHCoefficients, SHBufferView& shView,
  size_t pointCount, AttributeParameters& params, bool convertColourspace) {
  if (!witSHCoefficients) return false;
  if (params.gsdTransform2dPresentFlag_) {
    bool anyKlt = false;
    for (int flag : params.shKltFlag_) {
      if (flag == 1) { anyKlt = true; break; }
    }
    if (anyKlt) {
      const size_t shCoefCount = 16;
      if (shCoefCount > 1 && pointCount > 0) {
        double totalStd = 0.0;
        int count = 0;
        for (size_t l = 1; l < shCoefCount; ++l) {
          for (int ch = 0; ch < 3; ++ch) {
            double sum = 0.0, sumSq = 0.0;
            for (size_t p = 0; p < pointCount; ++p) {
              double val = shView.sh(p, l, ch);
              sum += val;
              sumSq += val * val;
            }
            double mean = sum / pointCount;
            double var = (sumSq / pointCount) - (mean * mean);
            totalStd += std::sqrt(std::max(0.0, var));
            count++;
          }
        }
        double meanStd = totalStd / count;
        if (meanStd < 0.05) {
          params.shTrTypeCol_.assign(16, 0);
          params.shTrTypeCol_[0] = 0;
          params.shTrGrpRow_ = { { 1, 5, 11 }, { 3, 7, 13 }, { 4, 10 }, { 8, 14 }, { 2, 6, 12 }, { 9, 15 } };
          params.shTrTypeRow_.assign(6, 2);
          params.shKltFlag_.assign(6, 0);
          params.dataAdaptiveShTransformEnabled_ = false;
        }
      }
    }
    return sphericalHarmonicsTransform2D(witSHCoefficients, shView, pointCount, params);
  }
  else {
    if (!convertColourspace) return true;
    
    for (size_t p = 0; p < pointCount; p++) {
      double sh[3];
      for (size_t m = 0; m < 16; m++) {
        sh[0] = shView.sh(p, m, 0);
        sh[1] = shView.sh(p, m, 1);
        sh[2] = shView.sh(p, m, 2);
        convertFromRgb(params.cicp_matrix_coefficients_idx, sh);
        shView.sh(p, m, 0) = sh[0];
        shView.sh(p, m, 1) = sh[1];
        shView.sh(p, m, 2) = sh[2];
      }
    }
  }
  return true;
}


//============================================================================

}  // namespace pcc
