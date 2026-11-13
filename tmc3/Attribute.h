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

#pragma once

#include <memory>
#include <map>
#include "hls.h"
#include "PayloadBuffer.h"
#include "PCCPointSet.h"
#include "entropy.h"

namespace pcc {

//============================================================================

class AttributeContexts;

//============================================================================

class AttributeDecoderIntf {
public:
  virtual ~AttributeDecoderIntf();

  virtual void decode(
    const SequenceParameterSet& sps,
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    int geom_num_points_minus1,
    int minGeomNodeSizeLog2,
    int minGeomNodeSizeLog2ForSH,
    const char* payload,
    size_t payloadLen,
    AttributeContexts& ctxtMem,
    PCCPointSet3& pointCloud,
    int rahtMaxAttributeGroupNum,
    std::unordered_map<int, std::vector<int>>* bufferRahtCoefficients) = 0;

    virtual void decodeMultiAttributeRAHT(
    const SequenceParameterSet& sps,
    std::map<int, AttributeParameterSet> _apss,
    std::vector<AttributeBrickHeader>& _allAttrAbh,
    int geom_num_points_minus1,
    int minGeomNodeSizeLog2,
    std::unordered_map<int, std::vector<int>>& bufferRahtCoefficients,
    PCCPointSet3& pointCloud,
    int maxAttributeNumberPerGroup) = 0;

  // Indicates if the attribute decoder can decode the given aps
  virtual bool isReusable(
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh) const = 0;

  virtual void outputPartialPoints(PCCPointSet3& pointCloud, PCCPointSet3& accumPointCloud) =0;
};

//----------------------------------------------------------------------------

std::unique_ptr<AttributeDecoderIntf> makeAttributeDecoder();

//============================================================================

class AttributeEncoderIntf {
public:
  virtual ~AttributeEncoderIntf();

  virtual void encode(
    const SequenceParameterSet& sps,
    const AttributeDescription& desc,
    const AttributeParameterSet& attr_aps,
    AttributeBrickHeader& abh,
    AttributeContexts& ctxtMem,
    PCCPointSet3& pointCloud,
    PayloadBuffer* payload,
    std::unordered_map<int, std::vector<int>> *bufferRahtCoefficients) = 0;

   virtual void encodeMultiAttributeRAHT(
    const SequenceParameterSet& sps,
    std::vector<AttributeParameterSet*> attr_aps,
    PCCPointSet3& pointCloud,
    std::vector<std::pair<std::string, int>> attributeIdxMap,
    int maxAttributeNumberPerGroup,
    std::unordered_map<int, std::vector<int>>* bufferRahtCoefficients) = 0;

  // Indicates if the attribute decoder can decode the given aps
  virtual bool isReusable(
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh) const = 0;
};

//----------------------------------------------------------------------------

std::unique_ptr<AttributeEncoderIntf> makeAttributeEncoder();

//============================================================================

int estimateDist2(
  const PCCPointSet3& cloud,
  int samplingPeriod,
  int searchRange,
  float percentileEstimate);

//============================================================================

void applyShAcInterAttributeScaling(
  pcc::PCCPointSet3& pointCloud,
  const std::vector<AttributeDescription>& attributeSets,
  bool isForwardTransform);

// Validate the attribute order required by SH AC inter-attribute scaling:
// opacity and scale must precede the SH coefficients attribute. Returns false
// (after reporting an error) when a required attribute is missing or ordered
// incorrectly.
bool validateShAcScalingAttributeOrder(
  const std::vector<AttributeDescription>& attributeSets);

// Returns true when SH AC scaling is enabled for the given attribute
// (i.e. it is an SH coefficient attribute with the scaling flag set).
bool isShAcScalingEnabled(const AttributeDescription& attributeSet);

//============================================================================

}  // namespace pcc
