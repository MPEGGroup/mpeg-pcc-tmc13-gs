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

#include <array>
#include <string>

#include "PCCPointSet.h"
#include "hls.h"

namespace pcc {
namespace ply {

  //============================================================================
  // This defines the the attribute names to be used when reading/writing ply
  // files.  The order of the names define the ordering used in the internal
  // point cloud representation.

  struct PropertyNameMap {
    // The names of the position attributes, typically {"x", "y", "z"}
    std::array<const char*, 3> position;
  };

  //============================================================================

  ///
  // Write @a pointCloud to a PLY file called @a fileName.
  // Each point position, pt, is converted prior to writing by:
  //  pt' = pt * positionScale + offset
  //
  // @param pointCloud  points to be written.
  // @param propertyNames  Describes ply property names of pointcloud attributes.
  // @param positionScale  scale factor for positions.
  // @param positionOffset  offset for positions (after scaling).
  // @param fileName  output filename.
  // @param asAscii  PLY writing format (true => ascii, false => binary).
  bool write(
    const PCCPointSet3& pointCloud,
    const PropertyNameMap& propertyNames,
    double positionScale,
    Vec3<double> positionOffset,
    const std::string& fileName,
    bool asAscii,
    const std::vector<AttributeDescription>& attrDescs,
    bool asConformance,
    bool externalPrePost
  );

  ///
  // Read @a pointCloud to a PLY file called @a fileName.
  // Point positions are scaled by positionScale and converted to integers.
  //
  bool read(
    const std::string& fileName,
    const PropertyNameMap& propertyNames,
    double positionScale,
    PCCPointSet3& cloud,
    std::vector<AttributeDescription>& attrDescs,
    bool externalPrePost
  );
  
  enum class Format { ASCII, BINARY };
  
  enum class PlyType {
    INT8, UINT8, INT16, UINT16, INT32, UINT32, FLOAT32, FLOAT64, UNKNOWN
  };

  struct Stats {
    // Min/Max (same behavior as MinMax)
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();

    // Accumulators for mean/stddev
    uint64_t count = 0;
    double sum = 0.0;
    double sumsq = 0.0;  // sum of v^2

    void update(double v)
    {
      if (!std::isfinite(v)) return; 
      if (v < min)
        min = v;
      if (v > max)
        max = v;
      ++count;
      sum += v;
      sumsq += v * v;
    }

    bool empty() const { return count == 0; }
    double mean() const
    {
      return empty() ? 0.0 : (sum / static_cast<double>(count));
    }

    // Population stddev (use this if you�re describing all data)
    double stddev_population() const
    {
      if (count == 0)
        return 0.0;
      double m = mean();
      double variance = (sumsq / static_cast<double>(count)) - m * m;
      return variance > 0.0 ? std::sqrt(variance) : 0.0;
    }

    // Sample stddev (use this if your data is a sample of a larger population)
    double stddev_sample() const
    {
      if (count <= 1)
        return 0.0;
      double n = static_cast<double>(count);
      double variance = (sumsq - (sum * sum) / n) / (n - 1.0);
      return variance > 0.0 ? std::sqrt(variance) : 0.0;
    }
  };  

  struct PropertyInfo {
    std::string name;
    KnownAttributeLabel attrLabel;
    int channel;
    PlyType type;
  };

  void stringToAttributeLabelAndChannel(const std::string& propName, KnownAttributeLabel& kAttrLabel, int& channelIdx);

  bool parseHeader(std::ifstream& ifs, PCCPointSet3& cloud, Format& format, std::vector<PropertyInfo>& properties);

  bool isRequestedAttribute(PropertyInfo& prop, const std::vector<AttributeDescription>& attrDescs);

  bool readVertices(std::ifstream& ifs, PCCPointSet3& cloud, const std::vector<AttributeDescription>& attrDescs, double positionScale, Format& format, std::vector<PropertyInfo>& properties);

  bool convertColor (PCCPointSet3& cloud, std::vector<AttributeDescription>& attrDescs);

  bool canonicalizeRotationAndScale (PCCPointSet3& cloud);

  bool analyzeMinMax (PCCPointSet3& cloud, std::vector<AttributeDescription>& attrDescs);

  bool quantizationLinear (PCCPointSet3& cloud, std::vector<AttributeDescription>& attrDescs);

  inline PlyType stringToPlyType(const std::string& typeStr) {
      if (typeStr == "char" || typeStr == "int8") return PlyType::INT8;
      if (typeStr == "uchar" || typeStr == "uint8") return PlyType::UINT8;
      if (typeStr == "short" || typeStr == "int16") return PlyType::INT16;
      if (typeStr == "ushort" || typeStr == "uint16") return PlyType::UINT16;
      if (typeStr == "int" || typeStr == "int32") return PlyType::INT32;
      if (typeStr == "uint" || typeStr == "uint32") return PlyType::UINT32;
      if (typeStr == "float" || typeStr == "float32") return PlyType::FLOAT32;
      if (typeStr == "double" || typeStr == "float64") return PlyType::FLOAT64;
      return PlyType::UNKNOWN;
  }

  template<typename T>
  inline float readAsciiData(std::ifstream& ifs) {
    T val;
    ifs >> val;
    return static_cast<float>(val);
  }

  inline float readTypedDataAscii(std::ifstream& ifs, PlyType type) {
    switch (type) {
      case PlyType::INT8:    return readAsciiData<int8_t>(ifs);
      case PlyType::UINT8:   return readAsciiData<uint16_t>(ifs);
      case PlyType::INT16:   return readAsciiData<int16_t>(ifs);
      case PlyType::UINT16:  return readAsciiData<uint16_t>(ifs);
      case PlyType::INT32:   return readAsciiData<int32_t>(ifs);
      case PlyType::UINT32:  return readAsciiData<uint32_t>(ifs);
      case PlyType::FLOAT32: return readAsciiData<float>(ifs);
      case PlyType::FLOAT64: return readAsciiData<double>(ifs);
      default: return 0.0f;
    }
  }
  
  template<typename T>
  inline float readBinaryData(std::ifstream& ifs) {
    T val;
    ifs.read(reinterpret_cast<char*>(&val), sizeof(T));
    return static_cast<float>(val);
  }

  inline float readTypedDataBinary(std::ifstream& ifs, PlyType type) {
    switch (type) {
      case PlyType::INT8:    return readBinaryData<int8_t>(ifs);
      case PlyType::UINT8:   return readBinaryData<uint8_t>(ifs);
      case PlyType::INT16:   return readBinaryData<int16_t>(ifs);
      case PlyType::UINT16:  return readBinaryData<uint16_t>(ifs);
      case PlyType::INT32:   return readBinaryData<int32_t>(ifs);
      case PlyType::UINT32:  return readBinaryData<uint32_t>(ifs);
      case PlyType::FLOAT32: return readBinaryData<float>(ifs);
      case PlyType::FLOAT64: return readBinaryData<double>(ifs);
      default: return 0.0f;
    }
    return 0.0f;
  }

  //============================================================================

}  // namespace ply
}  // namespace pcc