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

#include "TMC3.h"
#include "colourspace.h"

namespace pcc {

//============================================================================
// Supports both double and fixed-point data paths
class SHBufferView {
public:
    // Layout: for each band l and channel ch, buffer[l*3 + ch] points to array of pointCount doubles
    SHBufferView(double** buffers, int numBands, size_t pointCount)
        : buffers_(buffers), numBands_(numBands), pointCount_(pointCount) {}

    // Construct with fixed-point buffers for full FXP support
    SHBufferView(double** buffers, int64_t** fxpBuffers, int* fxpScale, int numBands, size_t pointCount)
        : buffers_(buffers), fxpBuffers_(fxpBuffers), fxpScale_(fxpScale), numBands_(numBands), pointCount_(pointCount) {}

    // Access SH coefficient at (point, band, channel) - non-const version (double)
    inline double& sh(size_t pointIdx, int band, int channel) {
        return buffers_[band * 3 + channel][pointIdx];
    }

    // Access SH coefficient at (point, band, channel) - const version (double)
    inline const double& sh(size_t pointIdx, int band, int channel) const {
        return buffers_[band * 3 + channel][pointIdx];
    }

    // Access for fixed-point data
    inline int64_t& shFxp(size_t pointIdx, int band, int channel) {
        return fxpBuffers_[band * 3 + channel][pointIdx];
    }
    inline const int64_t& shFxp(size_t pointIdx, int band, int channel) const {
        return fxpBuffers_[band * 3 + channel][pointIdx];
    }

    // Access/modify fixed-point scale for a band
    inline int& fxpScale(int band) { return fxpScale_[band]; }
    inline int fxpScale(int band) const { return fxpScale_[band]; }

    // Set fixed-point buffer pointers
    void setFxpBuffers(int64_t** fxpBuffers, int* fxpScale) {
        fxpBuffers_ = fxpBuffers;
        fxpScale_ = fxpScale;
    }

    // Check if fixed-point buffers are available
    bool hasFxpBuffers() const { return fxpBuffers_ != nullptr && fxpScale_ != nullptr; }

    size_t pointCount() const { return pointCount_; }
    int numBands() const { return numBands_; }

    // Get raw buffer pointer for a specific band and channel
    double* getBuffer(int band, int channel) { return buffers_[band * 3 + channel]; }
    const double* getBuffer(int band, int channel) const { return buffers_[band * 3 + channel]; }
    int64_t* getFxpBuffer(int band, int channel) { return fxpBuffers_[band * 3 + channel]; }
    const int64_t* getFxpBuffer(int band, int channel) const { return fxpBuffers_[band * 3 + channel]; }

private:
    double** buffers_;      // Array of (numBands * 3) buffer pointers
    int64_t** fxpBuffers_ = nullptr;  // Optional fixed-point buffers
    int* fxpScale_ = nullptr;         // Optional fixed-point scales per band
    int numBands_;
    size_t pointCount_;
};

//============================================================================
// GSPoint: unified SH coefficient storage for 2D transform / NQ pipeline

struct GSPoint {
    double sh[16][3] = {};
    int64_t shFxp[16][3] = {};
    int8_t fxpScale[16] = {};
};

//============================================================================
// ShQuantRegion: per-region parameters for non-uniform quantization

struct ShQuantRegion {
    std::vector<double> region_start_real;
    std::vector<int> tp_int;
    double min_offset;
    std::vector<int64_t> region_start_real_fxp;
    int64_t min_offset_fxp;
};

//============================================================================
// SH 2D transform and NQ functions (implemented in transform.cpp / ply.cpp)

void initializeShTransformParams(AttributeParameters& params);

// Original GSPoint-based interface (kept for NQ pipeline compatibility)
bool sphericalHarmonicsTransform(bool withSHCoefficients, std::vector<GSPoint>& gsBuf,
    size_t pointCount, AttributeParameters& params, bool convertColourspace);

bool sphericalHarmonicsInvTransform(bool withSHCoefficients, std::vector<GSPoint>& gsBuf,
    size_t pointCount, AttributeParameters& params, bool convertColourspace);

// Zero-copy SHBufferView interface (Option 3 - avoids gsBuf construction overhead)
bool sphericalHarmonicsTransform(bool withSHCoefficients, SHBufferView& shView,
    size_t pointCount, AttributeParameters& params, bool convertColourspace);

bool sphericalHarmonicsInvTransform(bool withSHCoefficients, SHBufferView& shView,
    size_t pointCount, AttributeParameters& params, bool convertColourspace);

// Original GSPoint-based interface (deprecated, kept for backward compatibility)
bool extractSHQuantRegions(bool withSHCoefficients, const std::vector<GSPoint>& gsBuf,
    size_t pointCount, std::vector<AttributeDescription>& attrDescs,
    const double dc_min, const double dc_max, const double ac_min, const double ac_max);

// Zero-copy SHBufferView interface (preferred)
bool extractSHQuantRegions(bool withSHCoefficients, SHBufferView& shView,
    size_t pointCount, std::vector<AttributeDescription>& attrDescs);

// Scale power transform functions (implemented in powertransform.cpp)
void setScalePowerTransformParams(
    AttributeDescription& desc, const std::vector<double>& inputAttrBuf,
    const std::vector<int>& indices, int pointCount, int64_t dCode);

void ensureScalePowerTransformParams(AttributeDescription& desc);

int quantizeScalePowerTransform(
    const AttributeDescription& desc, double value, int c, int64_t clipMax);

void quantizeScalePowerTransform(
    const AttributeDescription& desc, const double* src, int c, int64_t clipMax,
    std::vector<Vec3<attr_t>>& dst, int pointCount);

void inverseScalePowerTransform(
    const AttributeDescription& desc, const Vec3<attr_t>* src,
    Vec3<float>* dst, int pointCount);

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
    Stopwatch* clock,
    const PCCPointSet3& pointCloud,
    const PropertyNameMap& propertyNames,
    double positionScale,
    Vec3<double> positionOffset,
    const std::string& fileName,
    bool asAscii,
    const std::vector<AttributeDescription>& attrDescs,
    bool asConformance,
    bool convertColourspace,
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
    bool convertColourspace,
    bool externalPrePost,
    Stopwatch* clock
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
    int maxChannel;
    PlyType type;
    int attr_idx;
  };

  bool stringToAttributeLabelAndChannel(const std::string& propName, KnownAttributeLabel& kAttrLabel, int& channelIdx, int& maxChannelIdx);

  bool parseHeader(std::ifstream& ifs, PCCPointSet3& cloud, Format& format, std::vector<PropertyInfo>& properties);

  bool isRequestedAttribute(PropertyInfo& prop, const std::vector<AttributeDescription>& attrDescs);

  bool readVertices(
    std::ifstream& ifs, 
    PCCPointSet3& cloud, const std::vector<AttributeDescription>& attrDescs, double positionScale,
    Format& format, std::vector<PropertyInfo>& properties, std::vector<double>& inputAttrBuf);

  bool
  convertColor (
    int pointCount,
    std::vector<double>& inputAttrBuf,
    std::vector<AttributeDescription>& attrDescs, std::vector<PropertyInfo>& properties);

  bool
  convertRotationAndScale (
    int pointCount,
    std::vector<double>& inputAttrBuf,
    std::vector<AttributeDescription>& attrDescs, 
    std::vector<PropertyInfo>& properties);

  bool
  analyzeMinMax (
    int pointCount,
    std::vector<double>& inputAttrBuf,
    std::vector<AttributeDescription>& attrDescs, 
    std::vector<PropertyInfo>& properties);

  bool
  quantizationLinear (
    PCCPointSet3& cloud,
    std::vector<double>& inputAttrBuf,
    std::vector<AttributeDescription>& attrDescs, 
    std::vector<PropertyInfo>& properties);

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