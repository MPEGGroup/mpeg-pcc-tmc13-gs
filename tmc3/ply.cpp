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

#include <fstream>
#include <string>
#include <vector>

namespace pcc {

//============================================================================

static bool
compareSeparators(char aChar, const char* const sep)
{
  int i = 0;
  while (sep[i] != '\0') {
    if (aChar == sep[i])
      return false;
    i++;
  }
  return true;
}

//============================================================================

static bool
getTokens(
  const char* str, const char* const sep, std::vector<std::string>& tokens)
{
  if (!tokens.empty())
    tokens.clear();
  std::string buf = "";
  size_t i = 0;
  size_t length = ::strlen(str);
  while (i < length) {
    if (compareSeparators(str[i], sep)) {
      buf += str[i];
    } else if (buf.length() > 0) {
      tokens.push_back(buf);
      buf = "";
    }
    i++;
  }
  if (!buf.empty())
    tokens.push_back(buf);
  return !tokens.empty();
}

//============================================================================

void RGB2YUV(double (&src)[3])
{
  double R = src[0];
  double G = src[1];
  double B = src[2];
  src[0] = +0.29900*R + 0.58700*G + 0.11400*B;
  src[1] = -0.14713*R - 0.28886*G + 0.43600*B;
  src[2] = +0.61500*R - 0.51498*G - 0.10001*B;
}
void YUV2RGB(double (&src)[3])
{
  double Y = src[0];
  double U = src[1];
  double V = src[2];
  src[0] = Y             + 1.13983*V;
  src[1] = Y - 0.39465*U - 0.58060*V;
  src[2] = Y + 2.03211*U;
}

static void
permut_02(
  const double quat_in[4],
  const float scale_in[3],
  double quat_out[4],
  float scale_out[3])
{
  scale_out[0] = scale_in[2];
  scale_out[1] = scale_in[1];
  scale_out[2] = scale_in[0];
  quat_out[0] = quat_in[0] - quat_in[2];
  quat_out[1] = quat_in[1] - quat_in[3];
  quat_out[2] = quat_in[2] + quat_in[0];
  quat_out[3] = quat_in[3] + quat_in[1];
}

static void
permut_12(
  const double quat_in[4],
  const float scale_in[3],
  double quat_out[4],
  float scale_out[3])
{
  scale_out[0] = scale_in[0];
  scale_out[1] = scale_in[2];
  scale_out[2] = scale_in[1];
  quat_out[0] = quat_in[0] - quat_in[1];
  quat_out[1] = quat_in[1] + quat_in[0];
  quat_out[2] = quat_in[2] + quat_in[3];
  quat_out[3] = quat_in[3] - quat_in[2];
}

static void
canonicalize(double quat[4], float scale[3])
{
  // Test 6 scales permutation, and quaternion with maximum absolute value component
  double quats[6][4];
  float scales[6][3];

  std::copy_n(quat, 4, quats[0]);
  std::copy_n(scale, 3, scales[0]);
  permut_12(quats[0], scales[0], quats[1], scales[1]);
  permut_02(quats[1], scales[1], quats[2], scales[2]);
  permut_12(quats[2], scales[2], quats[3], scales[3]);
  permut_02(quats[3], scales[3], quats[4], scales[4]);
  permut_12(quats[4], scales[4], quats[5], scales[5]);

  double maxval = -1;
  int imax = 0;
  constexpr double mul[6] = {   // These values is to take into account that when
    1.00, 0.7071067811865476,   // the functions permut_02 or permut_12 are used,
    0.50, 0.3535533905932738,   // the resulting quaternion will have a norm equal
    0.25, 0.1767766952966369 }; // to the input quaternion norm times sqrt(2)
  for (int i = 0; i < 6; ++i) {
    double max = mul[i] * std::max(
      std::max(std::abs(quats[i][0]), std::abs(quats[i][1])),
      std::max(std::abs(quats[i][2]), std::abs(quats[i][3])));
    if (maxval < max) {
      maxval = max;
      imax = i;
    }
  }
  std::copy_n(scales[imax], 3, scale);
  if (maxval < 1e-6) {
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
    return;
  }
  
  // Put the maximum absolute value component of quaternion as the real part
  double w = std::abs(quats[imax][0]);
  double x = std::abs(quats[imax][1]);
  double y = std::abs(quats[imax][2]);
  double z = std::abs(quats[imax][3]);
  if (w >= x && w >= y && w >= z) {
    // copy
    quat[0] = +quats[imax][0];
    quat[1] = +quats[imax][1];
    quat[2] = +quats[imax][2];
    quat[3] = +quats[imax][3];
  } else if (x >= y && x >= z) {
    // rotate -pi radians about x
    quat[0] = +quats[imax][1];
    quat[1] = -quats[imax][0];
    quat[2] = -quats[imax][3];
    quat[3] = +quats[imax][2];
  } else if (y >= z) {
    // rotate -pi radians about y
    quat[0] = +quats[imax][2];
    quat[1] = +quats[imax][3];
    quat[2] = -quats[imax][0];
    quat[3] = -quats[imax][1];
  } else {
    // rotate -pi radians about z
    quat[0] = +quats[imax][3];
    quat[1] = -quats[imax][2];
    quat[2] = +quats[imax][1];
    quat[3] = -quats[imax][0];
  }

  // Scale quaternion
  double inorm = 1.0
    / std::sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2]
                + quat[3] * quat[3]);
  if (quat[0] < 0)
    inorm = -inorm;
  // Note: The real part of the quaternion is being ignored here, and will be
  // estimated at decoder side.
  quat[0] = 1.0;
  quat[1] *= inorm;
  quat[2] *= inorm;
  quat[3] *= inorm;
}

bool
ply::write(
  const PCCPointSet3& cloud,
  const PropertyNameMap& attributeNames,
  double positionScale,
  Vec3<double> positionOffset,
  const std::string& fileName,
  bool asAscii,
  const std::vector<AttributeDescription>& attrDescs)
{
  std::ofstream fout(fileName, std::ofstream::out);
  if (!fout.is_open()) {
    return false;
  }

  const size_t pointCount = cloud.getPointCount();
  fout << "ply" << std::endl;

  if (asAscii) {
    fout << "format ascii 1.0" << std::endl;
  } else {
    PCCEndianness endianess = PCCSystemEndianness();
    if (endianess == PCC_BIG_ENDIAN) {
      fout << "format binary_big_endian 1.0" << std::endl;
    } else {
      fout << "format binary_little_endian 1.0" << std::endl;
    }
  }
  fout << "element vertex " << pointCount << std::endl;
  if (asAscii) {
    fout << "property float " << attributeNames.position[0] << std::endl;
    fout << "property float " << attributeNames.position[1] << std::endl;
    fout << "property float " << attributeNames.position[2] << std::endl;
  } else {
    fout << "property float64 " << attributeNames.position[0] << std::endl;
    fout << "property float64 " << attributeNames.position[1] << std::endl;
    fout << "property float64 " << attributeNames.position[2] << std::endl;
  }

  if (cloud.hasColors()) {
    fout << "property uchar green" << std::endl;
    fout << "property uchar blue" << std::endl;
    fout << "property uchar red" << std::endl;
  }
  if (cloud.hasReflectances()) {
    fout << "property uint16 refc" << std::endl;
  }
  if (cloud.hasFrameIndex()) {
    fout << "property uint8 frameindex" << std::endl;
  }
  if (cloud.hasGSInfo()){
    auto gskeys = GSHelper::getGSKeys();
    for (auto key : gskeys){
      fout << "property float " << key << std::endl;
    }
  }
  fout << "element face 0" << std::endl;
  fout << "property list uint8 int32 vertex_index" << std::endl;
  fout << "end_header" << std::endl;
  if (asAscii) {
    //      fout << std::setprecision(std::numeric_limits<double>::max_digits10);
    fout << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < pointCount; ++i) {
      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      fout << position.x() << " " << position.y() << " " << position.z();
      if (cloud.hasColors()) {
        const Vec3<attr_t>& color = cloud.getColor(i);
        fout << " " << static_cast<int>(color[0]) << " "
             << static_cast<int>(color[1]) << " "
             << static_cast<int>(color[2]);
      }
      if (cloud.hasReflectances()) {
        fout << " " << static_cast<int>(cloud.getReflectance(i));
      }
      if (cloud.hasFrameIndex()) {
        fout << " " << static_cast<int>(cloud.getFrameIndex(i));
      }
      if (cloud.hasGSInfo()){
        for(int n=0; n < GSHelper::getGSKeysSize(); ++n){
          auto key = GSHelper::getGSKey(n);
          gs_t attr = gs_t(cloud.getGSInfo(i,key));
          float buf = float(attr) / GS_PROPERTY_SCALE;
          fout << " " << buf;
        }
      }
      fout << std::endl;
    }
  } else {
    fout.clear();
    fout.close();
    fout.open(fileName, std::ofstream::binary | std::ofstream::app);
    for (size_t i = 0; i < pointCount; ++i) {
      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      fout.write(
        reinterpret_cast<const char* const>(&position), sizeof(double) * 3);
      if (cloud.hasColors()) {
        const Vec3<attr_t>& c = cloud.getColor(i);
        Vec3<uint8_t> val8b{uint8_t(c[0]), uint8_t(c[1]), uint8_t(c[2])};
        fout.write(reinterpret_cast<const char*>(&val8b), sizeof(uint8_t) * 3);
      }
      if (cloud.hasReflectances()) {
        const attr_t& reflectance = cloud.getReflectance(i);
        fout.write(
          reinterpret_cast<const char*>(&reflectance), sizeof(uint16_t));
      }
      if (cloud.hasFrameIndex()) {
        const uint16_t& findex = cloud.getFrameIndex(i);
        fout.write(reinterpret_cast<const char*>(&findex), sizeof(uint16_t));
      }
      if (cloud.hasGSInfo()){
        double shdc[3];
        double shac[15][3];
        double rot[4] = {1.0, 0.0, 0.0, 0.0};
        float opacity;
        float scale[3];

        for(int n=0; n < GSHelper::getGSKeysSize(); ++n){
          auto key = GSHelper::getGSKey(n);
          int attr_instance_id = GSHelper::getGSInstanceId(key);
          auto desc = attrDescs[attr_instance_id];
          auto offset = desc.params.attr_offset;
          auto attr_frac_bits = desc.params.attr_frac_bits;
          auto attr_scale     = desc.params.attr_scale_minus1 + 1;
          gs_t attr = gs_t(cloud.getGSInfo(i,key));
          double buf = double(attr) * attr_scale / (1<<attr_frac_bits) + offset;
          if (key.find("f_dc") != std::string::npos) {
            int c = std::stoi(&key[5]);
            shdc[c] = buf;
          } else if (key.find("f_rest") != std::string::npos) {
            int c = std::stoi(&key[7]) / 15;
            int m = std::stoi(&key[7]) % 15;
            shac[m][c] = buf;
          } else if (key.find("opacity") != std::string::npos) {
            opacity = buf;
          } else if (key.find("rot") != std::string::npos) {
            int c = std::stoi(&key[4]);
            rot[c] = buf;
          } else if (key.find("scale") != std::string::npos) {
            int c = std::stoi(&key[6]);
            scale[c] = buf;
          }
        }

        YUV2RGB(shdc);
        for (size_t m =0 ; m < 15 ; m++)
          YUV2RGB(shac[m]);

        // Estimate the real part of the quaternion
        rot[0] = std::sqrt(std::min(std::max(
          1.0 - rot[1] * rot[1] - rot[2] * rot[2] - rot[3] * rot[3],
          0.0), 1.0));
        
        for(int n=0; n < GSHelper::getGSKeysSize(); ++n){
          auto key = GSHelper::getGSKey(n);
          float buf;
          if (key.find("f_dc") != std::string::npos) {
            int c = std::stoi(&key[5]);
            buf = shdc[c];
          } else if (key.find("f_rest") != std::string::npos) {
            int c = std::stoi(&key[7]) / 15;
            int m = std::stoi(&key[7]) % 15;
            buf = shac[m][c];
          } else if (key.find("opacity") != std::string::npos) {
            buf = opacity;
          } else if (key.find("rot") != std::string::npos) {
            int c = std::stoi(&key[4]);
            buf = rot[c];
          } else if (key.find("scale") != std::string::npos) {
            int c = std::stoi(&key[6]);
            buf = scale[c];
          }
          buf /= GS_PROPERTY_SCALE;
          fout.write(reinterpret_cast<char*>(&buf), sizeof(float));
        }
      }
    }
  }
  fout.close();
  return true;
}

//============================================================================

bool
ply::read(
  const std::string& fileName,
  const PropertyNameMap& attributeNames,
  double positionScale,
  PCCPointSet3& cloud,
  const std::vector<AttributeDescription>& attrDescs)
{
  std::ifstream ifs(fileName, std::ifstream::in | std::ifstream::binary);
  if (!ifs.is_open()) {
    return false;
  }
  enum AttributeType
  {
    ATTRIBUTE_TYPE_FLOAT64 = 0,
    ATTRIBUTE_TYPE_FLOAT32 = 1,
    ATTRIBUTE_TYPE_UINT64 = 2,
    ATTRIBUTE_TYPE_UINT32 = 3,
    ATTRIBUTE_TYPE_UINT16 = 4,
    ATTRIBUTE_TYPE_UINT8 = 5,
    ATTRIBUTE_TYPE_INT64 = 6,
    ATTRIBUTE_TYPE_INT32 = 7,
    ATTRIBUTE_TYPE_INT16 = 8,
    ATTRIBUTE_TYPE_INT8 = 9,
  };
  struct AttributeInfo {
    std::string name;
    AttributeType type;
    size_t byteCount;
    size_t index;
  };

  std::vector<AttributeInfo> attributesInfo;
  attributesInfo.reserve(16);
  const size_t MAX_BUFFER_SIZE = 4096;
  char tmp[MAX_BUFFER_SIZE];
  const char* sep = " \t\r";
  std::vector<std::string> tokens;

  ifs.getline(tmp, MAX_BUFFER_SIZE);
  getTokens(tmp, sep, tokens);
  if (tokens.empty() || tokens[0] != "ply") {
    std::cout << "Error: corrupted file!" << std::endl;
    return false;
  }
  bool isAscii = false;
  double version = 1.0;
  size_t pointCount = 0;
  bool isVertexProperty = true;
  while (1) {
    if (ifs.eof()) {
      std::cout << "Error: corrupted header!" << std::endl;
      return false;
    }
    ifs.getline(tmp, MAX_BUFFER_SIZE);
    getTokens(tmp, sep, tokens);
    if (tokens.empty() || tokens[0] == "comment") {
      continue;
    }
    if (tokens[0] == "format") {
      if (tokens.size() != 3) {
        std::cout << "Error: corrupted format info!" << std::endl;
        return false;
      }
      isAscii = tokens[1] == "ascii";
      version = atof(tokens[2].c_str());
    } else if (tokens[0] == "element") {
      if (tokens.size() != 3) {
        std::cout << "Error: corrupted element info!" << std::endl;
        return false;
      }
      if (tokens[1] == "vertex") {
        pointCount = atoi(tokens[2].c_str());
      } else {
        isVertexProperty = false;
      }
    } else if (tokens[0] == "property" && isVertexProperty) {
      if (tokens.size() != 3) {
        std::cout << "Error: corrupted property info!" << std::endl;
        return false;
      }
      const std::string& propertyType = tokens[1];
      const std::string& propertyName = tokens[2];
      const size_t attributeIndex = attributesInfo.size();
      attributesInfo.resize(attributeIndex + 1);
      AttributeInfo& attributeInfo = attributesInfo[attributeIndex];
      attributeInfo.name = propertyName;
      if (propertyType == "float64") {
        attributeInfo.type = ATTRIBUTE_TYPE_FLOAT64;
        attributeInfo.byteCount = 8;
      } else if (propertyType == "float" || propertyType == "float32") {
        attributeInfo.type = ATTRIBUTE_TYPE_FLOAT32;
        attributeInfo.byteCount = 4;
      } else if (propertyType == "uint64") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT64;
        attributeInfo.byteCount = 8;
      } else if (propertyType == "uint32") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT32;
        attributeInfo.byteCount = 4;
      } else if (propertyType == "uint16") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT16;
        attributeInfo.byteCount = 2;
      } else if (propertyType == "uchar" || propertyType == "uint8") {
        attributeInfo.type = ATTRIBUTE_TYPE_UINT8;
        attributeInfo.byteCount = 1;
      } else if (propertyType == "int64") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT64;
        attributeInfo.byteCount = 8;
      } else if (propertyType == "int32") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT32;
        attributeInfo.byteCount = 4;
      } else if (propertyType == "int16") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT16;
        attributeInfo.byteCount = 2;
      } else if (propertyType == "char" || propertyType == "int8") {
        attributeInfo.type = ATTRIBUTE_TYPE_INT8;
        attributeInfo.byteCount = 1;
      }
    } else if (tokens[0] == "end_header") {
      break;
    }
  }
  if (version != 1.0) {
    std::cout << "Error: non-supported version!" << std::endl;
    return false;
  }

  size_t indexX = PCC_UNDEFINED_INDEX;
  size_t indexY = PCC_UNDEFINED_INDEX;
  size_t indexZ = PCC_UNDEFINED_INDEX;
  size_t indexR = PCC_UNDEFINED_INDEX;
  size_t indexG = PCC_UNDEFINED_INDEX;
  size_t indexB = PCC_UNDEFINED_INDEX;
  size_t indexReflectance = PCC_UNDEFINED_INDEX;
  size_t indexFrame = PCC_UNDEFINED_INDEX;
  size_t indexNX = PCC_UNDEFINED_INDEX;
  size_t indexNY = PCC_UNDEFINED_INDEX;
  size_t indexNZ = PCC_UNDEFINED_INDEX;
  std::unordered_map<size_t,std::string> bytepos2gskey;
  const size_t attributeCount = attributesInfo.size();
  for (size_t a = 0; a < attributeCount; ++a) {
    auto& attributeInfo = attributesInfo[a];
    if (
      attributeInfo.name == attributeNames.position[0]
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexX = a;
    } else if (
      attributeInfo.name == attributeNames.position[1]
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexY = a;
    } else if (
      attributeInfo.name == attributeNames.position[2]
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexZ = a;
    } else if (attributeInfo.name == "red" && attributeInfo.byteCount == 1) {
      indexR = a;
    } else if (attributeInfo.name == "green" && attributeInfo.byteCount == 1) {
      indexG = a;
    } else if (attributeInfo.name == "blue" && attributeInfo.byteCount == 1) {
      indexB = a;
    } else if (
      (attributeInfo.name == "reflectance" || attributeInfo.name == "refc")
      && attributeInfo.byteCount <= 2) {
      indexReflectance = a;
    } else if (
      attributeInfo.name == "frameindex" && attributeInfo.byteCount <= 2) {
      indexFrame = a;
    } else if (
      attributeInfo.name == "nx"
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexNX = a;
    } else if (
      attributeInfo.name == "ny"
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexNY = a;
    } else if (
      attributeInfo.name == "nz"
      && (attributeInfo.byteCount == 8 || attributeInfo.byteCount == 4)) {
      indexNZ = a;
    } else if (GSHelper::isExist(attributeInfo.name)){
      bytepos2gskey[a] = attributeInfo.name;
    }
  }
  if (
    indexX == PCC_UNDEFINED_INDEX || indexY == PCC_UNDEFINED_INDEX
    || indexZ == PCC_UNDEFINED_INDEX) {
    std::cout << "Error: missing coordinates!" << std::endl;
    return false;
  }
  bool withColors = indexR != PCC_UNDEFINED_INDEX
    && indexG != PCC_UNDEFINED_INDEX && indexB != PCC_UNDEFINED_INDEX;
  bool withReflectances = indexReflectance != PCC_UNDEFINED_INDEX;
  bool withFrameIndex = indexFrame != PCC_UNDEFINED_INDEX;
  bool withGSInfo = bytepos2gskey.size() > 0;

  cloud.addRemoveAttributes(withColors, withReflectances, withGSInfo);
  if (withFrameIndex)
    cloud.addFrameIndex();
  else
    cloud.removeFrameIndex();

  cloud.resize(pointCount);
  if (isAscii) {
    size_t pointCounter = 0;
    while (!ifs.eof() && pointCounter < pointCount) {
      ifs.getline(tmp, MAX_BUFFER_SIZE);
      getTokens(tmp, sep, tokens);
      if (tokens.empty()) {
        continue;
      }
      if (tokens.size() < attributeCount) {
        return false;
      }
      auto& position = cloud[pointCounter];
      position[0] = atof(tokens[indexX].c_str()) * positionScale;
      position[1] = atof(tokens[indexY].c_str()) * positionScale;
      position[2] = atof(tokens[indexZ].c_str()) * positionScale;
      if (cloud.hasColors()) {
        auto& color = cloud.getColor(pointCounter);
        color[0] = atoi(tokens[indexG].c_str());
        color[1] = atoi(tokens[indexB].c_str());
        color[2] = atoi(tokens[indexR].c_str());
      }
      if (cloud.hasReflectances()) {
        cloud.getReflectance(pointCounter) =
          uint16_t(atoi(tokens[indexReflectance].c_str()));
      }
      if (cloud.hasFrameIndex()) {
        cloud.getFrameIndex(pointCounter) =
          uint8_t(atoi(tokens[indexFrame].c_str()));
      }
      ++pointCounter;
    }
  } else {
    double shdc[3];
    double shac[15][3];
    double rot[4];
    float opacity;
    float scale[3];
    for (size_t pointCounter = 0; pointCounter < pointCount && !ifs.eof();
         ++pointCounter) {
      auto& position = cloud[pointCounter];
      for (size_t a = 0; a < attributeCount && !ifs.eof(); ++a) {
        const auto& attributeInfo = attributesInfo[a];
        if (a == indexX) {
          if (attributeInfo.byteCount == 4) {
            float x;
            ifs.read(reinterpret_cast<char*>(&x), sizeof(float));
            position[0] = x * positionScale;
          } else {
            double x;
            ifs.read(reinterpret_cast<char*>(&x), sizeof(double));
            position[0] = x * positionScale;
          }
        } else if (a == indexY) {
          if (attributeInfo.byteCount == 4) {
            float y;
            ifs.read(reinterpret_cast<char*>(&y), sizeof(float));
            position[1] = y * positionScale;
          } else {
            double y;
            ifs.read(reinterpret_cast<char*>(&y), sizeof(double));
            position[1] = y * positionScale;
          }
        } else if (a == indexZ) {
          if (attributeInfo.byteCount == 4) {
            float z;
            ifs.read(reinterpret_cast<char*>(&z), sizeof(float));
            position[2] = z * positionScale;
          } else {
            double z;
            ifs.read(reinterpret_cast<char*>(&z), sizeof(double));
            position[2] = z * positionScale;
          }
        } else if (a == indexR && attributeInfo.byteCount == 1) {
          uint8_t val8b;
          ifs.read(reinterpret_cast<char*>(&val8b), sizeof(uint8_t));
          cloud.getColor(pointCounter)[2] = val8b;
        } else if (a == indexG && attributeInfo.byteCount == 1) {
          uint8_t val8b;
          ifs.read(reinterpret_cast<char*>(&val8b), sizeof(uint8_t));
          cloud.getColor(pointCounter)[0] = val8b;
        } else if (a == indexB && attributeInfo.byteCount == 1) {
          uint8_t val8b;
          ifs.read(reinterpret_cast<char*>(&val8b), sizeof(uint8_t));
          cloud.getColor(pointCounter)[1] = val8b;
        } else if (a == indexReflectance && attributeInfo.byteCount <= 2) {
          if (attributeInfo.byteCount == 1) {
            uint8_t reflectance;
            ifs.read(reinterpret_cast<char*>(&reflectance), sizeof(uint8_t));
            cloud.getReflectance(pointCounter) = reflectance;
          } else {
            auto& reflectance = cloud.getReflectance(pointCounter);
            ifs.read(reinterpret_cast<char*>(&reflectance), sizeof(uint16_t));
          }
        } else if (a == indexFrame && attributeInfo.byteCount <= 2) {
          if (attributeInfo.byteCount == 1) {
            auto& findex = cloud.getFrameIndex(pointCounter);
            ifs.read(reinterpret_cast<char*>(&findex), sizeof(uint8_t));
          } else {
            uint16_t findex;
            ifs.read(reinterpret_cast<char*>(&findex), sizeof(uint16_t));
            cloud.getFrameIndex(pointCounter) = uint8_t(findex);
          }
        } else if ( bytepos2gskey.find(a) != bytepos2gskey.end()){
          auto key = bytepos2gskey[a];
          assert(attributeInfo.name == key);
          double attr;
          if(attributeInfo.byteCount == 4){
            float buf;
            ifs.read(reinterpret_cast<char*>(&buf), sizeof(float));
            attr = buf * GS_PROPERTY_SCALE;
          } else if (attributeInfo.byteCount == 2){
            uint16_t buf;
            ifs.read(reinterpret_cast<char*>(&buf), sizeof(uint16_t));
            attr = buf * GS_PROPERTY_SCALE;
          }
          if (key.find("f_dc") != std::string::npos) {
            int c = std::stoi(&key[5]);
            shdc[c] = attr;
          } else if (key.find("f_rest") != std::string::npos) {
            int c = std::stoi(&key[7]) / 15;
            int m = std::stoi(&key[7]) % 15;
            shac[m][c] = attr;
          } else if (key.find("opacity") != std::string::npos) {
            opacity = attr;
          } else if (key.find("rot") != std::string::npos) {
            int c = std::stoi(&key[4]);
            rot[c] = attr;
          } else if (key.find("scale") != std::string::npos) {
            int c = std::stoi(&key[6]);
            scale[c] = attr;
          }
        } else {
          char buffer[128];
          ifs.read(buffer, attributeInfo.byteCount);
        }
      }
      // GS Preprocessing
      if (withGSInfo)
      {
        RGB2YUV(shdc);
        for (size_t m =0 ; m < 15 ; m++)
          RGB2YUV(shac[m]);

        canonicalize(rot, scale);

        for (auto key : GSHelper::getGSKeys())
        {
          int attr_instance_id = GSHelper::getGSInstanceId(key);
          auto desc = attrDescs[attr_instance_id];
          auto offset = desc.params.attr_offset;
          auto attr_frac_bits_scale = 1 << desc.params.attr_frac_bits;
          auto attr_scale     = desc.params.attr_scale_minus1 + 1;
          double attr;
          
          if (key.find("f_dc") != std::string::npos) {
            int c = std::stoi(&key[5]);
            attr = shdc[c];
          } else if (key.find("f_rest") != std::string::npos) {
            int c = std::stoi(&key[7]) / 15;
            int m = std::stoi(&key[7]) % 15;
            attr = shac[m][c];
          } else if (key.find("opacity") != std::string::npos) {
            attr = opacity;
          } else if (key.find("rot") != std::string::npos) {
            int c = std::stoi(&key[4]);
            attr = rot[c];
          } else if (key.find("scale") != std::string::npos) {
            int c = std::stoi(&key[6]);
            attr = scale[c];
          }
          cloud.getGSInfo(pointCounter, key) = gs_t(PCCClip(round((attr - offset) * attr_frac_bits_scale / attr_scale), 0.0f, attr_frac_bits_scale-1));
        }
      }
    }
  }
  return true;
}

//============================================================================

}  // namespace pcc
