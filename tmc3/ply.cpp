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
#include <string>
#include <unordered_map>

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

void RGB2YUV(double (&src)[3])
{
  double R = src[0];
  double G = src[1];
  double B = src[2];
  src[0] = +0.2126 * R + 0.7152 * G + 0.0722 * B;
  src[1] = -0.1146 * R - 0.3854 * G + 0.5 * B;
  src[2] = +0.5 * R - 0.4542 * G - 0.0458 * B;
  
}
void YUV2RGB(double (&src)[3])  
{
  double Y = src[0];
  double U = src[1];
  double V = src[2];
  src[0] = Y + 1.5748 * V;
  src[1] = Y - 0.1873 * U - 0.4681 * V;
  src[2] = Y + 1.8556 * U;
  
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
  // Apply the 6 scales permutation, and search for quaternion with maximum
  // absolute value component
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
  const std::vector<AttributeDescription>& attrDescs,
  bool asConformance)
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
    if (!asConformance) {
      fout << "property float " << attributeNames.position[0] << std::endl;
      fout << "property float " << attributeNames.position[1] << std::endl;
      fout << "property float " << attributeNames.position[2] << std::endl;
    }
    else {
      fout << "property float64 " << attributeNames.position[0] << std::endl;
      fout << "property float64 " << attributeNames.position[1] << std::endl;
      fout << "property float64 " << attributeNames.position[2] << std::endl;
    }
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
  if (cloud.hasSHCoefficients()) {
	  fout << "property float f_dc_0" << std::endl;
	  fout << "property float f_dc_1" << std::endl;
	  fout << "property float f_dc_2" << std::endl;

	  for (int i = 0; i < cloud.getNumSHCoeffs() - 1; i++) {
		  fout << "property float f_rest_" << 3 * i << std::endl;
		  fout << "property float f_rest_" << 3 * i + 1 << std::endl;
		  fout << "property float f_rest_" << 3 * i + 2 << std::endl;
	  }
  }
  if (cloud.hasOpacities()) {
	  fout << "property float opacity" << std::endl;
  }
  if (cloud.hasScales()) {
	  fout << "property float scale_0" << std::endl;
	  fout << "property float scale_1" << std::endl;
	  fout << "property float scale_2" << std::endl;
  }
  if (cloud.hasRotations()) {
	  fout << "property float rot_0" << std::endl;
	  fout << "property float rot_1" << std::endl;
	  fout << "property float rot_2" << std::endl;
	  fout << "property float rot_3" << std::endl;
  }
  if (asConformance) {
    fout << "element face 0" << std::endl;
    fout << "property list uint8 int32 vertex_index" << std::endl;
  }
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
			if (cloud.hasSHCoefficients()) {
				const Vec3<attr_t>& SHcoefficient_dc = cloud.getSHCoefficient(0, i);
				fout << " " << SHcoefficient_dc[0] 
					<< " "  << SHcoefficient_dc[1] 
					<< " "  << SHcoefficient_dc[2];

				for (int k = 0; k < 3; k++) {
					for (int m = 1; m < cloud.getNumSHCoeffs(); m++) {
						const Vec3<int32_t>& SHcoefficient_ac = cloud.getSHCoefficient(m, i);
						fout << " " << SHcoefficient_ac[k];
					}
				}
			}
			if (cloud.hasOpacities()) {
				fout << " " << cloud.getOpacity(i);
			}
			if (cloud.hasScales()) {
				const Vec3<attr_t>& scale = cloud.getScale(i);
				fout << " " << scale[0]
					<< " " << scale[1]
					<< " " << scale[2];
			}
			if (cloud.hasRotations()) {
				const pcc::Vec3<attr_t>& rotation_xyz = cloud.getRotationVec(i);
        fout << " " << cloud.getRotationScalar(i)
					<< " " << rotation_xyz[0]
					<< " " << rotation_xyz[1]
					<< " " << rotation_xyz[2];
			}
      fout << std::endl;
    }
  } else {
    fout.clear();
    fout.close();
    fout.open(fileName, std::ofstream::binary | std::ofstream::app);
    for (size_t i = 0; i < pointCount; ++i) {
      Vec3<double> position = cloud[i] * positionScale + positionOffset;
      if (!asConformance) {
        Vec3<float> positionFloat;
        positionFloat[0] = static_cast<float>(position[0]);
        positionFloat[1] = static_cast<float>(position[1]);
        positionFloat[2] = static_cast<float>(position[2]);
        fout.write(
          reinterpret_cast<const char* const>(&positionFloat), sizeof(float) * 3);
      }
      else 
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

      //for GS descs
      int shcoeffidx = (int)cloud.hasSHCoefficients(); 
      int opacityidx = cloud.getNumSHCoeffs() - 1 + (int)cloud.hasOpacities();
      int scaleidx = opacityidx + (int)cloud.hasScales();
      int rotationidx = scaleidx + (int)cloud.hasRotations();

			if (cloud.hasSHCoefficients()) {
				// SH DC
        double shdc[3];
        auto desc = attrDescs[0];
        auto offset = desc.params.attr_offset;
        auto attr_frac_bits = desc.params.attr_frac_bits;
        auto attr_scale     = desc.params.attr_scale_minus1 + 1;
        for (int k = 0; k < 3; k++) {
          int32_t attr = int32_t(cloud.getSHCoefficient(0, i)[k]);
          double buf = double(attr) * attr_scale / (1<<attr_frac_bits) + offset;
          shdc[k]= buf;
        }
        YUV2RGB(shdc);
				Vec3<float> val32b_dc{ float(shdc[0])/GS_PROPERTY_SCALE, float(shdc[1])/GS_PROPERTY_SCALE, float(shdc[2])/GS_PROPERTY_SCALE };
				fout.write(reinterpret_cast<char*>(&val32b_dc), sizeof(float) * 3);

				// SH AC
        double shac[15][3];
        desc = attrDescs[shcoeffidx];
        offset = desc.params.attr_offset;
        attr_frac_bits = desc.params.attr_frac_bits;
        attr_scale     = desc.params.attr_scale_minus1 + 1;
				for (int k = 0; k < 3; k++) {
					for (int m = 1; m < cloud.getNumSHCoeffs(); m++) {
            int32_t attr = int32_t(cloud.getSHCoefficient(m, i)[k]);
            double buf = double(attr) * attr_scale / (1<<attr_frac_bits) + offset;
            shac[m-1][k] = buf;
					}
				}
        for (size_t m =0 ; m < cloud.getNumSHCoeffs()-1 ; m++)
          YUV2RGB(shac[m]);
        for (int k = 0; k < 3; k++) {
					for (int m = 1; m < cloud.getNumSHCoeffs(); m++) {
						float val32b_ac = float(shac[m-1][k])/GS_PROPERTY_SCALE;
						fout.write(reinterpret_cast<char*>(&val32b_ac), sizeof(float));
					}
				}
			}
			if (cloud.hasOpacities()) {
        auto desc = attrDescs[opacityidx];
        auto offset = desc.params.attr_offset;
        auto attr_frac_bits = desc.params.attr_frac_bits;
        auto attr_scale     = desc.params.attr_scale_minus1 + 1;
				const int32_t& opacity = cloud.getOpacity(i);
        double _opacity = double(opacity) * attr_scale / (1<<attr_frac_bits) + offset;
				float buf = (float)_opacity/GS_PROPERTY_SCALE;
				fout.write(reinterpret_cast<char*>(&buf), sizeof(float));
			}
			if (cloud.hasScales()) {
        float scale[3];
        auto desc = attrDescs[scaleidx];
        auto offset = desc.params.attr_offset;
        auto attr_frac_bits = desc.params.attr_frac_bits;
        auto attr_scale     = desc.params.attr_scale_minus1 + 1;
        for (int k = 0; k < 3; k++) {
          int32_t attr = int32_t(cloud.getScale(i)[k]);
          double buf = double(attr) * attr_scale / (1<<attr_frac_bits) + offset;
          scale[k]= buf;
        }
				Vec3<float> val32b{ float(scale[0])/GS_PROPERTY_SCALE, float(scale[1])/GS_PROPERTY_SCALE, float(scale[2])/GS_PROPERTY_SCALE };
				fout.write(reinterpret_cast<char*>(&val32b), sizeof(float) * 3);
			}
			if (cloud.hasRotations()) {
        double rot[4] = {1.0, 0.0, 0.0, 0.0};
        auto desc = attrDescs[rotationidx];
        auto offset = desc.params.attr_offset;
        auto attr_frac_bits = desc.params.attr_frac_bits;
        auto attr_scale     = desc.params.attr_scale_minus1 + 1;
        for (int k = 0; k < 3; k++) {
          int32_t attr = int32_t(cloud.getRotationVec(i)[k]);
          double buf = double(attr) * attr_scale / (1<<attr_frac_bits) + offset;
          rot[k+1]= buf;
        }

        // Estimate the real part of the quaternion
        rot[0] = std::sqrt(std::min(std::max(
          1.0 - rot[1] * rot[1] - rot[2] * rot[2] - rot[3] * rot[3],
          0.0), 1.0));

				const auto rotation_w = rot[0];
				float _buf = (float)rotation_w/GS_PROPERTY_SCALE;
				fout.write(reinterpret_cast<char*>(&_buf), sizeof(float));

				Vec3<float> val32b{ float(rot[1])/GS_PROPERTY_SCALE, float(rot[2])/GS_PROPERTY_SCALE, float(rot[3])/GS_PROPERTY_SCALE };
				fout.write(reinterpret_cast<char*>(&val32b), sizeof(float) * 3);
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
  std::vector<AttributeDescription>& attrDescs)
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

  struct GSPoint {
    double shdc[3] = {0.0, 0.0, 0.0};
    double shac[15][3] = {};  // zero-initialized
    float opacity = 0.0f;
    double rot[4] = {1.0, 0.0, 0.0, 0.0};
    double scale[3] = {1.0, 1.0, 1.0 };
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
  // Gaussian Splat Compression
  size_t indexOpacity = PCC_UNDEFINED_INDEX;
  size_t indexScaleX = PCC_UNDEFINED_INDEX;
  size_t indexScaleY = PCC_UNDEFINED_INDEX;
  size_t indexScaleZ = PCC_UNDEFINED_INDEX;
  size_t indexRotationW = PCC_UNDEFINED_INDEX;
  size_t indexRotationX = PCC_UNDEFINED_INDEX;
  size_t indexRotationY = PCC_UNDEFINED_INDEX;
  size_t indexRotationZ = PCC_UNDEFINED_INDEX;
  std::vector<std::vector<size_t>> indexSHCoefficients;
  int numSHCoeffs = cloud.getNumSHCoeffs();
  indexSHCoefficients.resize(numSHCoeffs);
  for (int i = 0; i < numSHCoeffs; i++)
	  indexSHCoefficients[i].resize(3, PCC_UNDEFINED_INDEX);

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
    }	
		// Gaussian Splat Compression : Opacity
		else if (attributeInfo.name == "opacity" && attributeInfo.byteCount <= 4) {	
			indexOpacity = a;
		}
		// Gaussian Splat Compression : Scale
		else if (attributeInfo.name == "scale_0" && attributeInfo.byteCount <= 4) {
			indexScaleX = a;
		}
		else if (attributeInfo.name == "scale_1" && attributeInfo.byteCount <= 4) {
			indexScaleY = a;
		}
		else if (attributeInfo.name == "scale_2" && attributeInfo.byteCount <= 4) {
			indexScaleZ = a;
		}
		// Gaussian Splat Compression : Rotation
		else if (attributeInfo.name == "rot_0" && attributeInfo.byteCount <= 4) {
			indexRotationW = a;
		}
		else if (attributeInfo.name == "rot_1" && attributeInfo.byteCount <= 4) {
			indexRotationX = a;
		}
		else if (attributeInfo.name == "rot_2" && attributeInfo.byteCount <= 4) {
			indexRotationY = a;
		}
		else if (attributeInfo.name == "rot_3" && attributeInfo.byteCount <= 4) {
			indexRotationZ = a;
		}
		// Gaussian Splat Compression : SH DC coefficients
		else if (attributeInfo.name == "f_dc_0" && attributeInfo.byteCount <= 4) {
			indexSHCoefficients[0][0] = a;
		}
		else if (attributeInfo.name == "f_dc_1" && attributeInfo.byteCount <= 4) {
			indexSHCoefficients[0][1] = a;
		}
		else if (attributeInfo.name == "f_dc_2" && attributeInfo.byteCount <= 4) {
			indexSHCoefficients[0][2] = a;
		}
		// Gaussian Splat Compression : SH AC coefficients
		else if (!attributeInfo.name.compare(0, 7, "f_rest_") && attributeInfo.byteCount <= 4 
			&& numSHCoeffs > 1) {
			std::string suffix = attributeInfo.name.substr(7);
			int f_rest_index = std::stoi(suffix);

			int SHCoeff_index = f_rest_index % (numSHCoeffs - 1) + 1;
      int yuv_index = f_rest_index / (numSHCoeffs - 1);

			indexSHCoefficients[SHCoeff_index][yuv_index] = a;
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
  bool withOpacity = indexOpacity != PCC_UNDEFINED_INDEX;
  bool withScale = indexScaleX != PCC_UNDEFINED_INDEX
	  && indexScaleY != PCC_UNDEFINED_INDEX
	  && indexScaleZ != PCC_UNDEFINED_INDEX;
  bool withRotation = indexRotationW != PCC_UNDEFINED_INDEX
	  && indexRotationX != PCC_UNDEFINED_INDEX
	  && indexRotationY != PCC_UNDEFINED_INDEX
	  && indexRotationZ != PCC_UNDEFINED_INDEX;
  bool witSHCoefficients = indexSHCoefficients[0][0] != PCC_UNDEFINED_INDEX
	  && indexSHCoefficients[0][1] != PCC_UNDEFINED_INDEX
	  && indexSHCoefficients[0][2] != PCC_UNDEFINED_INDEX;

  cloud.addRemoveAttributes(withColors, withReflectances, withOpacity, withScale, withRotation, witSHCoefficients);

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
	  // Input format with Ascii is not supported in the current version of Gaussian Splat compression. 
      ++pointCounter;
    }
  } else {
    double shdc[3];
    double shac[15][3];
    double rot[4];
    float opacity;
    float scale[3];

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

    
    std::unordered_map<std::string, Stats> globalStats;  // keys: "f_dc","f_rest","opacity
    globalStats.reserve(5);
    std::vector<GSPoint> gsBuf;
    
    globalStats["f_dc"] = Stats{};
    globalStats["f_rest"] = Stats{};
    globalStats["opacity"] = Stats{};
    globalStats["rot"] = Stats{};
    globalStats["scale"] = Stats{};

    gsBuf.resize(pointCount);

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
        }

				// Gaussian Splat Compression
				else if (a == indexOpacity && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float _opacity;
						ifs.read(reinterpret_cast<char*>(&_opacity), sizeof(float));
            opacity = _opacity*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t _opacity;
						ifs.read(reinterpret_cast<char*>(&_opacity), sizeof(uint16_t));
            opacity = _opacity*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexScaleX && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float scaleX;
						ifs.read(reinterpret_cast<char*>(&scaleX), sizeof(float));
            scale[0] = scaleX*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t scaleX;
						ifs.read(reinterpret_cast<char*>(&scaleX), sizeof(uint16_t));
            scale[0] = scaleX*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexScaleY && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float scaleY;
						ifs.read(reinterpret_cast<char*>(&scaleY), sizeof(float));
            scale[1] = scaleY*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t scaleY;
						ifs.read(reinterpret_cast<char*>(&scaleY), sizeof(uint16_t));
            scale[1] = scaleY*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexScaleZ && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float scaleZ;
						ifs.read(reinterpret_cast<char*>(&scaleZ), sizeof(float));
            scale[2] = scaleZ*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t scaleZ;
						ifs.read(reinterpret_cast<char*>(&scaleZ), sizeof(uint16_t));
            scale[2] = scaleZ*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexRotationW && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float rotationW;
						ifs.read(reinterpret_cast<char*>(&rotationW), sizeof(float));
            rot[0] = rotationW*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t rotationW;
						ifs.read(reinterpret_cast<char*>(&rotationW), sizeof(uint16_t));
            rot[0] = rotationW*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexRotationX && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float rotationX;
						ifs.read(reinterpret_cast<char*>(&rotationX), sizeof(float));
            rot[1] = rotationX*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t rotationX;
						ifs.read(reinterpret_cast<char*>(&rotationX), sizeof(uint16_t));
            rot[1] = rotationX*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexRotationY && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float rotationY;
						ifs.read(reinterpret_cast<char*>(&rotationY), sizeof(float));
            rot[2] = rotationY*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t rotationY;
						ifs.read(reinterpret_cast<char*>(&rotationY), sizeof(uint16_t));
            rot[2] = rotationY*GS_PROPERTY_SCALE;
					}
				}
				else if (a == indexRotationZ && attributeInfo.byteCount <= 4) {
					if (attributeInfo.byteCount == 4) {
						float rotationZ;
						ifs.read(reinterpret_cast<char*>(&rotationZ), sizeof(float));
            rot[3] = rotationZ*GS_PROPERTY_SCALE;
					}
					else if (attributeInfo.byteCount == 2) {
						uint16_t rotationZ;
						ifs.read(reinterpret_cast<char*>(&rotationZ), sizeof(uint16_t));
            rot[3] = rotationZ*GS_PROPERTY_SCALE;
					}
				}
				else {
					bool isSHCoefficient = false;
					for (int i = 0; i < cloud.getNumSHCoeffs(); i++) {
						for (int j = 0; j < 3; j++) {
							if (a == indexSHCoefficients[i][j] && attributeInfo.byteCount <= 4) {
								if (attributeInfo.byteCount == 4) {
									float sh_coefficient;
									ifs.read(reinterpret_cast<char*>(&sh_coefficient), sizeof(float));
                  if (i==0)
                    shdc[j] = sh_coefficient*GS_PROPERTY_SCALE;
                  else
                    shac[i-1][j] = sh_coefficient*GS_PROPERTY_SCALE;

								}
								else if (attributeInfo.byteCount == 2) {
									uint16_t sh_coefficient;
									ifs.read(reinterpret_cast<char*>(&sh_coefficient), sizeof(uint16_t));          
                  if (i==0)
                    shdc[j] = sh_coefficient*GS_PROPERTY_SCALE;
                  else
                    shac[i-1][j] = sh_coefficient*GS_PROPERTY_SCALE;
								}
								isSHCoefficient = true;
								break;
							}
						}
					}

					if (!isSHCoefficient) {
						char buffer[128];
						ifs.read(buffer, attributeInfo.byteCount);
					}
        }
      }
      // GS Preprocessing
      if (witSHCoefficients || withOpacity || withScale || withRotation)
      {
        RGB2YUV(shdc);
        for (size_t m =0 ; m < 15 ; m++)
          RGB2YUV(shac[m]);

        canonicalize(rot, scale);

        
        // ---- Cache postprocessed GS values for this point ----
        GSPoint& g = gsBuf[pointCounter];

        // shdc (3)
        for (int c = 0; c < 3; ++c)
            g.shdc[c] = shdc[c];

        // shac (15x3)
        for (int m = 0; m < 15; ++m)
            for (int c = 0; c < 3; ++c)
              g.shac[m][c] = shac[m][c];

        // opacity
        g.opacity = opacity;

        // rot and scale (4 each)
        for (int c = 0; c < 4; ++c) {
            g.rot[c] = rot[c];
        }

        for (int c = 0; c < 3; ++c) {
            g.scale[c] = scale[c];
        }

        // ---- Update global min/max (single min/max per group) ----
        auto& st_dc = globalStats["f_dc"];
        auto& st_rest = globalStats["f_rest"];
        auto& st_op = globalStats["opacity"];
        auto& st_rot = globalStats["rot"];
        auto& st_scale = globalStats["scale"];

        // f_dc: 3 components share one min/max
        for (int c = 0; c < 3; ++c)
            st_dc.update(g.shdc[c]);

        // f_rest: 15x3 coefficients share one min/max
        for (int m = 0; m < 15; ++m)
            for (int c = 0; c < 3; ++c)
              st_rest.update(g.shac[m][c]);

        // opacity: scalar
        st_op.update(static_cast<double>(g.opacity));

        // rot: 4 components share one min/max
        for (int c = 0; c < 4; ++c)
            st_rot.update(g.rot[c]);

        // scale: 3 components share one min/max
        for (int c = 0; c < 3; ++c)
            st_scale.update(g.scale[c]);      


    }
  } 

    // second pass
    for (size_t pointCounter = 0; pointCounter < pointCount; ++pointCounter) {
        const GSPoint& g = gsBuf[pointCounter];

        //for descs
        int shcoeffidx = (int)cloud.hasSHCoefficients(); 
        int opacityidx = cloud.getNumSHCoeffs() - 1 + (int)cloud.hasOpacities();
        int scaleidx = opacityidx + (int)cloud.hasScales();
        int rotationidx = scaleidx + (int)cloud.hasRotations();

        // --- DC coefficients ---
        {
          auto& desc = attrDescs[0];
          auto attr_frac_bits_scale = 1 << desc.params.attr_frac_bits;
          const auto& st = globalStats["f_dc"];
          double min_val = st.min, max_val = st.max;
          int rounded_min = static_cast<int>(std::lround(min_val));
          int rounded_max = static_cast<int>(std::lround(max_val));

          // If both are zero, override to [-1, 1]
          if (rounded_min == 0.0) rounded_min = -1.0;
          if (rounded_max == 0.0) rounded_max = 1.0;

          int offset = rounded_min;
          int attr_scale = rounded_max - rounded_min;
          if (attr_scale <= 0)
              attr_scale = 1;

          desc.params.attr_offset = offset;
          desc.params.attr_scale_minus1 = attr_scale - 1;

          for (int c = 0; c < 3; ++c) {
              double attr = g.shdc[c];
              cloud.getSHCoefficient(0, pointCounter)[c] = attr_t(PCCClip(
                round((attr - offset) * attr_frac_bits_scale / attr_scale),
                0.0f, (1 << desc.bitdepth) - 1));
          }          
        }

        // --- AC coefficients ---
        {
          auto& desc = attrDescs[shcoeffidx];
          auto attr_frac_bits_scale = 1 << desc.params.attr_frac_bits;   
          const auto& st = globalStats["f_rest"];
          double min_val = st.min, max_val = st.max;
          int rounded_min = static_cast<int>(std::lround(min_val));
          int rounded_max = static_cast<int>(std::lround(max_val));

          // If both are zero, override to [-1, 1]
          if (rounded_min == 0.0) rounded_min = -1.0;
          if (rounded_max == 0.0) rounded_max = 1.0;

          int offset = rounded_min;
          int attr_scale = rounded_max - rounded_min;
          if (attr_scale <= 0)
              attr_scale = 1;

          desc.params.attr_offset = offset;
          desc.params.attr_scale_minus1 = attr_scale - 1;    

          for (int i = 0; i < cloud.getNumSHCoeffs() - 1; i++) {
						for (int j = 0; j < 3; j++) {
									double attr = g.shac[i][j];
									cloud.getSHCoefficient(i+1, pointCounter)[j] = attr_t(PCCClip(
                    round((attr - offset) * attr_frac_bits_scale / attr_scale),
                    0.0f, (1 << desc.bitdepth) - 1));
								}
          }
        }

        // --- Opacity ---
        {
          auto& desc = attrDescs[opacityidx];
          auto attr_frac_bits_scale = 1 << desc.params.attr_frac_bits;
          const auto& st = globalStats["opacity"];
          double min_val = st.min, max_val = st.max;
          int rounded_min = static_cast<int>(std::lround(min_val));
          int rounded_max = static_cast<int>(std::lround(max_val));

          if (rounded_min == 0.0) rounded_min = -1.0;
          if (rounded_max == 0.0) rounded_max = 1.0;

          int offset = rounded_min;
          int attr_scale = rounded_max - rounded_min;
          if (attr_scale <= 0)
              attr_scale = 1;

          desc.params.attr_offset = offset;
          desc.params.attr_scale_minus1 = attr_scale - 1;

          double attr = g.opacity;
          cloud.getOpacity(pointCounter) = attr_t(PCCClip(
              round((attr - offset) * attr_frac_bits_scale / attr_scale),
              0.0f, (1 << desc.bitdepth) - 1));
        }

        // --- Scale ---
        {
          auto& desc = attrDescs[scaleidx];
          auto attr_frac_bits_scale = 1 << desc.params.attr_frac_bits;
          const auto& st = globalStats["scale"];
          double min_val = st.min, max_val = st.max;
          int rounded_min = static_cast<int>(std::lround(min_val));
          int rounded_max = static_cast<int>(std::lround(max_val));

          if (rounded_min == 0.0) rounded_min = -1.0;
          if (rounded_max == 0.0) rounded_max = 1.0;

          int offset = rounded_min;
          int attr_scale = rounded_max - rounded_min;
          if (attr_scale <= 0)
              attr_scale = 1;

          desc.params.attr_offset = offset;
          desc.params.attr_scale_minus1 = attr_scale - 1;

          for (int c = 0; c < 3; ++c) {
            double attr = g.scale[c];
            cloud.getScale(pointCounter)[c] = attr_t(PCCClip(
                round((attr - offset) * attr_frac_bits_scale / attr_scale),
                0.0f, (1 << desc.bitdepth) - 1));
          }
        }

        // --- Rotation ---
        {
          auto& desc = attrDescs[rotationidx];
          auto attr_frac_bits_scale = 1 << desc.params.attr_frac_bits;
          const auto& st = globalStats["rot"];
          double min_val = st.min, max_val = st.max;
          int rounded_min = static_cast<int>(std::lround(min_val));
          int rounded_max = static_cast<int>(std::lround(max_val));

          if (rounded_min == 0.0) rounded_min = -1.0;
          if (rounded_max == 0.0) rounded_max = 1.0;

          int offset = rounded_min;
          int attr_scale = rounded_max - rounded_min;
          if (attr_scale <= 0)
              attr_scale = 1;

          desc.params.attr_offset = offset;
          desc.params.attr_scale_minus1 = attr_scale - 1;

          cloud.getRotationScalar(pointCounter) = attr_t(PCCClip(
              round((g.rot[0] - offset) * attr_frac_bits_scale / attr_scale),
              0.0f, (1 << desc.bitdepth) - 1));

          for (int c = 0; c < 3; ++c) {
            double attr = g.rot[c + 1];
            cloud.getRotationVec(pointCounter)[c] = attr_t(PCCClip(
                round((attr - offset) * attr_frac_bits_scale / attr_scale),
                0.0f, (1 << desc.bitdepth) - 1));
          }
        }
    }
  }

  return true;
}

//============================================================================

}  // namespace pcc
