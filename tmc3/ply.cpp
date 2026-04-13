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
#include <sstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <limits>

namespace pcc {

//============================================================================

void RGB2YUV(float (&src)[3])
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
  bool asConformance,
  bool externalPrePost=true)
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

        if (!externalPrePost)
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
        if (!externalPrePost)
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
        
        if (externalPrePost)
        rot[0] = double(int32_t(cloud.getRotationScalar(i))) * attr_scale / (1<<attr_frac_bits) + offset;
        
        for (int k = 0; k < 3; k++) {
          int32_t attr = int32_t(cloud.getRotationVec(i)[k]);
          double buf = double(attr) * attr_scale / (1<<attr_frac_bits) + offset;
          rot[k+1]= buf;
        }

        // Estimate the real part of the quaternion
        if (!externalPrePost)
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

void ply::stringToAttributeLabelAndChannel(const std::string& propName, KnownAttributeLabel& kAttrLabel, int& channelIdx) 
{
  if (propName == "x")
    channelIdx = 0;
  else if (propName == "y")
    channelIdx = 1;
  else if (propName == "z")
    channelIdx = 2;
  else if (propName == "red") {
    kAttrLabel = KnownAttributeLabel::kColour;
    channelIdx = 0;
  } 
  else if (propName == "green") {
    kAttrLabel = KnownAttributeLabel::kColour;
    channelIdx = 1;
  } 
  else if (propName == "blue") {
    kAttrLabel = KnownAttributeLabel::kColour;
    channelIdx = 2;
  } 
  else if (propName == "reflectance") {
    kAttrLabel = KnownAttributeLabel::kReflectance;
    channelIdx = 0;
  } 
  else if (propName == "opacity") {
    kAttrLabel = KnownAttributeLabel::kOpacity;
    channelIdx = 0;
  } 
  else if (propName == "frame_index") {
    kAttrLabel = KnownAttributeLabel::kFrameIndex;
    channelIdx = 0;
  } 
  else if (propName == "frame_number") {
    kAttrLabel = KnownAttributeLabel::kFrameNumber;
    channelIdx = 0;
  } 
  else if (propName == "material_id") {
    kAttrLabel = KnownAttributeLabel::kMaterialId;
    channelIdx = 0;
  } 
  else if (propName == "nx") {
    kAttrLabel = KnownAttributeLabel::kNormal;
    channelIdx = 0;
  } 
  else if (propName == "ny") {
    kAttrLabel = KnownAttributeLabel::kNormal;
    channelIdx = 1;
  } 
  else if (propName == "nz") {
    kAttrLabel = KnownAttributeLabel::kNormal;
    channelIdx = 2;
  } 
  else if (propName.find("scale_") == 0) {
    kAttrLabel = KnownAttributeLabel::kScale;
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1));
  }
  else if (propName.find("rot_") == 0) {
    kAttrLabel = KnownAttributeLabel::kRotation;
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1));
  }
  else if (propName.find("f_dc_") == 0) {
    kAttrLabel = KnownAttributeLabel::kSH_coefficients;
    channelIdx = 16 * std::stoi(propName.substr(propName.find_last_of('_') + 1));
  } else if (propName.find("f_rest_") == 0) {
    kAttrLabel = KnownAttributeLabel::kSH_coefficients;
    channelIdx = std::stoi(propName.substr(propName.find_last_of('_') + 1));
    channelIdx = 16 * (channelIdx / 15)  + (channelIdx % 15) + 1;
  } 
  else {
    kAttrLabel = KnownAttributeLabel::kOid;
    channelIdx = 0;
  }
}

bool ply::parseHeader(
  std::ifstream& ifs, PCCPointSet3& cloud, 
  Format& format, std::vector<PropertyInfo>& properties)
{
  std::string line, token;
  properties.clear();
  bool inVertexElement = false;
  int vertexCount = 0;
  while (std::getline(ifs, line)) {
    std::istringstream iss(line);
    iss >> token;

    if (token == "format") {
      std::string formatType;
      iss >> formatType;
      format = (formatType == "ascii") ? Format::ASCII : Format::BINARY;
    } 
    else if (token == "element") {
      std::string type;
      iss >> type;
      inVertexElement = (type == "vertex");
      if (inVertexElement) {
        iss >> vertexCount;
      }
    }
    else if (token == "property" && inVertexElement) {
      std::string typeStr, nameStr;
      iss >> typeStr >> nameStr;
      
      KnownAttributeLabel attrLabel;
      int channel;
      stringToAttributeLabelAndChannel(nameStr, attrLabel, channel);
      properties.push_back({nameStr, attrLabel, channel, stringToPlyType(typeStr)});
      if (nameStr != "x" && nameStr != "y" && nameStr != "z")
        cloud.setAttribute(attrLabel);
    }
    else if (token == "end_header") {
      cloud.resize(vertexCount);
      return vertexCount > 0 && !properties.empty();
    }
  }
  return false;
}

bool ply::isRequestedAttribute(PropertyInfo& prop, const std::vector<AttributeDescription>& attrDescs) {
  if (prop.name == "x" || prop.name == "y" || prop.name == "z")
    return true;
  
  for (auto &desc: attrDescs) {
    if (desc.attributeLabel.known_attribute_label == prop.attrLabel)
    {
      if (desc.attributeLabel.known_attribute_label == KnownAttributeLabel::kSH_coefficients)
      {
        if (int(desc.params.sh_coefficient_label) == (prop.channel % 16))
          return true;
      } else {
        return true;
      }
    }
  }
  return false;
}

bool ply::readVertices(
  std::ifstream& ifs, 
  PCCPointSet3& cloud, const std::vector<AttributeDescription>& attrDescs, double positionScale,
  Format& format, std::vector<PropertyInfo>& properties)
{
  std::vector<bool> keepProp(properties.size(), false);
  for (int i = 0; i < properties.size(); ++i) {
    PropertyInfo& prop = properties[i];
    if (isRequestedAttribute(prop, attrDescs)) {
      keepProp[i] = true;
    }
  }

  for (size_t i = 0; i < cloud.getPointCount(); ++i) {
    for (size_t j = 0; j < properties.size(); ++j) {
      PropertyInfo& prop = properties[j];
      float value;
      if (format == Format::ASCII) {
        value = readTypedDataAscii(ifs, prop.type);
      } else {
        value = readTypedDataBinary(ifs, prop.type);
      }
      if (keepProp[j]) {
        if (prop.name == "x" || prop.name == "y" || prop.name == "z") {
          cloud[i][prop.channel] = value * positionScale;
        } else {
          attr_t& attr = cloud.getAttributes(i, prop.attrLabel, prop.channel);
          reinterpret_cast<float*>(&attr)[0] = value;
        }
      }
    }
  }

  return true;
}


bool
ply::convertColor (
  PCCPointSet3& cloud,
  std::vector<AttributeDescription>& attrDescs) 
{
  int pointCount = cloud.getPointCount();

  for (auto& desc : attrDescs)
  {
    if (desc.attr_num_dimensions_minus1 == 2 
        && (desc.attributeLabel == KnownAttributeLabel::kColour 
        || desc.attributeLabel == KnownAttributeLabel::kSH_coefficients)) {

      auto& attrs = cloud.getAttrVec(desc);
      for (int i = 0; i < pointCount; i++)
      {
        float attr_flt[3];
        std::memcpy(attr_flt, &attrs[i][0], sizeof(float) * 3); 
        
        RGB2YUV(attr_flt);

        std::memcpy(&attrs[i][0], attr_flt, sizeof(float) * 3); 
      }
    }
  }

  return true;
}

bool
ply::canonicalizeRotationAndScale (
  PCCPointSet3& cloud) 
{
  if(!cloud.hasRotations() || !cloud.hasScales())
    return true;

  int pointCount = cloud.getPointCount();
  std::vector<Vec3<attr_t>>& scale_xyz = cloud.getScale();
  std::vector<Vec3<attr_t>>& rot_xyz = cloud.getRotationVec();
  std::vector<attr_t>& rot_w = cloud.getRotationScalar();

  double quat[4];
  float scale[3];
  for (int i = 0; i < pointCount; i++)
  {
    float rot[4];
    
    std::memcpy(&rot[0], &rot_w[i], sizeof(float));
    std::memcpy(&rot[1], &rot_xyz[i][0], sizeof(float) * 3);

    quat[0] = rot[0];
    quat[1] = rot[1];
    quat[2] = rot[2];
    quat[3] = rot[3];
    std::memcpy(&scale[0], &scale_xyz[i][0], sizeof(float) * 3);

    canonicalize(quat, scale);

    rot[0] = quat[0];
    rot[1] = quat[1];
    rot[2] = quat[2];
    rot[3] = quat[3];

    std::memcpy(&rot_w[i], &rot[0], sizeof(float));
    std::memcpy(&rot_xyz[i][0], &rot[1], sizeof(float) * 3);
    std::memcpy(&scale_xyz[i][0], &scale[0], sizeof(float) * 3);
  }

  return true;
}


bool
ply::analyzeMinMax (
  PCCPointSet3& cloud,
  std::vector<AttributeDescription>& attrDescs) 
{
  int pointCount = cloud.getPointCount();  
  std::vector<Stats> stats(attrDescs.size());

  for (int d = 0; d < attrDescs.size(); d++)
  {
    AttributeDescription& desc = attrDescs[d];
    Stats& stat = stats[d]; 
    
    if (desc.attr_num_dimensions_minus1 == 0) {
      auto& attrs = cloud.getAttrScalar(desc);
      for (int i = 0; i < pointCount; i++)
      {
        float attr_flt;
        std::memcpy(&attr_flt, &attrs[i], sizeof(float));
        stat.update(attr_flt);
      }
    } else {
      auto& attrs = cloud.getAttrVec(desc);
      for (int i = 0; i < pointCount; i++)
      {
        float attr_flt[3];
        std::memcpy(attr_flt, &attrs[i][0], sizeof(float) * 3); 
        for (int c = 0; c < 3; c++)
        {
          stat.update(attr_flt[c]);
        }
      }
    }
    
    int rounded_min = static_cast<int>(std::lround(stat.min));
    int rounded_max = static_cast<int>(std::lround(stat.max));

    if (rounded_min == 0) rounded_min = -1;
    if (rounded_max == 0) rounded_max = 1;

    int offset = rounded_min;
    int attr_scale = rounded_max - rounded_min;
    if (attr_scale <= 0)
    {
      attr_scale = 1;
      offset = 0;
    }

    desc.params.attr_offset = offset;
    desc.params.attr_scale_minus1 = attr_scale - 1;
  }
  return true;
}

bool
ply::quantizationLinear (
  PCCPointSet3& cloud,
  std::vector<AttributeDescription>& attrDescs) 
{
  int pointCount = cloud.getPointCount();  
  for (auto& desc : attrDescs)
  {
    int offset = desc.params.attr_offset;
    int attr_scale = desc.params.attr_scale_minus1 + 1;
    auto attr_frac_bits_scale = 1ULL << desc.params.attr_frac_bits;
    int clip_max = (1 << desc.bitdepth) - 1;

    if (desc.attr_num_dimensions_minus1 == 0) {
      auto& attrs = cloud.getAttrScalar(desc);
      for (int i = 0; i < pointCount; i++)
      {
        float attr_flt;
        std::memcpy(&attr_flt, &attrs[i], sizeof(float));

        attrs[i] = attr_t(PCCClip<int>(
          static_cast<int>(std::round((attr_flt - offset) * attr_frac_bits_scale / attr_scale)),
          0, clip_max));
      }
    } else {
      auto& attrs = cloud.getAttrVec(desc);
      for (int i = 0; i < pointCount; i++)
      {
        float attr_flt[3];
        std::memcpy(attr_flt, &attrs[i][0], sizeof(float) * 3); 

        for (int c = 0; c < 3; c++)
        {
          attrs[i][c] = attr_t(PCCClip<int>(
            static_cast<int>(std::round((attr_flt[c] - offset) * attr_frac_bits_scale / attr_scale)),
            0, clip_max));
        }
      }
    }
  }
  return true;
}

bool
ply::read(
  const std::string& fileName,
  const PropertyNameMap& attributeNames,
  double positionScale,
  PCCPointSet3& cloud,
  std::vector<AttributeDescription>& attrDescs, bool externalPrePost)
{
  std::ifstream ifs(fileName,std::ifstream::in | std::ios::binary);
  if (!ifs.is_open()) {
      std::cerr << "Error: can not open file (" << fileName << ")" << std::endl;
      return false;
  }

  Format format;
  std::vector<PropertyInfo> properties;

  if (!parseHeader(ifs, cloud, format, properties)) {
    return false;
  }

  readVertices(ifs, cloud, attrDescs, positionScale, format, properties);

  // Color conversion
  convertColor(cloud, attrDescs);

  // Rotation and Scale Canonicalization  
  canonicalizeRotationAndScale(cloud);

  // Analyze Min-Max
  analyzeMinMax(cloud, attrDescs);

  // Attributes Linear Quantization   
  quantizationLinear(cloud, attrDescs);

  return true;
}

//============================================================================

}  // namespace pcc
