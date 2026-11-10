#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <TMC3.h>
#include <PCCTMC3Encoder.h>
#include <PCCTMC3Decoder.h>
#include <io_tlv.h>
#include <ply.h>

#include <pthread.h>

using namespace emscripten;
using namespace std;
using namespace pcc;

enum class OutputSystem
{
  // Output after global scaling, don't convert to external system
  kConformance = 0,

  // Scale output to external coordinate system
  kExternal = 1,
};

//----------------------------------------------------------------------------

struct Parameters {
  bool isDecoder;

  // command line parsing should adjust dist2 values according to PQS
  bool positionQuantizationScaleAdjustsDist2;

  // Scale factor to apply when loading the ply before integer conversion.
  // Eg, If source point positions are in fractional metres converting to
  // millimetres will allow some fidelity to be preserved.
  double inputScale;

  // Length of the output point clouds unit vectors.
  double outputUnitLength;

  // output mode for ply writing (binary or ascii)
  bool outputBinaryPly;

  // Fractional fixed-point bits retained in conformance output
  int outputFpBits;

  // Output coordinate system to use
  OutputSystem outputSystem;

  // Output header type to use
  bool outputHeaderType;

  // when true, configure the encoder as if no attributes are specified
  bool disableAttributeCoding;

  // Frame number of first file in input sequence.
  int firstFrameNum;

  // Number of frames to process.
  int frameCount;

  std::string uncompressedDataPath;
  std::string compressedStreamPath;
  std::string reconstructedDataPath;

  // Filename for saving recoloured point cloud (encoder).
  std::string postRecolorPath;

  // Filename for saving pre inverse scaled point cloud (decoder).
  std::string preInvScalePath;

  pcc::EncoderParams encoder;
  pcc::DecoderParams decoder;

  // perform attribute colourspace conversion on ply input/output.
  bool convertColourspace;

  // resort the input points by azimuth angle
  bool sortInputByAzimuth;

  // Maximum number of attribute groups that can be processed in one RAHT encoding or decoding pass.
  int maxAttributeNumberPerGroup;
};

//----------------------------------------------------------------------------

class SequenceCodec {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceCodec(Parameters* params) : params(params) {}

  // Perform conversions and write output point cloud
  //  \params cloud  a mutable copy of reconFrame.cloud
  void writeOutputFrame(
    Stopwatch* clock,
    const std::string& postInvScalePath,
    const std::string& preInvScalePath,
    const CloudFrame& reconFrame,
    PCCPointSet3& cloud);

  // determine the output ply scale factor
  double outputScale(const CloudFrame& cloud) const;

  // the output ply origin, scaled according to output coordinate system
  Vec3<double> outputOrigin(const CloudFrame& cloud) const;

  void scaleAttributesForInput(
    const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

  void scaleAttributesForOutput(
    const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

protected:
  Parameters* params;
};

//----------------------------------------------------------------------------

class SequenceEncoder
  : public SequenceCodec
  , PCCTMC3Encoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceEncoder(Parameters* params);

  int compress(Stopwatch* clock);

protected:
  int compressOneFrame(Stopwatch* clock);

  void onOutputBuffer(const PayloadBuffer& buf) override;
  void onPostRecolour(const PCCPointSet3& cloud) override;

private:
  ply::PropertyNameMap _plyAttrNames;

  // The raw origin used for input sorting
  Vec3<int> _angularOrigin;

  PCCTMC3Encoder3 encoder;

  std::ofstream bytestreamFile;

  int frameNum;
};

//----------------------------------------------------------------------------

class SequenceDecoder
  : public SequenceCodec
  , PCCTMC3Decoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceDecoder(Parameters* params);

  int decompress(Stopwatch* clock);

  PCCTMC3Decoder3& getDecoder() {return decoder;}
  PCCTMC3Decoder3::Callbacks* callback() {return this;}

protected:
  void onOutputCloud(const CloudFrame& cloud) override;

private:
  PCCTMC3Decoder3 decoder;

  std::ofstream bytestreamFile;

  Stopwatch* clock;
};

//============================================================================

//============================================================================
// SH remapping, as in SequenceDecoder::onOutputCloud() (TMC3.cpp)

inline void filterInvalidShRemapPairs(
  std::vector<std::pair<int, int>>& pairs, const int numShCoeffs)
{
  pairs.erase(std::remove_if(pairs.begin(), pairs.end(), [numShCoeffs](const std::pair<int, int>& p) {
      return (p.first & 15) >= numShCoeffs || (p.second & 15) >= numShCoeffs;
  }), pairs.end());
}

inline void applyGsShRemap(
  PCCPointSet3& pointCloud, const std::vector<std::pair<int, int>>& pairs)
{
  const size_t voxelCount = pointCloud.getPointCount();
  for (const auto& p : pairs) {
    const int c0 = p.first >> 4;
    const int m0 = p.first & 15;
    const int c1 = p.second >> 4;
    const int m1 = p.second & 15;
    auto& sh0 = pointCloud.getSHCoefficient(m0);
    auto& sh1 = pointCloud.getSHCoefficient(m1);
    for (size_t n = 0; n < voxelCount; ++n)
      std::swap(sh0[n][c0], sh1[n][c1]);
  }
}

//============================================================================
// The decoded frame is written by the PLY writer of tmc3 (ply::write) in the
// in-memory file system, then read back: the output is the same as the one of
// "tmc3 --mode=1 --convertPlyColourspace=1 --outputBinaryPly=1" (inverse SH
// transform, SH non-uniform dequantisation, colour space conversion, scale
// power transform, rotation formats), without duplicating that code here.

static const char* kOutputPly = "/tmc13_decoded.ply";

class TMC13Decoder {
public:
    TMC13Decoder() {}

    bool decode(emscripten::val data) {
        unsigned int length = data["length"].as<unsigned int>();
        if (length == 0) {
            std::cout << "Warning: Length is 0 inside C++" << std::endl;
            return false;
        }

        std::vector<uint8_t> bitstream(length);

        emscripten::val view = emscripten::val(typed_memory_view(length, bitstream.data()));
        view.call<void>("set", data);

        // same values as the defaults of the tmc3 command line (decoder)
        Parameters params{};
        params.isDecoder = 1;
        params.convertColourspace = 1;
        params.outputBinaryPly = 1;
        params.outputUnitLength = 1.0;  // srcUnitLength
        params.outputSystem = OutputSystem::kExternal;
        params.outputFpBits = -1;
        params.firstFrameNum = 0;
        params.maxAttributeNumberPerGroup = 1;
        params.encoder.externalPrePost = false;
        params.decoder.minGeomNodeSizeLog2 = 0;
        params.decoder.minGeomNodeSizeLog2forSH = 0;
        params.decoder.decodeMaxPoints = 0;
        params.decoder.outputFpBits = params.outputFpBits;
        params.decoder.maxAttributeNumberPerGroup = params.maxAttributeNumberPerGroup;
        params.decoder.targetSHDegree = -1;
        params.decoder.gs_ext = 1;
        SequenceDecoder seqDec = SequenceDecoder(&params);

        int readSize = 0;
        while (true) {
            pcc::PayloadBuffer buf;

            readTlv(bitstream, &buf, readSize);

            if (seqDec.getDecoder().decompress(&buf, seqDec.callback())) {
              std::cout << "Error: can't decompress point cloud!" << std::endl;
              return false;
            }

            if (readSize >= (int)length)
              break;
        }

        if (seqDec.getDecoder().decompress(nullptr, seqDec.callback())) {
          std::cout << "Error: can't decompress point cloud!" << std::endl;
          return false;
        }

        auto& outframe = seqDec.getDecoder().getOutputFrame();
        {
          // copy of the point cloud modified according to the output options
          PCCPointSet3 pointCloud(outframe.cloud);
          if (outframe.gs_sh_remap_enabled_flag && pointCloud.hasSHCoefficients()) {
            std::vector<std::pair<int, int>> gs_sh_remap_pairs =
              {{ 5, 17}, { 6, 18},  {7, 19}, { 9, 20}, {10, 36},
               {11, 33}, {12, 34}, {13, 35}, {14, 24}, {15, 40},
               {26, 41}, {27, 37}, {28, 38}, {29, 39}, {31, 46}};
            filterInvalidShRemapPairs(gs_sh_remap_pairs, pointCloud.getNumSHCoeffs());
            applyGsShRemap(pointCloud, gs_sh_remap_pairs);
          }
          outframe.cloud.clear();

          pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock;
          seqDec.writeOutputFrame(&clock, kOutputPly, {}, outframe, pointCloud);
        }

        const bool ok = readPly(kOutputPly);
        std::remove(kOutputPly);
        if (!ok) {
          std::cout << "Error: can't read the decoded point cloud!" << std::endl;
          return false;
        }

        std::cout << " G-PCC finish !!!" << std::endl;

        return true;
    }

    emscripten::val getAttrVals(size_t idx) {
        if (idx >= m_attrs.size())
            return emscripten::val::undefined();

        auto& v = m_attrs[idx];

        return emscripten::val(
            emscripten::typed_memory_view(v.size(), v.data())
        );
    }


    size_t getPointCount() const {
        return m_pointCount;
    }

private:
    // Reads the binary PLY written by ply::write and fills the 59 attribute arrays:
    // x y z, f_dc_0..2, (f_rest_m, f_rest_15+m, f_rest_30+m) for m = 0..14,
    // scale_0..2, rot_0..3 (w x y z), opacity. Missing attributes are set to 0.
    bool readPly(const char* fileName) {
        std::ifstream fin(fileName, std::ios::binary);
        if (!fin.is_open())
            return false;
        std::string line;
        size_t N = 0;
        std::vector<std::string> props;
        bool binaryLE = false;
        while (std::getline(fin, line)) {
            std::istringstream ls(line);
            std::string w0, w1, w2;
            ls >> w0 >> w1 >> w2;
            if (w0 == "format")
                binaryLE = (w1 == "binary_little_endian");
            else if (w0 == "element" && w1 == "vertex")
                N = std::stoul(w2);
            else if (w0 == "property") {
                if (w1 != "float")  // non-conformance output: float properties only
                    return false;
                props.push_back(w2);
            } else if (w0 == "end_header")
                break;
        }
        if (!binaryLE || props.empty())
            return false;

        auto index = [](const std::string& name) -> int {
            if (name == "x") return 0;
            if (name == "y") return 1;
            if (name == "z") return 2;
            if (name.rfind("f_dc_", 0) == 0) return 3 + std::stoi(name.substr(5));
            if (name.rfind("f_rest_", 0) == 0) {
                int k = std::stoi(name.substr(7));
                return 6 + (k % 15) * 3 + k / 15;
            }
            if (name.rfind("scale_", 0) == 0) return 51 + std::stoi(name.substr(6));
            if (name.rfind("rot_", 0) == 0) return 54 + std::stoi(name.substr(4));
            if (name == "opacity") return 58;
            return -1;
        };
        std::vector<int> dst(props.size());
        for (size_t p = 0; p < props.size(); p++)
            dst[p] = index(props[p]);

        m_pointCount = N;
        m_attrs.assign(59, std::vector<float>(N, 0.f));
        std::vector<float> row(props.size());
        for (size_t i = 0; i < N; ++i) {
            if (!fin.read(reinterpret_cast<char*>(row.data()), sizeof(float) * row.size()))
                return false;
            for (size_t p = 0; p < props.size(); p++)
                if (dst[p] >= 0 && dst[p] < 59)
                    m_attrs[dst[p]][i] = row[p];
        }
        return true;
    }

    std::vector<std::vector<float>> m_attrs;  // positoin, SHs, scales, rots, opacity
    size_t m_pointCount = 0;
};

EMSCRIPTEN_BINDINGS(my_module) {
    class_<TMC13Decoder>("TMC13Decoder")
        .constructor<>()
        .function("decode"       , &TMC13Decoder::decode)
        .function("getAttrVals"  , &TMC13Decoder::getAttrVals)
        .function("getPointCount", &TMC13Decoder::getPointCount);
}

extern "C" int main() {
  return 0;
}