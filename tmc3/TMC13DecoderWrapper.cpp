#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>

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
};

//----------------------------------------------------------------------------

class SequenceCodec {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceCodec(Parameters* params) : params(params) {}

  // Perform conversions and write output point cloud
  //  \params cloud  a mutable copy of reconFrame.cloud
  void writeOutputFrame(
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

inline void YUV2RGB(float& y, float& u, float& v)
{
    const float Y = y;
    const float U = u;
    const float V = v;

    y = Y + 1.5748f * V;                 // R
    u = Y - 0.1873f * U - 0.4681f * V;   // G
    v = Y + 1.8556f * U;                 // B
}

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

        Parameters params;
        params.isDecoder = 1;
        params.convertColourspace = 0;
        params.outputBinaryPly = 1;
        params.outputUnitLength = 1.0;
        params.outputSystem = OutputSystem::kExternal;
        params.decoder.minGeomNodeSizeLog2 = 0;
        params.decoder.gs_ext = 1;
        params.decoder.maxAttributeNumberPerGroup = 1;
        SequenceDecoder seqDec = SequenceDecoder(&params);
        std::vector<pcc::PayloadBuffer*> attr_bufs;

        int readSize = 0;
        bool isFirstADU=true;
        while (true) {
            pcc::PayloadBuffer* buf_ptr = new pcc::PayloadBuffer;

            readTlv(bitstream, buf_ptr, readSize);

            if (seqDec.getDecoder().decompress(buf_ptr, seqDec.callback())) {
              std::cout << "Error: can't decompress point cloud!" << std::endl;
              return false;
            }

            if (readSize == length)
              break;
        }

        if (seqDec.getDecoder().decompress(nullptr, seqDec.callback())) {
          std::cout << "Error: can't decompress point cloud!" << std::endl;
          return false;
        }

        auto& outframe = seqDec.getDecoder().getOutputFrame();
        PCCPointSet3& pc = outframe.cloud;

        const size_t N = m_pointCount = pc.getPointCount();

        // 3 + 48 + 4 + 3 + 1 = 59
        m_attrs.resize(59);
        for (int i = 0; i < 59; i++)
          m_attrs[i].resize(N);


        int attr_count = 0;
        // positions: m_attrs[0~2]
        {
            auto plyScale  = seqDec.outputScale(outframe) / (1 << outframe.outputFpBits);
            auto plyOrigin = seqDec.outputOrigin(outframe);

            const float scale  = plyScale;

            const float offset_x = plyOrigin[0];
            const float offset_y = plyOrigin[1];
            const float offset_z = plyOrigin[2];

            std::vector<float>& pos_x = m_attrs[attr_count++];
            std::vector<float>& pos_y = m_attrs[attr_count++];
            std::vector<float>& pos_z = m_attrs[attr_count++];



            for (size_t i = 0; i < N; ++i) {
                const auto& p = pc[i];  // pc[i][0..2]

                pos_x[i] = p[0] * scale + offset_x;
                pos_y[i] = p[1] * scale + offset_y;
                pos_z[i] = p[2] * scale + offset_z;
            }
        }

        {
          auto& desc = outframe.attrDesc[0];
          const float scale  = float(desc.params.attr_scale_minus1 + 1) / float(1 << desc.params.attr_frac_bits);
          const float offset = desc.params.attr_offset;
          
          // src (quantized)
          const std::vector<Vec3<attr_t>>& srcDc =  pc.getSHCoefficient(0);

          // dst (float)
          std::vector<float>& dstR = m_attrs[attr_count++];
          std::vector<float>& dstG = m_attrs[attr_count++];
          std::vector<float>& dstB = m_attrs[attr_count++];

          for (size_t i = 0; i < N; ++i) {
            float y = srcDc[i][0] * scale + offset;
            float u = srcDc[i][1] * scale + offset;
            float v = srcDc[i][2] * scale + offset;
            
            YUV2RGB(y, u, v);

            dstR[i] = y;
            dstG[i] = u;
            dstB[i] = v;
          }

          // m = SH / coefficient index
          for (int m = 1; m < 16; ++m) {
            auto& desc = outframe.attrDesc[m];
            const float scale  = float(desc.params.attr_scale_minus1 + 1) / float(1 << desc.params.attr_frac_bits);
            const float offset = desc.params.attr_offset;

            const std::vector<Vec3<attr_t>>& srcAc =  pc.getSHCoefficient(m);

            std::vector<float>& dstR = m_attrs[attr_count++];
            std::vector<float>& dstG = m_attrs[attr_count++];
            std::vector<float>& dstB = m_attrs[attr_count++];

            for (size_t i = 0; i < N; ++i) {
              float y = srcAc[i][0] * scale + offset;
              float u = srcAc[i][1] * scale + offset;
              float v = srcAc[i][2] * scale + offset;

              YUV2RGB(y, u, v);

              dstR[i] = y;
              dstG[i] = u;
              dstB[i] = v;
            }
          }
        }

        // scale: m_attrs[51, 52, 53]
        {
          auto& desc = outframe.attrDesc[17];
          const float scale  = float(desc.params.attr_scale_minus1 + 1) / float(1 << desc.params.attr_frac_bits);
          const float offset = desc.params.attr_offset;

            const std::vector<Vec3<attr_t>>& srcS =  pc.getScale();
            std::vector<float>& dstSx = m_attrs[attr_count++];
            std::vector<float>& dstSy = m_attrs[attr_count++];
            std::vector<float>& dstSz = m_attrs[attr_count++];

            for (size_t i = 0; i < N; ++i) {
              dstSx[i] = srcS[i][0] * scale + offset;
              dstSy[i] = srcS[i][1] * scale + offset;
              dstSz[i] = srcS[i][2] * scale + offset;
            }
        }

        // rot: m_attrs[54, 55, 56, 57]
        {
          auto& desc = outframe.attrDesc[18];
          const float scale  = float(desc.params.attr_scale_minus1 + 1) / float(1 << desc.params.attr_frac_bits);
          const float offset = desc.params.attr_offset;
          const std::vector<Vec3<attr_t>>& srcR_XYZ = pc.getRotationVec();
          std::vector<float>& dstR_W = m_attrs[attr_count++];
          std::vector<float>& dstR_X = m_attrs[attr_count++];
          std::vector<float>& dstR_Y = m_attrs[attr_count++];
          std::vector<float>& dstR_Z = m_attrs[attr_count++];
          for (size_t i = 0; i < N; ++i) {
            dstR_X[i] = srcR_XYZ[i][0] * scale + offset;
            dstR_Y[i] = srcR_XYZ[i][1] * scale + offset;
            dstR_Z[i] = srcR_XYZ[i][2] * scale + offset;
            dstR_W[i] =
              std::sqrt(std::min(std::max(
                1.0 - dstR_X[i] * dstR_X[i] - dstR_Y[i] * dstR_Y[i] - dstR_Z[i] * dstR_Z[i],
                0.0), 1.0));
          }
        }

        // opacity: m_attrs[58]
        {
            auto& desc = outframe.attrDesc[16];
            const float scale  = float(desc.params.attr_scale_minus1 + 1) / float(1 << desc.params.attr_frac_bits);
            const float offset = desc.params.attr_offset;
            const std::vector<attr_t>& srcO = pc.getOpacity();
            std::vector<float>& dstO = m_attrs[attr_count++];
            for (size_t i = 0; i < N; ++i) {
              dstO[i] = srcO[i] * scale + offset;
            }
        }

        outframe.cloud.clear();

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
    std::vector<std::vector<float>> m_attrs;  // positoin, SHs, scales, rots, opacity
    int m_pointCount;
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