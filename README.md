# TMC13

## Building

### OSX
- mkdir build
- cd build
- cmake .. -G Xcode 
- open the generated xcode project and build it

### Linux
- mkdir build
- cd build
- cmake .. 
- make

### Windows
- md build
- cd build
- cmake .. -G "Visual Studio 15 2017 Win64"
- open the generated visual studio solution and build it


## Running

This TMC13 codec implementation encodes frame sequences.  A single binary
contains the encoder and decoder implementation, with selection using
the `--mode` option.  Documentation of options is provided via the
`--help` command line option.

### Runtime configuration and configuration files

All command line parameters may be specified in a configuration file.
A set of configuration file templates compliant with the current Common
Test Conditions is provided in the cfg/ directory.

### Example

To generate the configuration files, run the gen-cfg.sh script:

```console
mpeg-pcc-tmc13/cfg$ ../scripts/gen-cfg.sh --all
```

An example script (`scripts/Makefile.tmc13-step`) demonstrates how
to launch the encoder, decoder and metric software for a single
input frame.  The VERBOSE=1 make variable shows the detailed command
execution sequence.  Further documentation of the parameters are
contained within the script.

The following example encodes and decodes frame 0100 of the sequence
`Ford_01_q_1mm`, making use of the configuration file
`cfg/lossy-geom-no-attrs/ford_01_q1mm/r01/encoder.cfg` and storing
the intermediate results in the output directory
`experiment/lossy-geom-no-attrs/ford_01_q1mm/r01/`.

```console
mpeg-pcc-tmc13$ make -f $PWD/scripts/Makefile.tmc13-step \
    -C experiment/lossy-geom-no-attrs/ford_01_q1mm/r01/ \
    VPATH=$PWD/cfg/octree-predlift/lossy-geom-no-attrs/ford_01_q1mm/r01/ \
    ENCODER=$PWD/build/tmc3/tmc3 \
    DECODER=$PWD/build/tmc3/tmc3 \
    PCERROR=/path/to/pc_error \
    SRCSEQ=/path/to/Ford_01_q_1mm/Ford_01_vox1mm-0100.ply \
    NORMSEQ=/path/to/Ford_01_q_1mm/Ford_01_vox1mm-0100.ply

  [encode]  Ford_01_vox1mm-0100.ply.bin <- /path/to/Ford_01_q_1mm/Ford_01_vox1mm-0100.ply
  [md5sum]  Ford_01_vox1mm-0100.ply.bin.md5
  [md5sum]  Ford_01_vox1mm-0100.ply.bin.ply.md5
  [decode]  Ford_01_vox1mm-0100.ply.bin.decoded.ply <- Ford_01_vox1mm-0100.ply.bin
  [md5sum]  Ford_01_vox1mm-0100.ply.bin.decoded.ply.md5
  [metric]  Ford_01_vox1mm-0100.ply.bin.decoded.ply.pc_error <- Ford_01_vox1mm-0100.ply.bin.decoded.ply
```

## 3D Gaussian Coding (3DGC) support by KDDI
### Samples of cfg files
Samples of the cfg files are stored in the folder `cfg_3dgs_sample/`.

#### encoder.cfg
- The `3DGS` option shall be set to 1 for 3D Gaussian coding.
- For INRIA's pre-train data, we set `inputScale: 1000` because geometrical information is in meter order. 
- Additionally, qp, bit-depth, and attribute names shall be set for each attribute.

#### decoder.cfg
- Attribute name shall be set for each attribute.
- For INRIA's pre-train data, we set `outputUnitLength: 1000`. 

### Pre-processing in the software
#### Scaling
Input values are scaled by a fixed value (i.e., 1000) in `ply::read()`. 

```cpp
          ifs.read(reinterpret_cast<char*>(&buf), sizeof(float));
          gs_t attr = buf * GS_PROPERTY_SCALE;
          cloud.getGSInfo(pointCounter,key) = gs_t(attr);
```

### Add Offset
Offset values are derived as the minimum value for each attribute. Offset values are sent to the decoder as attr_offset.
```cpp
  if (codeGSInfo)
 {
    for(auto& desc : params->encoder.sps.attributeSets){
      auto key = GSHelper::getGSKey(desc.attr_instance_id);
      auto& gsinfo = pointCloud.getGSInfo(key);
      gs_t min = *std::min_element(gsinfo.begin(),gsinfo.end());
      desc.params.scalingParametersPresent = true;
      desc.params.attr_offset = min;
    }
 }
```

### Encoding/Decoding
All attributes are encoded/decoded as single-channel information. Encoding/decoding processes for reflectance are used.

### Limitation
- KDDI has run the software only in the case of Octree-RAHT lossless-geom-lossy-attr.
- KDDI used [INRIA's pre-trained 3DGS data](https://repo-sam.inria.fr/fungraph/3d-gaussian-splatting/datasets/pretrained/models.zip) such as Bicycle etc.
- Parameters (QP, bit-depth, scaling factor, and offset value) are not tuned. It would be appreciated if you tuned those parameters.
