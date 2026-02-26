# JEE 6.2 & JEE6.6 software for G-PCC-based Anchor Generation 

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

### Notes

#### Build for Release
With single-target generators (e.g. Makefile on linux), cmake may generate
a build system configured for debugging, by default.  To generate a software
more suited for runtime comparisons, it may be a good practice to use a build
system configured for Release. You may add the option
` -D CMAKE_BUILD_TYPE=Release` to the cmake command line, or set the variable
`CMAKE_BUILD_TYPE` to `Release` if you are using a frontend for cmake, for
instance.

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
mpeg-pcc-tmc13-gs/cfg$ ../scripts/gen-cfg.sh --gaussiansplat --all
```

The following example encodes and decodes a GS file named `point_cloud.ply` by setting `r01`.

- `-C`: folder path for output result files. The folder should be created before the encoding.
- `VPATH`: folder path for config files.
- `ENCODER` and `DECODER`: path for the executable file of the software.
- `SRCSEQ`: path for pre-processed input GS file.

```console
mpeg-pcc-tmc13-gs$  make -f $PWD/scripts/Makefile.tmc13-step \
    -C results/r01/ \
    VPATH=$PWD/cfg_3dgs_sample/octree-raht/lossless-geom-lossy-attrs/r01/ \
    ENCODER=$PWD/build/tmc3/tmc3 \
    DECODER=$PWD/build/tmc3/tmc3 \
    SRCSEQ=$PWD/../pre_post_processing/point_cloud.ply
```

When encoding/decoding is done by the above command, the decoded GS file is output as `results/r01/point_cloud.ply.bin.decoded.ply`.

## Evaluation

The specific procedure to evaluate results is provided in Common test conditions for Gaussian splat coding.

## Backward Compatibility

Backward compatibility for mpeg-pcc-tmc13 v12.10 can be tested by running the `../scripts/backward_compatibility/run_verification.py` after setting the appropriate paths within the python script.