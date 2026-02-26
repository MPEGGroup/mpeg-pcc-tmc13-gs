# The copyright in this software is being made available under the BSD
# Licence, included below.  This software may be subject to other third
# party and contributor rights, including patent rights, and no such
# rights are granted under this licence.
#
# Copyright (c) 2025, ISO/IEC
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the ISO/IEC nor the names of its contributors
#   may be used to endorse or promote products derived from this
#   software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import glob
import argparse
import os

# ===== Constants =====
SH_B0_KEYS = ["SH_dc_coefficients_0"]
SH_B1_KEYS = [f"SH_ac_coefficients_{m}" for m in range(1,4)]
SH_B2_KEYS = [f"SH_ac_coefficients_{m}" for m in range(4,9)]
SH_B3_KEYS = [f"SH_ac_coefficients_{m}" for m in range(9,16)]
ROT_KEYS   = ["rotations_XYZ", "rotations_W"]
SCALE_KEYS = ["scales"]

ATTR_KEYS = SH_B0_KEYS + SH_B1_KEYS + SH_B2_KEYS + SH_B3_KEYS + ROT_KEYS + SCALE_KEYS + ["opacities"]
ATTR_KEY_ORDER = ["positions","SH_dc_coefficients_0","SH_ac_coefficients_1","SH_ac_coefficients_10","SH_ac_coefficients_11","SH_ac_coefficients_12","SH_ac_coefficients_13","SH_ac_coefficients_14","SH_ac_coefficients_15","SH_ac_coefficients_2","SH_ac_coefficients_3","SH_ac_coefficients_4","SH_ac_coefficients_5","SH_ac_coefficients_6","SH_ac_coefficients_7","SH_ac_coefficients_8","SH_ac_coefficients_9","opacities","rotations_XYZ","rotations_W","scales"]

# ===== Results structure =====
results = {
    "TMC13Info": {},
    "MetricInfo": {
        "encMaxRSS": [],
        "decMaxRSS": [],
        "RGB-PSNR": [],
        "YUV-PSNR": [],
        "YUV-SSIM": [],
        "YUV-IVSSIM": [],
        "GT-RGB-PSNR": []
    }
}

def register_keys(keys, results_dict):
    """Initialize TMC13Info entries for given keys."""
    for key in keys:
        results_dict["TMC13Info"][key] = {"byte": [], "encT": [], "decT": []}

# Register all required keys
register_keys(["positions"], results)
register_keys(ATTR_KEYS, results)
register_keys(["Total"], results)

def avg(values):
    """Return average of list rounded to 12 decimals, or 0.0 if empty."""
    return round(sum(values) / len(values), 12) if values else 0.0

ATTR_KEY_SIZE=len(ATTR_KEY_ORDER)
def parse_enc_log(file_path):
    """Parse encoder log file for byte size and encoding times."""
    with open(file_path, "rt") as fp:
        attr_id=0
        for line in fp:
            if "Total bitstream size" in line:
                token = line.strip().split()    
                results["TMC13Info"]["Total"]["byte"].append(int(token[3]))
            elif "bitstream size" in line:
                token = line.strip().split()
                results["TMC13Info"][ATTR_KEY_ORDER[attr_id%ATTR_KEY_SIZE]]["byte"].append(int(token[3]))
            elif "processing time (user)" in line:
                token = line.strip().split()
                results["TMC13Info"][ATTR_KEY_ORDER[attr_id%ATTR_KEY_SIZE]]["encT"].append(float(token[-2]))
                attr_id+=1
            elif "Processing time (user)" in line:
                token = line.strip().split()
                results["TMC13Info"]["Total"]["encT"].append(float(token[-2]))

def parse_dec_log(file_path):
    """Parse decoder log file for decoding times."""
    with open(file_path, "rt") as fp:
        attr_id=0
        for line in fp:
            if "processing time (user)" in line:
                token = line.strip().split()
                results["TMC13Info"][ATTR_KEY_ORDER[attr_id % ATTR_KEY_SIZE]]["decT"].append(float(token[-2]))
                attr_id+=1
            elif "Processing time (user)" in line:
                token = line.strip().split()
                results["TMC13Info"]["Total"]["decT"].append(float(token[-2]))

def parse_rss_file(file_path, key):
    """Parse RSS file for max resident set size."""
    with open(file_path, "rt") as fp:
        token = fp.readline().strip().split()[-1]
        max_rss = int(token.split("maxresident)k")[0])
        results["MetricInfo"][key].append(max_rss)

def parse_metric_file(file_path):
    """Parse metric file for various quality metrics."""
    # Initialize default values
    results["MetricInfo"]["YUV-IVSSIM"].append(0)
    results["MetricInfo"]["GT-RGB-PSNR"].append(0)

    metric_map = {
        "Psnr RGB (avg)": "RGB-PSNR",
        "Psnr YUV (avg)": "YUV-PSNR",
        "SSIM (avg)": "YUV-SSIM"
    }

    with open(file_path, "rt") as fp:
        for line in fp:
            for search_str, metric_key in metric_map.items():
                if search_str in line and "OM-" not in line:
                    results["MetricInfo"][metric_key].append(float(line.strip().split()[-1]))

def fmt(v):
    if isinstance(v, int):
        return str(v)                 
    elif isinstance(v, float):
        return f"{v:.12f}"            
    else:
        return str(v)                 

def main():
    parser = argparse.ArgumentParser(description="TMC13 & MPEG-GSC-METRIC results parser")
    parser.add_argument("-d", "--dir", required=True, help="Directory containing log/metric files")
    parser.add_argument("-s", "--seq", required=True, help="Sequence name")
    parser.add_argument("-r", "--rate", required=True, help="Rate string (e.g., r01)")
    parser.add_argument("-f", "--flag_1f", type=int, default=0, help="Flag for single frame processing")
    args = parser.parse_args()

    # Collect file lists
    enc_logs = sorted(glob.glob(f"{args.dir}/*.bin.log"))
    dec_logs = sorted(glob.glob(f"{args.dir}/*.bin.decoded.log"))
    enc_rss  = sorted(glob.glob(f"{args.dir}/*.bin.time"))
    dec_rss  = sorted(glob.glob(f"{args.dir}/*.bin.decoded.time"))
    metrics  = sorted(glob.glob(f"{args.dir}/*.metric"))

    # Sanity check
    if len(enc_logs) != len(dec_logs):
        print(f"[Error] Encoder/decoder log count mismatch in {args.dir}")
        exit(1)

    # Parse files
    for f in enc_logs: parse_enc_log(f)
    for f in dec_logs: parse_dec_log(f)
    for f in enc_rss:  parse_rss_file(f, "encMaxRSS")
    for f in dec_rss:  parse_rss_file(f, "decMaxRSS")
    for f in metrics:  parse_metric_file(f)

    # Prepare output
    out = [
        args.seq,
        args.rate,
        sum(results["TMC13Info"]["Total"]["byte"]),
        sum(results["TMC13Info"]["positions"]["byte"]),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SH_B0_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SH_B1_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SH_B2_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SH_B3_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in ROT_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SCALE_KEYS),
        sum(results["TMC13Info"]["opacities"]["byte"]),
        avg(results["MetricInfo"]["RGB-PSNR"]),
        avg(results["MetricInfo"]["YUV-PSNR"]),
        avg(results["MetricInfo"]["YUV-SSIM"]),
        avg(results["MetricInfo"]["YUV-IVSSIM"]),
        avg(results["MetricInfo"]["GT-RGB-PSNR"]),
        
        min(results["MetricInfo"]["RGB-PSNR"], default=0),
        min(results["MetricInfo"]["YUV-PSNR"], default=0),
        min(results["MetricInfo"]["YUV-SSIM"], default=0),
        
        max(results["MetricInfo"]["RGB-PSNR"], default=0),
        max(results["MetricInfo"]["YUV-PSNR"], default=0),
        max(results["MetricInfo"]["YUV-SSIM"], default=0),

        sum(results["TMC13Info"]["Total"]["encT"]),
        sum(results["TMC13Info"]["positions"]["encT"]),
        sum(sum(results["TMC13Info"][k]["encT"]) for k in ATTR_KEYS),

        sum(results["TMC13Info"]["Total"]["decT"]),
        sum(results["TMC13Info"]["positions"]["decT"]),
        sum(sum(results["TMC13Info"][k]["decT"]) for k in ATTR_KEYS),
        
        max(results["MetricInfo"]["encMaxRSS"], default=0),
        max(results["MetricInfo"]["decMaxRSS"], default=0),
    ]

    # Output CSV line
    print(",".join(fmt(v) for v in out))

if __name__ == "__main__":
    main()