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
SHDC_KEYS  = [f"f_dc_{c}" for c in range(3)]
SHAC0_KEYS = [f"f_rest_{15*c + m}" for c in range(3) for m in range(0, 3)]
SHAC1_KEYS = [f"f_rest_{15*c + m}" for c in range(3) for m in range(3, 8)]
SHAC2_KEYS = [f"f_rest_{15*c + m}" for c in range(3) for m in range(8, 15)]
ROT_KEYS   = [f"rot_{c}" for c in range(4)]
SCALE_KEYS = [f"scale_{c}" for c in range(3)]

ATTR_KEYS = SHDC_KEYS + SHAC0_KEYS + SHAC1_KEYS + SHAC2_KEYS + ROT_KEYS + SCALE_KEYS + ["opacity"]
KEY_ORDER = ["positions","f_dc_0","f_dc_1","f_dc_2","f_rest_0","f_rest_1","f_rest_10","f_rest_11","f_rest_12","f_rest_13","f_rest_14","f_rest_15","f_rest_16","f_rest_17","f_rest_18","f_rest_19","f_rest_2","f_rest_20","f_rest_21","f_rest_22","f_rest_23","f_rest_24","f_rest_25","f_rest_26","f_rest_27","f_rest_28","f_rest_29","f_rest_3","f_rest_30","f_rest_31","f_rest_32","f_rest_33","f_rest_34","f_rest_35","f_rest_36","f_rest_37","f_rest_38","f_rest_39","f_rest_4","f_rest_40","f_rest_41","f_rest_42","f_rest_43","f_rest_44","f_rest_5","f_rest_6","f_rest_7","f_rest_8","f_rest_9","opacity","rot_0","rot_1","rot_2","rot_3","scale_0","scale_1","scale_2","Total"]

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
        "RGB-LPIPS": [],
        "GT-RGB-PSNR": [],
        "renT": []
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
    """Return average of list rounded to 8 decimals, or 0.0 if empty."""
    return round(sum(values) / len(values), 8) if values else 0.0

def parse_enc_log(file_path):
    """Parse encoder log file for byte size and encoding times."""
    with open(file_path, "rt") as fp:
        attr_id=0
        for line in fp:
            if "bitstream size" in line:
                token = line.strip().split()
                results["TMC13Info"][KEY_ORDER[attr_id]]["byte"].append(int(token[3]))
            elif "processing time (user)" in line:
                token = line.strip().split()
                results["TMC13Info"][KEY_ORDER[attr_id]]["encT"].append(float(token[-2]))
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
                results["TMC13Info"][KEY_ORDER[attr_id]]["decT"].append(float(token[-2]))
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
    results["MetricInfo"]["YUV-IVSSIM"].append(0.0)
    results["MetricInfo"]["RGB-LPIPS"].append(0.0)
    results["MetricInfo"]["GT-RGB-PSNR"].append(0.0)

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
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SHDC_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SHAC0_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SHAC1_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SHAC2_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in ROT_KEYS),
        sum(sum(results["TMC13Info"][k]["byte"]) for k in SCALE_KEYS),
        sum(results["TMC13Info"]["opacity"]["byte"]),
        "",
        sum(results["TMC13Info"]["Total"]["encT"]),
        sum(results["TMC13Info"]["Total"]["decT"]),
        sum(results["TMC13Info"]["positions"]["encT"]),
        sum(results["TMC13Info"]["positions"]["decT"]),
        sum(sum(results["TMC13Info"][k]["encT"]) for k in ATTR_KEYS),
        sum(sum(results["TMC13Info"][k]["decT"]) for k in ATTR_KEYS),
        max(results["MetricInfo"]["encMaxRSS"], default=0),
        max(results["MetricInfo"]["decMaxRSS"], default=0),
        avg(results["MetricInfo"]["RGB-PSNR"]),
        avg(results["MetricInfo"]["YUV-PSNR"]),
        avg(results["MetricInfo"]["YUV-SSIM"]),
        avg(results["MetricInfo"]["YUV-IVSSIM"]),
        avg(results["MetricInfo"]["RGB-LPIPS"]),
        avg(results["MetricInfo"]["GT-RGB-PSNR"]),
        sum(results["MetricInfo"]["renT"])
    ]

    # Adjust bit counts if not single frame
    if args.flag_1f == 0:
        for n in range(2, 11):
            out[n] = out[n] * 8 * 30 / len(enc_logs) / 1000

    # Output CSV line
    print(",".join(str(v) for v in out))

if __name__ == "__main__":
    main()