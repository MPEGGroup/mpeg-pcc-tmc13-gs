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

import argparse
from pathlib import Path
import openpyxl

START_COL = 11
DEBUG_MODE = False

def write_values(ws, csv_tokens, row):
    """Write a range of CSV values into the worksheet."""
    for col, val in enumerate(csv_tokens):
        ws.cell(
            row=row,
            column=START_COL + col,
            value=float(val) if val else None
        )


def get_offset_dict(ws, flag_1f):
    """Return a mapping from sequence name to row index based on flag."""
    target_flag = "1F" if flag_1f else "NF"
    offset_dict = {}

    for r in range(1, 183):
        f_flag = ws.cell(row=r, column=2).value
        seqname = ws.cell(row=r, column=3).value
        if f_flag == target_flag and seqname:
            offset_dict[seqname] = r
    if DEBUG_MODE:
        print(f"[INFO] Loaded {len(offset_dict)} sequences for flag={target_flag}")
    return offset_dict


def csv_to_xlsm(ws, csv_path, flag_1f):
    """Read a CSV file and write its data into the given worksheet."""
    offset_dict = get_offset_dict(ws, flag_1f)

    with csv_path.open("r", encoding="utf-8") as fp:
        next(fp)  # Skip header
        for line in fp:
            csv_tokens = line.strip().split(",")
            if len(csv_tokens) <= 3:
                continue

            seq = csv_tokens[0]
            if seq not in offset_dict:
                print(f"[WARN] Sequence not found in sheet: {seq}")
                continue
            
            rate_id = 0 if csv_tokens[1] == "lossless" else int(csv_tokens[1].split("r")[-1]) - 1
            row = offset_dict[seq] + rate_id
            
            if DEBUG_MODE:
                print(f"[INFO] Writing seq={seq}, rate_id={rate_id}, row={row}")
            
            write_values(ws, csv_tokens[2:], row)


def process_csv(wb, csv_path, sheet_name, flag_1f):
    """Load a CSV into the given workbook sheet."""
    if not csv_path.exists():
        print(f"[ERROR] CSV file not found: {csv_path}")
        exit(1)

    ws = wb[sheet_name]
    print(f"[INFO] Processing {sheet_name} sheet with CSV: {csv_path.name}")
    
    csv_to_xlsm(ws, csv_path, flag_1f)


def main():
    parser = argparse.ArgumentParser(description="Generate XLSM file from CSV data.")
    parser.add_argument("-f", "--flag_1f", type=int, default=0,
                        help="Flag for single frame processing (1: 1F, 0: NF)")
    parser.add_argument("-t", "--template", required=True, help="Path to XLSM template file")
    parser.add_argument("-x", "--xlsm", required=True, help="Path to output XLSM file")
    parser.add_argument("--aC1", required=False, help="Path to Anchor.C1 CSV file")
    parser.add_argument("--pC1", required=False, help="Path to Proposal.C1 CSV file")
    parser.add_argument("--aCW", required=False, help="Path to Anchor.CW CSV file")
    parser.add_argument("--pCW", required=False, help="Path to Proposal.CW CSV file")
    args = parser.parse_args()

    template_path = Path(args.template)
    wb = openpyxl.load_workbook(template_path, read_only=False, keep_vba=True)
    
    if args.aC1 is not None:
        process_csv(wb, Path(args.aC1), "Anchor.C1", args.flag_1f)
    
    if args.pC1 is not None:
        process_csv(wb, Path(args.pC1), "Proposal.C1", args.flag_1f)
    
    if args.aCW is not None:
        process_csv(wb, Path(args.aCW), "Anchor.CW", args.flag_1f)
    
    if args.pCW is not None:
        process_csv(wb, Path(args.pCW), "Proposal.CW", args.flag_1f)
    
    output_path = Path(args.xlsm)
    wb.save(output_path)
    print(f"[INFO] XLSM file saved: {output_path}")


if __name__ == "__main__":
    main()