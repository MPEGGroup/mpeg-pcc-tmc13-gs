import os
import subprocess
import sys

# ==============================================================================
# CONFIGURATION SECTION
# Please modify the paths below according to your environment.
# ==============================================================================

# Path to the V1 executable
V1_EXEC_PATH = "/path/to/tmc13-v12.10/build/tmc3/tmc3"

# Path to the Target executable
TARGET_EXEC_PATH = "/path/to/tmc13-gs/build/tmc3/tmc3"

# Root directory containing the bitstream files
BITSTREAM_ROOT_DIR = "/path/to/conformance/G-PCCv1/"

# Path to the configuration file required by the decoder
DECODER_CONFIG_PATH = "./cfg/decoder.cfg"

# Path to the text file containing the list of bitstreams and unit lengths
BITSTREAM_LIST_FILE = "bitstreams.txt"

# Command to generate MD5 (e.g., "md5sum", "md5", etc.)
MD5_COMMAND = "md5sum"

# Output directories
OUTPUT_BASE_DIR = "./output/base"
OUTPUT_TARGET_DIR = "./output/target"

# ==============================================================================
# SCRIPT LOGIC
# ==============================================================================

def ensure_directory(path):
    """
    Creates the directory if it does not exist.
    """
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            print(f"[Error] Failed to create directory {path}: {e}")
            sys.exit(1)

def run_command(cmd):
    """
    Executes a shell command.
    """
    try:
        # Using shell=True to support complex command strings
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False

def read_md5_hash(md5_file_path):
    """
    Reads the hash string from the generated .md5 file.
    Assumes standard md5sum output format: "hash  filename"
    """
    try:
        with open(md5_file_path, 'r') as f:
            content = f.read().strip()
            # Split by space and take the first part (the hash)
            return content.split()[0]
    except Exception as e:
        print(f"[Error] Could not read MD5 file {md5_file_path}: {e}")
        return None

def main():
    # 1. Prepare output directories
    ensure_directory(OUTPUT_BASE_DIR)
    ensure_directory(OUTPUT_TARGET_DIR)

    # 2. Read the bitstream list file
    if not os.path.exists(BITSTREAM_LIST_FILE):
        print(f"[Error] Bitstream list file not found: {BITSTREAM_LIST_FILE}")
        sys.exit(1)

    with open(BITSTREAM_LIST_FILE, 'r') as f:
        lines = f.readlines()

    total_count = 0
    match_count = 0
    mismatch_count = 0
    error_count = 0

    print(f"{'='*80}")
    print(f"G-PCC Amd1 for GSC: Starting Backward Compatibility Check")
    print(f"{'='*80}")

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Parse bitstream path and outputUnitLength
        parts = line.split()
        if len(parts) < 2:
            print(f"[Warning] Invalid line format: {line}")
            continue

        rel_bitstream_path = parts[0]
        output_unit_length = parts[1]
        
        # Construct full bitstream path
        full_bitstream_path = os.path.join(BITSTREAM_ROOT_DIR, rel_bitstream_path)
        
        # Extract filename for output naming (e.g., EBS_A_Panasonic_3.bit)
        filename = os.path.basename(rel_bitstream_path)
        
        # Define output filenames
        decoded_ply_name = f"{filename}.decoded.ply"
        md5_name = f"{filename}.decoded.md5"

        # Define full paths for V1 (Base)
        v1_output_ply = os.path.join(OUTPUT_BASE_DIR, decoded_ply_name)
        v1_output_md5 = os.path.join(OUTPUT_BASE_DIR, md5_name)

        # Define full paths for Target
        target_output_ply = os.path.join(OUTPUT_TARGET_DIR, decoded_ply_name)
        target_output_md5 = os.path.join(OUTPUT_TARGET_DIR, md5_name)

        print(f"Processing: {filename} ... ", end="", flush=True)

        # 3. Construct Decoder Commands
        # Command format: executable --config=... --compressedStreamPath=... --reconstructedDataPath=... --outputBinaryPly=1 --outputUnitLength=...
        
        cmd_v1 = (
            f"{V1_EXEC_PATH} "
            f"--config={DECODER_CONFIG_PATH} "
            f"--compressedStreamPath={full_bitstream_path} "
            f"--reconstructedDataPath={v1_output_ply} "
            f"--outputBinaryPly=1 "
            f"--outputUnitLength={output_unit_length}"
        )

        cmd_target = (
            f"{TARGET_EXEC_PATH} "
            f"--config={DECODER_CONFIG_PATH} "
            f"--compressedStreamPath={full_bitstream_path} "
            f"--reconstructedDataPath={target_output_ply} "
            f"--outputBinaryPly=1 "
            f"--outputUnitLength={output_unit_length}"
        )

        # 4. Execute Decoders
        if not run_command(cmd_v1):
            print("FAIL (V1 Execution Error)")
            error_count += 1
            continue
        
        if not run_command(cmd_target):
            print("FAIL (Target Execution Error)")
            error_count += 1
            continue

        # 5. Generate MD5 Checksums
        # Command: md5sum output.ply > output.md5
        cmd_md5_v1 = f"{MD5_COMMAND} {v1_output_ply} > {v1_output_md5}"
        cmd_md5_target = f"{MD5_COMMAND} {target_output_ply} > {target_output_md5}"

        if not run_command(cmd_md5_v1) or not run_command(cmd_md5_target):
            print("FAIL (MD5 Generation Error)")
            error_count += 1
            continue

        # 6. Compare MD5 Results
        hash_v1 = read_md5_hash(v1_output_md5)
        hash_target = read_md5_hash(target_output_md5)

        if hash_v1 and hash_target:
            if hash_v1 == hash_target:
                print("MATCH")
                match_count += 1
            else:
                print("MISMATCH : "+hash_v1+", "+hash_target)
                mismatch_count += 1
        else:
            print("FAIL (MD5 Read Error)")
            error_count += 1
        
        total_count += 1

    # 7. Print Final Summary
    print(f"{'='*80}")
    print(f"Summary Report")
    print(f"{'='*80}")
    print(f"Total Bitstreams Processed : {total_count}")
    print(f"Matches (Pass)             : {match_count}")
    print(f"Mismatches (Fail)          : {mismatch_count}")
    if error_count > 0:
        print(f"Execution Errors           : {error_count}")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()
