"""
######################################################################
# PLY read/write utilities for 3D Gaussian Splatting (INRIA format). #
######################################################################

Provides:
    • read3DG_ply()
            parses INRIA-style 3DGS PLY files (positions, SH-DC/AC,
            opacity, scales, rotation).
    • write3DG_ply()
            writes INRIA-style 3DGS datasets to PLY.

Also includes:
    • write_quantization_metadata()
            writes min-max abd bit-depth quantization metadata.
    • read_dequantization_metadata()
            loads metadata for reconstruction.
"""

from pathlib import Path
import numpy as np
import struct
import json

def _tqdm(input):
    return input

def read3DG_ply(file, tqdm=_tqdm):
    """
    Reads a PLY file containing a 3D Gaussian Splat (3DG) dataset encoded in the INRIA format.

    Parameters
    ----------
    file : Path
        Path to the input PLY file.
    tqdm : callable
        Progress-bar function used to display loading progress.

    Returns
    -------
    positions : (N, 3) ndarray
        Gaussian center coordinates.
    sh : (N, M, 3) ndarray
        Spherical harmonics coefficients, where M is the number of SH coefficients.
    opacity : (N,) or (N, 1) ndarray
        Per-Gaussian opacity values.
    scale : (N, 3) ndarray
        Scaling parameters for each Gaussian.
    rotation : (N, 4) ndarray
        Rotation quaternions defining Gaussian orientation.
    comment : str or sequence of str
        Comment line(s) in the PLY header.
    """

    if not file.exists():
        raise Exception(f"Unable to open {file}")
    
    vertex = 0
    faces = 0
    file_format = "format ascii 1.0"
    f_rest_numel = 0
    comments = ''
    with file.open("rb") as fd:
        header_index = []
        while True:
            line = fd.readline().decode("utf-8")

            if "element vertex" in line:
                vertex = int(line.split()[2])
            elif "element face" in line:
                faces = int(line.split()[2])
                while not "end_header" in line:
                    line = fd.readline().decode("utf-8")
                break
            elif "format" in line:
                file_format = line.replace('\n','').replace('\r','')
            elif "property" in line:
                data = line.split()
                if data[1]=="float" or data[1]=="float32":
                    format_type = ['f', 4]
                elif data[1]=="double" or data[1]=="float64":
                    format_type = ['d', 8]
                elif data[1]=="char" or data[1]=="int8":
                    format_type = ['b', 1]
                elif data[1]=="uchar" or data[1]=="uint8":
                    format_type = ['B', 1]
                elif data[1]=="short" or data[1]=="int16":
                    format_type = ['h', 2]
                elif data[1]=="ushort" or data[1]=="uint16":
                    format_type = ['H', 2]
                elif data[1]=="int" or data[1]=="int32":
                    format_type = ['i', 4]
                elif data[1]=="uint" or data[1]=="uint32":
                    format_type = ['I', 4]

                if data[2]=="x":
                    header_index.append(["pos", 0, format_type[0], format_type[1], False])
                elif data[2]=="y":
                    header_index.append(["pos", 1, format_type[0], format_type[1], False])
                elif data[2]=="z":
                    header_index.append(["pos", 2, format_type[0], format_type[1], False])
                elif data[2]=="nx":
                    header_index.append(["normal", 0, format_type[0], format_type[1], False])
                elif data[2]=="ny":
                    header_index.append(["normal", 1, format_type[0], format_type[1], False])
                elif data[2]=="nz":
                    header_index.append(["normal", 2, format_type[0], format_type[1], False])
                elif "f_dc" in data[2]:
                    index = int(data[2].split('_')[-1])
                    header_index.append(["sh_dc", index, format_type[0], format_type[1], False])
                elif "f_rest" in data[2]:
                    index = int(data[2].split('_')[-1])
                    header_index.append(["sh_ac", index, format_type[0], format_type[1], False])
                    f_rest_numel = max(f_rest_numel, index+1)
                elif data[2]=="opacity":
                    header_index.append(["opacity", 0, format_type[0], format_type[1], False])
                elif "scale" in data[2]:
                    index = int(data[2].split('_')[-1])
                    header_index.append(["scale", index, format_type[0], format_type[1], False])
                elif "rot" in data[2]:
                    index = int(data[2].split('_')[-1])
                    header_index.append(["rot", index, format_type[0], format_type[1], False])
                else:
                    header_index.append(["", 0, format_type[0], format_type[1], False])
            elif "comment" in line:
                comments += line
            elif "end_header" in line:
                break

        pos = np.zeros(shape=(vertex,3))
        sh_dc = np.zeros(shape=(vertex,3))
        sh_ac = np.zeros(shape=(vertex,f_rest_numel))
        opacity = np.zeros(shape=(vertex,1))
        scale = np.zeros(shape=(vertex,3))
        rot = np.zeros(shape=(vertex,4))
        
        for n in range(len(header_index)):
            if header_index[n][0]=="pos":
                header_index[n][0] = pos
            elif header_index[n][0]=="sh_dc":
                header_index[n][0] = sh_dc
            elif header_index[n][0]=="sh_ac":
                header_index[n][0] = sh_ac
            elif header_index[n][0]=="scale":
                header_index[n][0] = scale
            elif header_index[n][0]=="rot":
                header_index[n][0] = rot
            elif header_index[n][0]=="opacity":
                header_index[n][0] = opacity
            else:
                header_index[n][4] = True
            if file_format == "format ascii 1.0":
                header_index[n][2] = header_index[n][2]=='f' or header_index[n][2]=='d'
        
        if file_format == "format binary_little_endian 1.0":
            for n in tqdm(range(vertex)):
                for array, k, data_type, byte_count, skip in header_index:
                    value = struct.unpack(data_type, fd.read(byte_count))[0]
                    if skip:
                        continue
                    array[n,k] = value
        elif file_format == "format ascii 1.0":
            for n in tqdm(range(vertex)):
                line = fd.readline().decode("utf-8")
                values = line.split()
                i = 0
                for array, k, isfloat, _, skip in header_index:
                    if skip:
                        i += 1
                        continue
                    if isfloat:
                        value = float(values[i])
                    else:
                        value = int(values[i])
                    array[n,k] = value
        else:
            raise Exception(f"Unsupported {file_format}")
    
    Ncoef = int(f_rest_numel/3)+1
    sh = np.zeros(shape=(vertex, Ncoef, 3))
    for i in range(3):
        sh[:,0,i] = sh_dc[:,i]
        for j in range(1, Ncoef):
            sh[:,j,i] = sh_ac[:,j-1+i*(Ncoef-1)]

    return pos, sh, opacity, scale, rot, comments

def write3DG_ply(pos, sh, opacity, scale, rot, flag_isfloat, file, tqdm=_tqdm, comments=''):
    """
    Writes a PLY file representing a 3D Gaussian Splat (3DG) in th INRIA format.

    Parameters
    ----------
    pos : (N, 3) ndarray
        3DG center positions.
    sh : (N, M, 3) ndarray
        Spherical harmonics coefficients, where M is the number of coefficients.
    opacity : (N,) or (N, 1) ndarray
        Per-Gaussian opacity values.
    scale : (N, 3) ndarray
        Anisotropic scale parameters [s0, s1, s2] for each Gaussian.
    rot : (N, 4) ndarray
        Rotation quaternions defining Gaussian orientation.
    flag_isfloat : bool
        If True, 3DG attributes are written as floating-point values;
        otherwise they are quantized and written as uint16.
    file : Path
        Path to the output PLY file.
    tqdm : callable
        Progress-bar function used to display write progress.
    comment : str or sequence of str
        Comment line(s) in the PLY header.
    """
    
    vertex = pos.shape[0]
    sh_numel = sh.shape[1]*sh.shape[2]

    if flag_isfloat:
        data_type = "float"
    else:
        data_type = "uint16"

    with file.open("wb") as fd:
        fd.write(str.encode("ply\n"))
        if comments != '':
            fd.write(str.encode(comments))
        fd.write(str.encode("format binary_little_endian 1.0\n"))
        fd.write(str.encode(f"element vertex {vertex}\n"))
        for p in ["x", "y", "z"]:
            fd.write(str.encode(f"property float {p}\n"))
        if flag_isfloat:
            for p in ["nx", "ny", "nz"]:
                fd.write(str.encode(f"property float {p}\n"))
        for i in range(sh.shape[2]):
            fd.write(str.encode(f"property {data_type} f_dc_{i}\n"))
        for i in range((sh.shape[1]-1)*sh.shape[2]):
            fd.write(str.encode(f"property {data_type} f_rest_{i}\n"))
        fd.write(str.encode(f"property {data_type} opacity\n"))
        for i in range(3):
            fd.write(str.encode(f"property {data_type} scale_{i}\n"))
        for i in range(4):
            fd.write(str.encode(f"property {data_type} rot_{i}\n"))
        fd.write(str.encode("end_header\n"))

        if flag_isfloat:
            for n in tqdm(range(vertex)):
                for i in range(3):
                    fd.write(struct.pack('f', pos[n,i]))
                for i in range(3):
                    fd.write(struct.pack('f', 0))
                for i in range(sh.shape[2]):
                    fd.write(struct.pack('f', sh[n,0,i]))
                for i in range(sh.shape[2]):
                    for j in range(1, sh.shape[1]):
                        fd.write(struct.pack('f', sh[n,j,i]))
                fd.write(struct.pack('f', opacity[n,0]))
                for i in range(3):
                    fd.write(struct.pack('f', scale[n,i]))
                for i in range(4):
                    fd.write(struct.pack('f', rot[n,i]))
        else:
            for n in tqdm(range(vertex)):
                for i in range(3):
                    fd.write(struct.pack('f', pos[n,i]))
                for i in range(sh.shape[2]):
                    fd.write(struct.pack('H', int(sh[n,0,i])))
                for i in range(sh.shape[2]):
                    for j in range(1, sh.shape[1]):
                        fd.write(struct.pack('H', int(sh[n,j,i])))
                fd.write(struct.pack('H', int(opacity[n,0])))
                for i in range(3):
                    fd.write(struct.pack('H', int(scale[n,i])))
                for i in range(4):
                    fd.write(struct.pack('H', int(rot[n,i])))

##################################
## min-max read/write functions ##
##################################

def write_quantization_metadata(bits, minmax_dict, metadata_filepath):
    """
    Writes quantization metadata JSON file.
    """

    (bits_pos, bits_sh_dc, bits_sh_ac, bits_opacity, bits_scale, bits_rot) = bits

    metadata = {
        "geometry_bitdepth": bits_pos,
        "min_bound_xyz":     minmax_dict["pos_min"],
        "max_bound_xyz":     minmax_dict["pos_max"],

        "sh_dc_bitdepth":    bits_sh_dc,
        "min_bound_sh_dc":   minmax_dict["sh_dc_min"],
        "max_bound_sh_dc":   minmax_dict["sh_dc_max"],

        "sh_ac_bitdepth":    bits_sh_ac,
        "min_bound_sh_ac":   minmax_dict["sh_ac_min"],
        "max_bound_sh_ac":   minmax_dict["sh_ac_max"],

        "opacity_bitdepth":  bits_opacity,
        "min_opacity":       minmax_dict["op_min"],
        "max_opacity":       minmax_dict["op_max"],

        "scales_bitdepth":   bits_scale,
        "min_bound_scales":  minmax_dict["scale_min"],
        "max_bound_scales":  minmax_dict["scale_max"],

        "rotation_bitdepth": bits_rot,
        "min_bound_rotation": minmax_dict["rot_min"],
        "max_bound_rotation": minmax_dict["rot_max"],
    }

    print(f"Saving metadata to {metadata_filepath}...")
    with open(metadata_filepath, "w") as json_file:
        json.dump(metadata, json_file, indent=2)


def read_dequantization_metadata(metadata_filepath):
    """
    Load quantization metadata from JSON file.
    """
    metadata_file_path = Path(metadata_filepath)
    with metadata_file_path.open("r") as fd:
        metadata = json.load(fd)
    return metadata
