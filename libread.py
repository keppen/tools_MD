def readline_2_list(file_content):
    return file_content.readline().strip().split()


def slice_2_type(slice: str, new_type):
    if new_type == "f":
        return float(slice.strip())
    if new_type == "i":
        return int(slice.strip())
    if new_type == "s":
        return slice.strip()
    raise f"ERROR: Invalid type: {new_type}"


# def find_by_item(dictionary, item):
#     """Helper function. Searches list of dictionaries to find dict by id value"""
#     for k, i in dictionary.items():
#         if i == item:
#             return k


def str2float(list):
    return [float(item) for item in list]


def str2int(list):
    return [float(item) for item in list]


def read_xyz(file_content):
    while True:
        n = 0
        atom_coordinates = {}
        atom_names = {}
        while True:
            line = readline_2_list(file_content)
            if len(line) == 1:
                atom_number = int(line[0])
            if len(line) == 2:
                step = int(line[1])
            if len(line) == 4:
                n += 1
                atom_coordinates[str(n)] = str2float(line[1:])
                atom_names[str(n)] = line[0]
                if n == atom_number:
                    break
            if not line:
                break
        yield [step, atom_coordinates, atom_names] if step else [atom_coordinates, atom_names]


# def set_dict_cube(x, y, z):
#     from itertools import product
#     verticies = product(range(x), range(y), range(z))
#     return {key: None for key in verticies}


def read_cube(file_content):
    import numpy as np
    comment = [file_content.readline() for _ in range(2)]
    atom_number, *center = readline_2_list(file_content)
    xvoxels_number, *xaxis_vector = readline_2_list(file_content)
    yvoxels_number, *yaxis_vector = readline_2_list(file_content)
    zvoxels_number, *zaxis_vector = readline_2_list(file_content)
    atom, charge, *coordinates = readline_2_list(file_content)
    cube_dictionary = {
        "xvoxels number": xvoxels_number,
        "yvoxels number": yvoxels_number,
        "zvoxels number": zvoxels_number,
        "xaxis vector": xaxis_vector,
        "yaxis vector": yaxis_vector,
        "zaxis vector": zaxis_vector,
        "comment": comment,
        "atom_number": atom_number,
        "atom": atom,
        "charge": charge,
        "coordinates": coordinates,
    }
    print("Creating dictionary")
    volumetric_data = np.zeros((
        int(xvoxels_number),
        int(yvoxels_number),
        int(zvoxels_number)
    ))
    # volumetric_data = set_dict_cube(int(xvoxels_number),
    #                                 int(yvoxels_number),
    #                                 int(zvoxels_number))
    print("Dictionary Created")
    nx = 0
    ny = 0
    nz = 0
    print("Reading cube...")
    while True:
        line = str2float(readline_2_list(file_content))
        if not line:
            break
        for item in line:
            # print(nx, ny, nz)
            volumetric_data[nx, ny, nz] = item
            nz += 1
            if nz == int(xvoxels_number):
                ny += 1
                nz = 0
            if ny == int(yvoxels_number):
                nx += 1
                ny = 0
    # print(find_by_item(volumetric_data, min(volumetric_data.values())))
    # min = np.where(volumetric_data == np.min(volumetric_data))
    # print(list(zip(min)))
    return volumetric_data, cube_dictionary


def read_pdb(file_content):
    '''Import pdb file into DataFrame. Pdb is output of graph_mol.Molecule class
    It may be truncated to a specified number of models

    Each unique frame is sparated by "MODEL" tag.
    Coordinates tagged by "ATOM" tag and resname "UNK" are collected.
    It is a rigid condition!

    It stores all data in list objects, it is converted into DataFrame object afterwards.
    pandas.df.concate method efficency decreases when volume of DF increases.
    Desipite that one of DF is a single row'''
    # from utils import _printProgress

    # models = []
    atom_serialnumber = []
    atom_name = []
    alternate_location = []
    residue_name = []
    chain_location = []
    residue_seqnumber = []
    code_insertions_residue = []
    x = []
    y = []
    z = []
    occupancy = []
    temperature_factor = []
    segment_identifier = []
    element_symbol = []
    charge = []

    dictionary_data = {

        # "models": models,
        "atom_serialnumber": atom_serialnumber,
        "atom_name": atom_name,
        "alternate_location": alternate_location,
        "residue_name": residue_name,
        "chain_location": chain_location,
        "residue_seqnumber": residue_seqnumber,
        "code_insertions_residue": code_insertions_residue,
        "x": x,
        "y": y,
        "z": z,
        "occupancy": occupancy,
        "temperature_factor": temperature_factor,
        "segment_identifier": segment_identifier,
        "element_symbol": element_symbol,
        "charge": charge
    }

    model = 0

    print('INFO: Reading pdb file and generating DataFrame')
    while True:
        line = file_content.readline()

        if not line:
            break

        if line[:5] == "MODEL":
            model += 1

        # if truncate and model > truncate:
        #     print(
        #         f'INFO: Procedure is truncated to {truncate} models')
        #     break

        # _printProgress(model)

        if line[:4] == 'ATOM':
            # models.append(model)
            atom_serialnumber.append(slice_2_type(line[6:11], 'i'))
            atom_name.append(slice_2_type(line[12:16], 's'))
            alternate_location.append(slice_2_type(line[16], 's'))
            residue_name.append(slice_2_type(line[17:20], 's'))
            chain_location.append(slice_2_type(line[21], 's'))
            residue_seqnumber.append(slice_2_type(line[22:26], 'i'))
            code_insertions_residue.append(slice_2_type(line[27], 's'))
            x.append(slice_2_type(line[30:38], 'f'))
            y.append(slice_2_type(line[38:46], 'f'))
            z.append(slice_2_type(line[46:54], 'f'))
            occupancy.append(slice_2_type(line[54:60], 'f'))
            temperature_factor.append(slice_2_type(line[60:66], 'f'))
            segment_identifier.append(slice_2_type(line[72:77], 's'))
            element_symbol.append(slice_2_type(line[76:78], 's'))
            charge.append(slice_2_type(line[78:80], 's'))

        if line.strip() in ["ENDMDL", "END"]:
            yield model, dictionary_data
            for data in dictionary_data.values():
                data.clear()

    print(f'Collecting model: {model - 1}', '\tINFO: Done', sep='\n')


if __name__ == '__main__':
    from libmath import NumericAnalysis
    import numpy as np
    file = '/home/keppen/CPMD/C3_ISO/META/V.final.cube'
    file_content = open(file, 'r')
    vd, cd = read_cube(file_content)
    file_content.close()
    NumAnal = NumericAnalysis(vd,
                              float(cd["xaxis vector"][0]),
                              np.array([5.102, 3.203, 4.105]),
                              len(cd['xaxis vector']))
    # print(NumAnal.linar_interpolation(NumAnal.point))
    print(NumAnal.steepest_descent())
