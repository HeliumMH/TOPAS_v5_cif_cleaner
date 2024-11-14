import os, read_cif
from collections import namedtuple
from read_cif import __loadCIF__, is_equal, read_number, is_exist_bond, output_loop

Bond = namedtuple("Bond", ["atom1", "atom2", "length", "loc1", "loc2"])
Angle = namedtuple("Angle", ["atom1", "atom2", "atom3", "angle", "loc1", "loc2", "loc3"])



def sort(cif_dir):
    print('\nSorting blocks now:\n')
    if '.cif' not in cif_dir:
        raise TypeError('not a cif file')
    print("Input_File: " + cif_dir)

    if os.path.split(cif_dir)[0] != '':
        cif_new_dir = os.path.join(os.path.split(cif_dir)[0],
                                   os.path.split(cif_dir)[1].split('.cif')[0] + '_sorted.cif')
    else:
        cif_new_dir = cif_dir.split('.cif')[0] + '_sorted.cif'
    print("Output_File: " + cif_new_dir)

    data = __loadCIF__(cif_dir)
    cif_new_content = ""
    for blockName in data.keys():
        block = data[blockName]
        atom_list_n = []
        bond_list_n = []
        bond_list = []
        angle_list_n = []
        atom_dict = dict()
        for i in range(len(block["atom_site_label"])):
            if not is_equal(read_number(block["atom_site_occupancy"][i]), 0):
                atom_dict[block["atom_site_label"][i]] = block["atom_site_type_symbol"][i]
                atom_list_n.append(i)
        for j in range(len(block["geom_bond_atom_site_label_1"])):
            temp_bond = Bond(block["geom_bond_atom_site_label_1"][j], block["geom_bond_atom_site_label_2"][j],
                             read_number(block["geom_bond_distance"][j]), block["geom_bond_site_symmetry_1"][j],
                             block["geom_bond_site_symmetry_2"][j])
            try:
                if read_cif.checkBond(atom_dict[temp_bond.atom1], atom_dict[temp_bond.atom2], temp_bond.length):
                    bond_list_n.append(j)
                    bond_list.append(temp_bond)
            except KeyError:
                continue
        for t in range(len(block["geom_angle_atom_site_label_1"])):
            try:
                temp_angle = Angle(block["geom_angle_atom_site_label_1"][t], block["geom_angle_atom_site_label_2"][t],
                                   block["geom_angle_atom_site_label_3"][t], read_number(block["geom_angle"][t]),
                                   block["geom_angle_site_symmetry_1"][t], block["geom_angle_site_symmetry_2"][t],
                                   block["geom_angle_site_symmetry_3"][t])
            except ValueError:
                print(block["geom_angle_atom_site_label_1"][t], block["geom_angle_atom_site_label_2"][t],
                      block["geom_angle_atom_site_label_3"][t], block["geom_angle"][t],
                      block["geom_angle_site_symmetry_1"][t], block["geom_angle_site_symmetry_2"][t],
                      block["geom_angle_site_symmetry_3"][t])
            if is_exist_bond(temp_angle.atom1, temp_angle.atom2, bond_list) and temp_angle.angle != 0:
                if is_exist_bond(temp_angle.atom2, temp_angle.atom3, bond_list):
                    angle_list_n.append(t)

        block = read_cif.tidy_up(block)
        cif_new_content += blockName + "\n\n"
        for key in block.keys():
            if key == "symmetry_equiv_pos_as_xyz":
                temp_str = "\nloop_\n\t_symmetry_equiv_pos_as_xyz\n"
                for item in block["symmetry_equiv_pos_as_xyz"]:
                    if item[-1] == " ":
                        temp_str += "\t\"" + item[:-1] + "\"\n"
                    else:
                        temp_str += "\t\"" + item + "\"\n"
                temp_str += "\n"
                cif_new_content += temp_str
            elif key == "atom_site_label":
                variables = ["atom_site_label", "atom_site_type_symbol", "atom_site_symmetry_multiplicity",
                             "atom_site_fract_x", "atom_site_fract_y", "atom_site_fract_z", "atom_site_occupancy",
                             "atom_site_B_iso_or_equiv"]
                cif_new_content += output_loop(variables, block, atom_list_n)
            elif key == "geom_bond_atom_site_label_1":
                variables = ["geom_bond_atom_site_label_1", "geom_bond_atom_site_label_2", "geom_bond_distance",
                             "geom_bond_site_symmetry_1", "geom_bond_site_symmetry_2"]
                cif_new_content += output_loop(variables, block, bond_list_n)
            elif key == "geom_angle_atom_site_label_1":
                variables = ["geom_angle_atom_site_label_1", "geom_angle_atom_site_label_2",
                             "geom_angle_atom_site_label_3",
                             "geom_angle", "geom_angle_site_symmetry_1", "geom_angle_site_symmetry_2",
                             "geom_angle_site_symmetry_3"]
                cif_new_content += output_loop(variables, block, angle_list_n)
            elif len(block[key]) <= 3: # any other general keys
                if " " in block[key][0] and ';' not in block[key][0]:
                    cif_new_content += "_%-35s\"%s\"\n" % (key, block[key][0])
                elif ';' in block[key][0]:
                    cif_new_content += ";\n"
                else:
                    cif_new_content += "_%-35s%s\n" % (key, block[key][0])
            elif key == 'pd_meas_2theta_scan': # refined pattern key
                variables = ['pd_meas_2theta_scan', 'pd_proc_intensity_total',
                             'pd_calc_intensity_total', 'pd_proc_ls_weight']
                cif_new_content += output_loop(variables, block, range(len(block[key])))


    with open(cif_new_dir, 'w') as f:
        f.writelines(cif_new_content)
        f.close()

    print("Program finished normally.")


if __name__ == "__main__":
    intro = 'Enter dir of cif file.\n'
    cif_dir = input(intro).replace('"', '')
    sort(cif_dir)
