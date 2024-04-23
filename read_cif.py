#
# cif.py
#
# Python CIF parser: https://gitlab.com/pdbjapan/tools/cif-parsers
# 
# By Gert-Jan Bekker
# License: MIT
#   See https://gitlab.com/pdbjapan/tools/cif-parsers/blob/master/LICENSE
#
import pandas as pd

# column 0-8
# atomic_number	symbol	name	empirical 	Calculated	van_der_Waals	Covalent(single_bond)	Covalent (triple bond)	Metallic
element_table_path = r'Covalent_radii_data.csv'
element_table = pd.read_csv(element_table_path)
df = pd.DataFrame(element_table)

wt_dict = {'H': 1.008, 'D': 2.014, 'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011, 'N': 14.007, 'O': 15.999,
           'F': 18.998, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.066,
           'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942,
           'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
           'Ga': 69.723, 'Ge': 72.631, 'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 84.798, 'Rb': 84.468,
           'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.906,
           'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.711, 'Sb': 121.760, 'Te': 126.7,
           'I': 126.904, 'Xe': 131.294, 'Cs': 132.905, 'Ba': 137.328, 'La': 138.905, 'Ce': 140.116,
           'Pr': 140.908, 'Nd': 144.243, 'Pm': 144.913, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.925,
           'Dy': 162.500, 'Ho': 164.930, 'Er': 167.259, 'Tm': 168.934, 'Yb': 173.055, 'Lu': 174.967, 'Hf': 178.49,
           'Ta': 180.948, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.085, 'Au': 196.967,
           'Hg': 200.592, 'Tl': 204.383, 'Pb': 207.2, 'Bi': 208.980}


def checkBond(element1, element2, length, min_single_bond_dist=0.9):
    element_pair = [element1, element2]
    max_single_bond_dist = bond_dist_max(element_pair)
    if float(length) > max_single_bond_dist:
        return 0  # not a bond
    elif float(length) < min_single_bond_dist:
        return 0  # not a bond
    else:
        return 1  # a bond


# might need to consider H-H bond as a special case
# as well as customized/user-defined atom label like Ow
def bond_dist_max(element_pair, calc_dist=0, bond_param=0.56):
    i = 0
    while i < len(element_pair):
        if element_pair[i] == 'D':
            element_pair[i] = 'H'
        i += 1
    if element_pair[0] == 'H' and element_pair[1] == 'H':  # [EVIL_FUNS] no H-H, only when not refining H2 adsorption
        max_single_dist = 0
        return max_single_dist
    atom1_idx = df.index[df['symbol'] == element_pair[0]].values[0]  # pd.iloc use 0-based indexing i.e. H-No.0 C-No.5,
    atom2_idx = df.index[df['symbol'] == element_pair[1]].values[0]  # not referring to "atomic number" in the file
    sum_radii = float(df.iloc[atom1_idx, 6]) / 100 + float(df.iloc[atom2_idx, 6]) / 100
    max_single_dist = bond_param + sum_radii
    return max_single_dist


def bond_dist_min_manual(element_pair):  # [EVIL_FUNCS]
    if element_pair[0] == 'C' and element_pair[1] == 'C':
        return 1.32
    else:
        return 0.9


str.partition


def partitionString(string, sep):
    return string.partition(sep)


class _loop:
    def __init__(self, parserObj):
        self.parserObj = parserObj
        self.length = 0
        self.refID = -1
        self.refList = []
        self.namesDefined = False

    def addName(self, name):
        catName = type(name) == str and partitionString(name, ".") or ["", "", ""]
        if catName[1]:
            if not catName[0] in self.parserObj.currentTarget[-2]:
                self.parserObj.currentTarget[-2][catName[0]] = {}
            if not catName[2] in self.parserObj.currentTarget[-2][catName[0]]:
                self.parserObj.currentTarget[-2][catName[0]][catName[2]] = []
                self.refList.append(self.parserObj.currentTarget[-2][catName[0]][catName[2]])
        else:
            if not catName[0] in self.parserObj.currentTarget[-2]:
                self.parserObj.currentTarget[-2][catName[0]] = []
            self.refList.append(self.parserObj.currentTarget[-2][catName[0]])
        self.length = len(self.refList)

    def pushValue(self, value):
        if not self.namesDefined:
            self.namesDefined = True
        target = self.nextTarget()
        if value == "stop_":
            return self.stopPush()
        target.append(value)

    def nextTarget(self):
        self.refID = (self.refID + 1) % self.length
        return self.refList[self.refID]

    def stopPush(self):
        self.refID = -1


def specialSplit(content):
    output = [["", False]]
    quote = False
    length = len(content)
    for c in range(length):
        isWS = content[c] == " " or content[c] == "\t"
        if (content[c] == "'" or content[c] == '"') and (
                c == 0 or content[c - 1] == " " or content[c - 1] == "\t" or c == length - 1 or content[c + 1] == " " or
                content[c + 1] == "\t"):
            quote = not quote
        elif not quote and isWS and output[-1][0] != "":
            output.append(["", False])
        elif not quote and content[c] == "#":
            break
        elif not isWS or quote:
            output[-1][0] += content[c]
            output[-1][1] = quote
    if output[-1][0] == "":
        output.pop()
    return output


class targetSetter:
    def __init__(self, obj, key):
        self.obj = obj
        self.key = key

    def setValue(self, value): self.obj[self.key] = value


class CIFparser:
    def __init__(self):
        self.data = {}
        self.currentTarget = None
        self.loopPointer = None

    def parseString(self, contents):
        multi_line_mode = False
        buffer = []
        for line in contents.splitlines():
            Z = line[:1]
            line = line.strip()
            if Z == ";":
                if multi_line_mode:
                    self.setDataValue("\n".join(buffer))
                else:
                    buffer = []
                multi_line_mode = not multi_line_mode
                line = line[1:].strip()
            if multi_line_mode:
                buffer.append(line)
            else:
                self.processContent(specialSplit(line))

    def parse(self, fileobj):
        multi_line_mode = False
        buffer = []
        for line in fileobj.readlines():
            Z = line[:1]
            line = line.strip()
            if Z == ";":
                if multi_line_mode:
                    self.setDataValue("\n".join(buffer))
                else:
                    buffer = []
                multi_line_mode = not multi_line_mode
                line = line[1:].strip()
            if multi_line_mode:
                buffer.append(line)
            else:
                self.processContent(specialSplit(line))

    def processContent(self, content):
        for c, quoted in content:
            if c == "global_" and not quoted:
                self.loopPointer = None
                self.selectGlobal()
            elif c[:5] == "data_" and not quoted:
                self.loopPointer = None
                self.selectData(c)
            elif c[:5] == "save_" and not quoted:
                self.loopPointer = None
                if c[5:]:
                    self.selectFrame(c)
                else:
                    self.endFrame()
            elif c == "loop_" and not quoted:
                self.loopPointer = _loop(self)
            elif c[:1] == "_" and not quoted:
                self.setDataName(c[1:])
            else:
                self.setDataValue(c)

    def setDataName(self, name):
        if self.loopPointer is not None:
            if self.loopPointer.namesDefined:
                self.loopPointer = None
            else:
                return self.loopPointer.addName(name)
        name = partitionString(name, ".")
        self.currentTarget.pop()
        if name[1]:
            if not name[0] in self.currentTarget[-1]:
                self.currentTarget[-1][name[0]] = {}
            self.currentTarget[-1][name[0]][name[2]] = ""
            self.currentTarget = self.currentTarget + [targetSetter(self.currentTarget[-1][name[0]], name[2])]
        else:
            self.currentTarget[-1][name[0]] = ""
            self.currentTarget = self.currentTarget + [targetSetter(self.currentTarget[-1], name[0])]

    def setDataValue(self, value):
        if self.loopPointer is not None:
            self.loopPointer.pushValue(value)
        else:
            self.currentTarget[-1].setValue([value])

    def selectGlobal(self):
        self.currentTarget = [self.data, self.data, None]

    def selectData(self, name):
        if not name in self.data:
            self.data[name] = {}
        self.currentTarget = [self.data, self.data[name], None]

    def selectFrame(self, name=""):
        if not name in self.currentTarget[1]:
            self.currentTarget[1][name] = {}
        self.currentTarget = self.currentTarget[:2] + [self.currentTarget[1][name], None]

    def endData(self):
        self.currentTarget = self.currentTarget[:2]

    def endFrame(self):
        self.currentTarget = self.currentTarget[:3]


def __loadCIF__(cifFile):
    parser = CIFparser()
    parser.parse(open(cifFile))
    return parser.data


##############################################################################################

def output_loop(variables, block, list_allow):
    cif_str = "\nloop_\n"
    for key in variables:
        cif_str += "\t_" + key + "\n"
    for i in range(len(block[variables[0]])):
        if i not in list_allow:
            continue
        str_line = "\t"
        for key in variables:
            str_line += block[key][i] + "  "
        str_line += "\n"
        cif_str += str_line
    cif_str += "\n"
    return cif_str


def read_number(str_n):
    if "(" in str_n:
        return float(str_n.split("(")[0])
    elif str_n == "Invalid":
        return 0
    else:
        return float(str_n)


def is_equal(float1, float2):
    return abs(float1 - float2) < 0.001


def is_exist_bond(atom1, atom2, list_b):
    for temp_bond in list_b:
        if atom1 == temp_bond.atom1 and atom2 == temp_bond.atom2:
            return True
        elif atom1 == temp_bond.atom2 and atom2 == temp_bond.atom1:
            return True
    return False


def tidy_up(block):
    # crystal_system
    if block["cell_angle_alpha"][0] == "90":
        if block["cell_angle_beta"][0] == "90":
            if block["cell_length_a"] == block["cell_length_b"]:
                if block["cell_length_a"] == block["cell_length_c"]:
                    crystal_system = "cubic"
                else:
                    crystal_system = "tetragonal"
            else:
                crystal_system = "orthorhombic"
        elif block["cell_angle_beta"][0] == "120":
            crystal_system = "hexagonal"
        else:
            crystal_system = "monoclinic"
    else:
        if block["cell_length_a"] == block["cell_length_b"]:
            crystal_system = "rhombohedral"
        else:
            crystal_system = "triclinic"
    block["space_group_crystal_system"][0] = crystal_system

    # formula_sum
    z_num = int(block["cell_formula_units_Z"][0])
    atom_n_dist = {}
    for i in range(len(block["atom_site_label"])):
        if block["atom_site_type_symbol"][i] in atom_n_dist.keys():
            temp_occ = int(block["atom_site_symmetry_multiplicity"][i]) * read_number(block["atom_site_occupancy"][i])
            atom_n_dist[block["atom_site_type_symbol"][i]] += temp_occ
        else:
            temp_occ = int(block["atom_site_symmetry_multiplicity"][i]) * read_number(block["atom_site_occupancy"][i])
            atom_n_dist[block["atom_site_type_symbol"][i]] = temp_occ
    atom_n_dist = dict(sorted(atom_n_dist.items()))
    formula = ""
    for key, value in atom_n_dist.items():
        if value/z_num == int(value/z_num):
            formula += "%s%d " % (key, value/z_num)
        else:
            formula += "%s%.2f " % (key, value/z_num)
    block["chemical_formula_sum"][0] = formula[:-1]
    if "chemical_formula_moiety" in block.keys():
        block["chemical_formula_moiety"][0] = formula[:-1]

    # formula_weight
    formula_wt = 0.0
    for key, value in atom_n_dist.items():
        formula_wt += value * wt_dict[key]
    block["chemical_formula_weight"][0] = "%.2f" % (formula_wt/z_num)

    # density
    density = formula_wt / 6.022 / read_number(block["cell_volume"][0]) * 10
    block["exptl_crystal_density_diffrn"][0] = "%.4f" % (density)

    return block

