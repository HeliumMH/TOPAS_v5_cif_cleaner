import pandas as pd
import os
import re
import math

# column 0-8
# atomic_number	symbol	name	empirical 	Calculated	van_der_Waals	Covalent(single_bond)	Covalent (triple bond)	Metallic
element_table_path = r'E:\PythonProjects\src\Covalent_radii_data.csv'
element_table = pd.read_csv(element_table_path)
df = pd.DataFrame(element_table)


# might need to consider H-H bond as a special case
# as well as customized/user-defined atom label like Ow
def bond_dist_max(element_pair, calc_dist=0, bond_param=0.56):
    i = 0
    if element_pair[0]=='H' and element_pair[1]=='H': #[EVIL_FUNS] no H-H, only when not refining H2 adsorption
        max_single_dist=0
        return max_single_dist
    while i < len(element_pair):
        if element_pair[i] == 'Ow':
            element_pair[i] = 'O'
        elif element_pair[i] == 'D':
            element_pair[i] = 'H'
        i += 1
    atom1_idx = df.index[df['symbol'] == element_pair[0]].values[0]  # pd.iloc use 0-based indexing i.e. H-No.0 C-No.5,
    atom2_idx = df.index[df['symbol'] == element_pair[1]].values[0]  # not referring to "atomic number" in the file
    sum_radii = float(df.iloc[atom1_idx, 6]) / 100 + float(df.iloc[atom2_idx, 6]) / 100
    max_single_dist = bond_param + sum_radii
    print(max_single_dist)
    return max_single_dist

def bond_dist_min_manual(element_pair): #[EVIL_FUNCS]
    if element_pair[0]=='C' and element_pair[1]=='C':
        return 1.32
    else:
        return 0.9

 #new_coord=intger_part+negative_fractional_part for negative value integer part +1 when fractional part >0.5
 #for positive value just remove their integer part
 #should conserve same digit
def translate_coords_to_onezero(x,y,z):
    new_coords = []
    esds=[]
    for i in x, y, z:
        if re.match(r"([\-]{0,1}[0-9]+[\.]{0,1}[0-9]*)([\(]{0,1}[0-9]*[\)]{0,1})", i) != None:
            if re.match(r"([\-]{0,1}[0-9]+[\.][0-9]*)([\(]{0,1}[0-9]*[\)]{0,1})", i) != None: #has fraction part
                if re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)(\([0-9]*\))", i) != None : # has negative coords with esd e.g. -1.5(1)
                    print('negative coords with esd')
                    int_part=math.modf(float(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)(\([0-9]+\))", i).group(1)))[1]
                    frac_part=math.modf(float(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)(\([0-9]+\))", i).group(1)))[0]
                    frac_digits=len(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)(\([0-9]+\))", i).group(1))-1-len(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)(\([0-9]+\))", i).group(1).split('.')[0])
                    print(frac_part)
                    print(frac_digits)
                    print(1+frac_part)
                    print(round(1+frac_part,frac_digits))
                    if abs(frac_part) > 0.1 :
                        coord=str(round(1+frac_part,frac_digits))
                        if len(coord.split('.')[1]) <frac_digits:
                            for ext_digit in range(frac_digits-len(coord.split('.')[1])):
                                coord+='0'
                        new_coords.append(coord)
                        esds.append(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]+)(\([0-9]+\))", i).group(2))
                        print(new_coords,esds)
                    else:
                        coord=str(round(frac_part,frac_digits))
                        if len(coord.split('.')[1]) <frac_digits:
                            for ext_digit in range(frac_digits-len(coord.split('.')[1])):
                                coord+='0'
                        new_coords.append(coord)
                        esds.append(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]+)(\([0-9]+\))", i).group(2))
                        print(new_coords,esds)
                elif re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)", i) != None and re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)(\([0-9]+\))", i) ==None: # has negative coords but no esd (wyckoff posit.)
                    print('negative coords but no esd')
                    int_part=math.modf(float(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)", i).group(1)))[1]
                    frac_part=math.modf(float(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)", i).group(1)))[0]
                    frac_digits=len(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)", i).group(1))-1-len(re.match(r"([\-][0-9]+[\.]{0,1}[0-9]*)", i).group(1).split('.')[0])
                    print(abs(frac_part))
                    if abs(frac_part) > 0.4 :
                        coord=str(round(1+frac_part,frac_digits))
                        if len(coord.split('.')[1]) <frac_digits:
                            for ext_digit in range(frac_digits-len(coord.split('.')[1])):
                                coord+='0'
                        new_coords.append(coord)
                        esds.append('')
                        print(new_coords,esds)
                    else:
                        coord=str(round(frac_part,frac_digits))
                        if len(coord.split('.')[1]) <frac_digits:
                            for ext_digit in range(frac_digits-len(coord.split('.')[1])):
                                coord+='0'
                        new_coords.append(coord)
                        esds.append('')
                        print(new_coords,esds)
                elif re.match(r"([0-9]+\.[0-9]*)(\([0-9]*\))", i) != None : #positive coords with esd
                    print('positive coords with esd')
                    int_part=math.modf(float(re.match(r"([0-9]+\.[0-9]*)", i).group(1)))[1]
                    frac_part=math.modf(float(re.match(r"([0-9]+\.[0-9]*)", i).group(1)))[0]
                    frac_digits=len(re.match(r"([0-9]+\.[0-9]*)", i).group(1))-1-len(re.match(r"([0-9]+\.[0-9]*)", i).group(1).split('.')[0])
                    esds.append(re.match(r"([0-9]+\.[0-9]*)(\([0-9]*\))", i).group(2))
                    coord = str(round(frac_part, frac_digits))
                    if len(coord.split('.')[1]) < frac_digits:
                        for ext_digit in range(frac_digits - len(coord.split('.')[1])):
                            coord += '0'
                    new_coords.append(coord)
                    print(new_coords, esds)
                elif   re.match(r"([0-9]+[\.]{0,1}[0-9]*)", i) != None  and re.match(r"([0-9]+[\.]{0,1}[0-9]*)(\([0-9]*\))", i) ==None: #positive coords with no esd (wyckoff posit.)
                    print('positive coords with no esd')
                    int_part=math.modf(float(re.match(r"([0-9]+[\.]{0,1}[0-9]*)", i).group(1)))[1]
                    frac_part=math.modf(float(re.match(r"([0-9]+[\.]{0,1}[0-9]*)", i).group(1)))[0]
                    frac_digits=len(re.match(r"([0-9]+[\.]{0,1}[0-9]*)", i).group(1))-1-len(re.match(r"([0-9]+[\.]{0,1}[0-9]*)", i).group(1).split('.')[0])
                    esds.append('')
                    coord = str(round(frac_part, frac_digits))
                    if len(coord.split('.')[1]) < frac_digits:
                        for ext_digit in range(frac_digits - len(coord.split('.')[1])):
                            coord += '0'
                    new_coords.append(coord)
                    print(new_coords, esds)
            else: #no fractional
                print('no fraction')
                esds.append('')
                new_coords.append(i)
        else:
            raise ValueError('Not a coordinate!')

    return new_coords[0]+esds[0],new_coords[1]+esds[1],new_coords[2]+esds[2]



if __name__ == '__main__':
    #print(bond_dist_max(['H', 'O']))
    #print(translate_coords_to_onezero('5.7502', '6.250(1)', '-7.6251(1)'))
    #print(translate_coords_to_onezero('-5.33', '-6.251(1)', '7.6251'))
    #print(translate_coords_to_onezero('0','0', '-5.678(6)'))
    print(translate_coords_to_onezero('-5.360(4)','-0.250000(12)', '5.178(6)'))


