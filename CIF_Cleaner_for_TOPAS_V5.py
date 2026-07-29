"""
Author: Meng He
Usage:
1.remove dummy point related redundant info during TOPAS refinement, e.g. "a0"
2.check bond distance based on covalent radii
3.translate (guest molecule) coordinates back within [0,1] [deprecated]
4.remove unreasonable angle (H-H-X) and formatting error from TOPAS "out_bonds_angles" macro, e.g., "invalid"
"""

from funcs import *
from sort import *
import itertools

intro = 'Enter dir of cif file.\n'
cif_dir = input(intro).replace('"', '')

verbose = False

esd_flag = input("Want to check the errors of bond angles/length? Y[yes] or N[no] \n")

if esd_flag in ['Y', 'y']:
    esd_check = True
elif esd_flag in ['N', 'n']:
    esd_check = False
else:
    print('Invalid input, defaulting to check esd.')
    esd_check = True

sort_flag = input('Options: [1] Clean and sort (default); [2] Clean only; [3] Sort only. You choose (number): \n')
if sort_flag in ['1','2']:
    clean_required = True
elif sort_flag == '3':
    clean_required = False
else:
    print('No option selected, processing with defaults: clean and sort.')
    sort_flag = '1'
    clean_required = True

if '.cif' not in cif_dir:
    raise TypeError('not a cif file')
else:
    print(f'loading cif file: {cif_dir} \n')

if os.path.split(cif_dir)[0] != '':
    cif_new_dir = os.path.join(os.path.split(cif_dir)[0], os.path.split(cif_dir)[1].split('.cif')[0] + '_cleaned.cif')
else:
    cif_new_dir = cif_dir.split('.cif')[0] + '_cleaned.cif'
#print(os.path.split(cif_dir)[0])
#print(cif_new_dir)

with open(cif_dir, 'r') as f:
    cif_content = f.readlines()

cif_new_content = []
ln = 0

reg_s = re.compile('\s+')  # split line and return a list of elements

if clean_required:
    print(f'\nCleaning the cif now:\n')
    print(f'Input_File: {cif_dir}')
    print(f'Output_File: {cif_new_dir}')

    while ln < len(cif_content):
        line_content = reg_s.split(cif_content[ln])
        if verbose: print(line_content)
        element_pair = []
        dummy_flag = False
        skip = False

        #0 other checks
        for i in line_content:
            if re.match(r'^(a)([0-9]{0,4})', i) is not None:
                skip = True
                if verbose: print('dummy atoms detected, line skipped')
                break

        if 'Invalid' in line_content:
            if verbose: print('Invalid value')
            skip = True

        if len(line_content) == 9 and line_content[-3] == '0':  # remove 0 occupancy sites
            if verbose: print('zero occupancy')
            skip = True

        if skip is True:
            if verbose: print('line skipped')
            ln = ln + 1
            continue

        #replace oxidation state with simple element symbol, e.g., Fe+2 -> Fe
        for i in range(len(line_content)):
            if re.match(r'^([A-Z][a-z]{0,1})([+-]{1}[0-9]{0,1})', line_content[i]) is not None and '_symmetry_space_group_name_H-M' not in line_content:
                line_content[i] = re.match(r'^([A-Z][a-z]{0,1})([+-]{1}[0-9]{0,1})', line_content[i]).group(1)

        atom_symbol = r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}'
        #1 bond dist check
        if re.match(atom_symbol, line_content[1])  and \
           re.match(atom_symbol, line_content[2])  and \
           re.match(atom_symbol, line_content[3]) is None:

            if esd_check:
                if verbose: print('check error of bond distance...')
                # skip large error due to software issue(should not output certain lines at special conditions)
                if '(' in line_content[3] and ')' in line_content[3]:
                    dist_digits = len(line_content[3].split('(')[0])
                    error_start = line_content[3].index('(')
                    error_end = line_content[3].index(')')
                    error_digits = error_end - error_start - 1
                    if error_digits > 2:
                        if verbose: print('large esd on distance, check if it is the model or the software issue')
                        ln = ln + 1
                        continue
                    elif dist_digits > 6:
                        if verbose: print(
                            'unrealistic accuracy on bond distance, check if it is the model or the software issue')
                        ln = ln + 1
                        continue
                    elif dist_digits - error_digits < 2:
                        if verbose: print('large esd on distance, check if it is the model or the software issue')
                        ln = ln + 1
                        continue
                else:
                    dist_digits = len(line_content[3].split('.')[1])
                    if dist_digits > 6:
                        if verbose: print(
                            'unrealistic accuracy on bond distance, check if it is the model or the software issue')
                        ln = ln + 1
                        continue

            if verbose: print('check bond distance value...')
            element_pair.append(
                re.match(atom_symbol, line_content[1]).group(1))
            element_pair.append(
                re.match(atom_symbol, line_content[2]).group(1))
            if verbose: print(element_pair[0], ' ', element_pair[1], ' ', line_content[3].split('(')[0])
            if element_pair[0] == element_pair[1] == 'H':  #carefully implement if no H2 gas involved
                if verbose: print('skip H-H bond, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[3].split('(')[0]) > bond_dist_max(element_pair):
                if verbose: print('too long, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[3].split('(')[0]) < bond_dist_min_manual(element_pair):
                if verbose: print('too short, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[3].split('(')[0]) == 0:
                if verbose: print('overlapped atom, line skipped')
                ln = ln + 1
                continue
            else:
                if verbose: print('reasonable bond distance, line kept')
                cif_new_content.append(cif_content[ln])
                ln = ln + 1
                continue

        #2 angle check
        if re.match(atom_symbol, line_content[1]) and \
           re.match(atom_symbol, line_content[2]) and \
           re.match(atom_symbol, line_content[3]) is not None:
            redundant = False
            if esd_check:
                if verbose: print('check error of bond angle...')
                # skip large error due to software issue(should not output certain lines at special conditions)
                if '(' in line_content[4] and ')' in line_content[4]:
                    angle_digits = len(line_content[4].split('(')[0])
                    error_start = line_content[4].index('(')
                    error_end = line_content[4].index(')')
                    error_digits = error_end - error_start - 1
                    if '.' in line_content[4]:
                        decimal_points = error_start - line_content[4].index('.') - 1
                    else:
                        decimal_points = 0
                    if verbose: print(f'num of decimal is {decimal_points}')
                    if error_digits > 2:
                        if verbose: print(
                            'large esd on angle, check if it is the model or the software issue, line skipped')
                        ln = ln + 1
                        continue
                    elif angle_digits > 6:
                        if verbose: print(
                            'unrealistic accuracy on angle, check if it is the model or the software issue, line skipped')
                        ln = ln + 1
                        continue
                    elif decimal_points == 0 and angle_digits - error_digits <= 1:
                        if verbose: print(
                            'large esd on angle, check if it is the model or the software issue, line skipped')
                        ln = ln + 1
                        continue

            if verbose: print('check bond angle value...')
            if float(line_content[4].split('(')[0]) == 0:
                if verbose: print('0 bond angle, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[4].split('(')[0]) < 85:  #remove low angle value
                if verbose: print('small angle, line skipped')
                ln = ln + 1
                continue
            elif re.match(atom_symbol, line_content[2]).groups()[0] == 'H' or \
                    re.match(atom_symbol, line_content[2]).groups()[
                        0] == 'D':  #[EVIL FUNCs]
                if verbose: print('angle with H in the middle, line skipped')
                ln = ln + 1
                continue
            elif len(line_content) != 9:  #[EVIL FUNCs]
                if verbose: print('num of line elements not match with loop_')
                ln = ln + 1
                continue
            for pair in itertools.permutations([line_content[1], line_content[2], line_content[3]],
                                               2):  # check same atom site in angle [EVIL FUNCs]
                if pair[0] == pair[1]:
                    redundant = True
                    break
            if redundant is True:
                if verbose: print('not A-B-C angle, atom sites in angle not unique, line skipped')
                ln = ln + 1
                continue
            else:
                if verbose: print('reasonable angle value, line kept')
                cif_new_content.append(cif_content[ln])
                ln = ln + 1
                continue

        cif_new_content.append(cif_content[ln])
        ln = ln + 1
        continue

if clean_required:
    print('New file written.')
    with open(cif_new_dir, 'wt') as f:
        f.writelines(cif_new_content)
        f.close()

#sort things out, the sorted new cif cannot be cleaned using this script again (format changed)
if sort_flag == '1':
    sort(cif_new_dir)
elif sort_flag == '3':
    sort(cif_dir)

