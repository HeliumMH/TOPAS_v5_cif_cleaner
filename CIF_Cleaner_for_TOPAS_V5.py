"""
Code name: cif_clean_plus2.0
Author: Meng He (PhD, University of Manchester)
Usage:
1.remove dummy point related redundant info during TOPAS refinement
2.check bond distance based on covalent radii
3.translate (guest molecule) coordinates back within [0,1] [deprecated]
4.remove unreasonable angle (H-H-X) and formatting error from TOPAS "out_bonds_angles" macro, e.g., "invalid"
4.call olex2 for quick check of cif (optional)

"""

import os
import re
from cif_clean_funcs import *
import subprocess
import itertools

__author__ = '2Vague'
__version__ = 'cif_clean_plus2.0'

intro = 'Enter dir of cif file.\n'
cif_dir = input(intro).replace('"', '')
#if '.cif' not in cif_dir:
#    raise TypeError('not a cif file')
print(cif_dir)

if os.path.split(cif_dir)[0] != '':
    cif_new_dir = os.path.join(os.path.split(cif_dir)[0], os.path.split(cif_dir)[1].split('.cif')[0]+'_cleaned.cif')
else:
    cif_new_dir = cif_dir.split('.cif')[0]+'_cleaned.cif'
print(os.path.split(cif_dir)[0])
print(cif_new_dir)

with open(cif_dir, 'r') as f:
    cif_content = f.readlines()

cif_new_content = []
ln = 0

reg_s = re.compile('\s+')  # split line and return a list of elements

while ln < len(cif_content):
    # print(cif_new_content)
    # line_content = []
    line_content = reg_s.split(cif_content[ln])
    print(line_content)
    element_pair = []
    if re.match(r'.*a\d.*', cif_content[ln]) != None:  # remove dummy points during refinement
        for i in line_content:
            if re.match(r'a\d.*', i) != None:
                ln=ln+1
                print('dummy atoms detected, line skipped')
                break
        continue
    elif 'Invalid' in cif_content[ln]:
        print('invalid')
        ln = ln + 1
        continue
    elif len(line_content)==9 and line_content[-3]=='0': # remove 0 occupancy sites
        print('zero occupancy')
        ln=ln+1
        continue

    #bond dist check
    elif re.match(r'[A-Z][a-z]{0,2}[0-9]{0,4}', line_content[1]) != None and \
            re.match(r'[A-Z][a-z]{0,2}[0-9]{0,4}',line_content[2]) != None and \
            re.match(r'[A-Z][a-z]{0,2}[0-9]{0,4}', line_content[3]) == None:
        print('dist. check')

        # skip large error due to software issue(should not output certain lines at special conditions)
        if '(' in line_content[3] and ')' in line_content[3]:
            print('esd check')
            dist_digits=len(line_content[3].split('(')[0])
            error_start=line_content[3].index('(')
            error_end=line_content[3].index(')')
            error_digits=error_end-error_start-1
            if error_digits>2:
                print('large esd on distance, check if it is the model or the software issue')
                ln = ln + 1
                continue
            elif dist_digits>6:
                print('unrealistic accuracy on bond distance, check if it is the model or the software issue')
                ln = ln + 1
                continue
            elif dist_digits-error_digits<2:
                print('large esd on distance, check if it is the model or the software issue')
                ln = ln + 1
                continue

        element_pair.append(re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}', line_content[1]).group(1))
        element_pair.append(re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}', line_content[2]).group(1))
        print(element_pair[0], ' ', element_pair[1], ' ', line_content[3].split('(')[0])
        if float(line_content[3].split('(')[0]) > bond_dist_max(element_pair):
            print('too long')
            ln = ln + 1
            continue
        elif float(line_content[3].split('(')[0])<bond_dist_min_manual(element_pair):
            print('too short')
            ln = ln + 1
            continue
        elif float(line_content[3].split('(')[0]) ==0:
            print('overlapped atom')
            ln = ln + 1
            continue
        elif '(' not in line_content[3]: # skip bond with no esd [EVIL FUNCs]
            print('missing esd, skipped')
            ln = ln + 1
            continue
        else:
            print('okay')
            cif_new_content.append(cif_content[ln])
            ln = ln + 1
            continue


    #angle check
    elif re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}', line_content[1]) and \
            re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}',line_content[2]) and \
            re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}',line_content[3]) !=None:
        print('angle check')
        redundant=False

        # skip large error due to software issue(should not output certain lines at special conditions)
        if '(' in line_content[4] and ')' in line_content[4]:
            angle_digits=len(line_content[4].split('(')[0])
            error_start=line_content[4].index('(')
            error_end=line_content[4].index(')')
            error_digits=error_end-error_start-1
            if error_digits>2:
                print('large esd on angle, check if it is the model or the software issue')
                ln = ln + 1
                continue
            elif angle_digits>6:
                print('unrealistic accuracy on angle, check if it is the model or the software issue')
                ln = ln + 1
                continue
            elif angle_digits-error_digits<3:
                print('large esd on angle, check if it is the model or the software issue')
                ln = ln + 1
                continue

        # angle before esd
        if float(line_content[4].split('(')[0]) ==0:
            print('0 bond angle')
            ln = ln + 1
            continue
        elif float(line_content[4].split('(')[0]) < 90: #remove low angle value
            print('small angle')
            ln = ln + 1
            continue
        elif '(' not in line_content[4]: # skip angle with no esd [EVIL FUNCs]
            print('missing esd, skipped')
            ln = ln + 1
            continue
        elif re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}',line_content[2]).groups()[0] is 'H' or \
                re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,3}',line_content[2]).groups()[0] is 'D': #[EVIL FUNCs]
            print('angle with H in the middle, line skipped')
            ln=ln+1
            continue
        elif len(line_content) !=9: #[EVIL FUNCs]
            print('num of line elements not match with loop_')
            ln=ln+1
            continue
        for pair in itertools.permutations([line_content[1],line_content[2],line_content[3]],2): # check same atom site in angle [EVIL FUNCs]
            if pair[0]==pair[1]:
                redundant=True
                break
        if redundant is True:
            print('not A-B-C angle, atom sites in angle not unique')
            ln=ln+1
            continue
        else:
            print('okay')
            cif_new_content.append(cif_content[ln])
            ln = ln + 1
            continue

    #normalise coordinates(deprecated, use the built-in function "normalize_FCs" in Topas)
    #elif len(line_content)>4 and re.match(r"([\-]{0,1}[0-9]+[\.]{0,1}[0-9]*)([\(]{0,1}[0-9]*[\)]{0,1})", line_content[3]) and re.match(r"([\-]{0,1}[0-9]+[\.]{0,1}[0-9]*)([\(]{0,1}[0-9]*[\)]{0,1})", line_content[4]) and re.match(r"([\-]{0,1}[0-9]+[\.]{0,1}[0-9]*)([\(]{0,1}[0-9]*[\)]{0,1})", line_content[5]) != None:
    #    line_content[3],line_content[4],line_content[5]=translate_coords_to_onezero(line_content[3],line_content[4],line_content[5])
    #    idx=0
    #    cif_content_modified=''
    #    while idx <len(line_content):
    #        cif_content_modified=cif_content_modified+line_content[idx]+' '
    #        idx+=1
    #    print(cif_content_modified)
    #    cif_new_content.append(cif_content_modified+'\n')
    #    ln=ln+1
    #    continue
    else:
        print('no error, line appended to new cif')
        cif_new_content.append(cif_content[ln])
        ln = ln + 1
        continue
with open(cif_new_dir, 'wt') as f:
    f.writelines(cif_new_content)
    f.close()

#run olex2 directly to check the cleaned cif
#make sure you have it and 'olex2' is callable from CMD

#subprocess.run([r"C:\Program Files\Olex2-1.5\olex2.exe",os.path.abspath(cif_new_dir)])
