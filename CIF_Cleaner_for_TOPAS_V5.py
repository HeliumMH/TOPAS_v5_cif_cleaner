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

from funcs import *
from sort import *
import itertools



intro = 'Enter dir of cif file.\n'
cif_dir = input(intro).replace('"', '')

esd_flag = input("Want to check the errors of bond angles/length? Y[yes] or N[no] \n")
esd_check = True
if esd_flag == 'Y' or 'y' or '':
    esd_check = True
elif esd_flag == 'N' or 'n':
    esd_check = False
else:
    raise TypeError('Invalid option')

sort_flag = input('Options: [1] Clean and sort (default); [2] Clean only; [3] Sort only. You choose (number): \n')
print(sort_flag)
if sort_flag == '1' :
    sort_required = True
    clean_required = True
elif sort_flag == '2':
    sort_required = False
    clean_required = True
elif sort_flag == '3':
    sort_required = True
    clean_required = False
else:
    print('No option selected, processing with default settings.')
    sort_required = True
    clean_required = True

if '.cif' not in cif_dir:
    raise TypeError('not a cif file')
else:
    print(f'loading cif file: {cif_dir} \n')

if os.path.split(cif_dir)[0] != '':
    cif_new_dir = os.path.join(os.path.split(cif_dir)[0], os.path.split(cif_dir)[1].split('.cif')[0]+'_cleaned.cif')
else:
    cif_new_dir = cif_dir.split('.cif')[0]+'_cleaned.cif'
#print(os.path.split(cif_dir)[0])
#print(cif_new_dir)

with open(cif_dir, 'r') as f:
    cif_content = f.readlines()

cif_new_content = []
ln = 0

reg_s = re.compile('\s+')  # split line and return a list of elements

if clean_required:
    print(f'\nCleaning the cif now:\n')
    print(f'Input_File:{cif_dir}')
    print(f'Input_File:{cif_new_dir}')

    while ln < len(cif_content):
        # print(cif_new_content)
        # line_content = []
        line_content = reg_s.split(cif_content[ln])
        print(line_content)
        element_pair = []
        dummy_flag=False
        skip=False

        for i in line_content:
            if re.match(r'^(a)([0-9]{0,4})(\_[A-Z]{0,1}[a-z]{0,1}[0-9]{0,4}){0,5}', i) != None:
                skip = True
                print('dummy atoms detected, line skipped')
                break
            elif i=='Invalid':
                skip = True
                print('Invalid value')
                break
        if skip is True:
            print('line skipped')
            ln=ln+1
            continue

        if len(line_content)==9 and line_content[-3]=='0': # remove 0 occupancy sites
            print('zero occupancy')
            ln=ln+1
            continue

        #bond dist check
        elif re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}', line_content[1]) != None and \
                re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}',line_content[2]) != None and \
                re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}', line_content[3]) == None:
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
            else:
                dist_digits = len(line_content[3].split('.')[1])
                if dist_digits > 6:
                    print('unrealistic accuracy on bond distance, check if it is the model or the software issue')
                    ln = ln + 1
                    continue
            element_pair.append(re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}', line_content[1]).group(1))
            element_pair.append(re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}', line_content[2]).group(1))
            print(element_pair[0], ' ', element_pair[1], ' ', line_content[3].split('(')[0])
            if element_pair[0]==element_pair[1]=='H': #carefully implement if no H2 gas involved
                print('skip H-H bond, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[3].split('(')[0]) > bond_dist_max(element_pair):
                print('too long, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[3].split('(')[0])<bond_dist_min_manual(element_pair):
                print('too short, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[3].split('(')[0]) ==0:
                print('overlapped atom, line skipped')
                ln = ln + 1
                continue
            elif '(' not in line_content[3]: # skip bond with no esd [EVIL FUNCs]
                if esd_check == True:
                    print('Warning: missing esd, line skipped')
                    ln = ln + 1
                    continue
                elif esd_check == False:
                    print('Warning: missing esd, ignored and continue')
                    cif_new_content.append(cif_content[ln])
                    ln = ln + 1
                    continue
            else:
                print('okay')
                cif_new_content.append(cif_content[ln])
                ln = ln + 1
                continue


        #angle check
        elif re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}', line_content[1]) and \
                re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}',line_content[2]) and \
                re.match(r'^([A-Z][a-z]{0,1})([0-9]{0,4})(\_[0-9]{0,4}){0,5}',line_content[3]) !=None:
            #print('angle check')
            redundant=False

            # skip large error due to software issue(should not output certain lines at special conditions)
            if '(' in line_content[4] and ')' in line_content[4]:
                angle_digits=len(line_content[4].split('(')[0])
                error_start=line_content[4].index('(')
                error_end=line_content[4].index(')')
                error_digits=error_end-error_start-1
                if '.' in line_content[4]:
                    decimal_points = error_start-line_content[4].index('.')-1
                else:
                    decimal_points = 0
                print(f'num of decimal is {decimal_points}')
                if error_digits>2:
                    print('large esd on angle, check if it is the model or the software issue, line skipped')
                    ln = ln + 1
                    continue
                elif angle_digits>6:
                    print('unrealistic accuracy on angle, check if it is the model or the software issue, line skipped')
                    ln = ln + 1
                    continue
                elif decimal_points == 0 and angle_digits-error_digits <= 1:
                    print('large esd on angle, check if it is the model or the software issue, line skipped')
                    ln = ln + 1
                    continue

            # angle before esd
            if float(line_content[4].split('(')[0]) ==0:
                print('0 bond angle, line skipped')
                ln = ln + 1
                continue
            elif float(line_content[4].split('(')[0]) < 85: #remove low angle value
                print('small angle, line skipped')
                ln = ln + 1
                continue
            elif '(' not in line_content[4]: # skip angle with no esd [EVIL FUNCs]
                if esd_check == True:
                    print('Warning: missing esd, line skipped')
                    ln = ln + 1
                    continue
                elif esd_check == False:
                    print('Warning: missing esd, ignored and continue')
                    cif_new_content.append(cif_content[ln])
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
                print('not A-B-C angle, atom sites in angle not unique, line skipped')
                ln=ln+1
                continue
            else:
                print('no error, line appended to new cif')
                cif_new_content.append(cif_content[ln])
                ln = ln + 1
                continue
        else:
            print('no error, line appended to new cif')
            cif_new_content.append(cif_content[ln])
            ln = ln + 1
            continue

if  sort_flag == '1' or sort_flag == '2':
    with open(cif_new_dir, 'wt') as f:
        f.writelines(cif_new_content)
        f.close()

#sort things out, after sort you can not run the previous cleaner because the format of cif change slightly
if sort_flag == '1':
    sort(cif_new_dir)
elif sort_flag == '3':
    sort(cif_dir)
