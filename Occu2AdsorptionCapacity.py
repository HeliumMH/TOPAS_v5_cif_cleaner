from datetime import datetime

# convert Rietveld refined occupancy to uptake
MOF_name=input("Name of the structure: \n")
molecule_per_uc=float(input("Number of molecules per unitcell from Riet result: \n"))
MW_uc=float(input("M.W. of desolvated framework: \n"))
uptake=molecule_per_uc/(MW_uc/1000)  # mmol/g
print(f"uptake calculated from refinement result is {uptake} mmol/g")
with open("uptake_log.txt", "a+") as log:
    time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    log.write(time+"\n")
    log.write("MOF\t\tmolecule_per_uc\t\tdesolvated_cell_mass\t\tuptake\n")
    log.write(MOF_name+"\t"+"\t"+str(molecule_per_uc)+"\t"+"\t"+str(MW_uc)+"\t"+"\t"+str(uptake)+"\n")
    log.write("\n")
    log.close()
