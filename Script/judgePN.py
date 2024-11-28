import os
import sys
import pymol
from pymol import cmd

def is_protein_present():
    """Check if the structure contains any protein"""
    protein_selection = cmd.select("protein_selection", "polymer.protein")
    return protein_selection > 0 # 如果选择的有东西的话，将选择print出来是一个大于零的数字，没有的话就是零。

def is_nucleic_acid_present():
    """Check if the structure contains any nucleic acid"""
    nucleic_acid_selection = cmd.select("nucleic_acid_selection", "polymer.nucleic")
    return nucleic_acid_selection > 0

def find_protein_nucleic_acid_structures(pdbs):
    """Find all structures in a directory that contain both protein and nucleic acid"""
    protein_nucleic_acid_structures = []
    print(pdbs)
    # Iterate over all PDB files in the directory
    for pdb_file in pdbs:
        pdb_path = os.path.join(pdb_file, pdb_file+".cif.gz")
        if os.path.exists(pdb_path):
            print(pdb_file)
            try:
                cmd.load(pdb_path, "structure")
            
            # Check for both protein and nucleic acid in the structure
                if is_protein_present() and is_nucleic_acid_present():
                    protein_nucleic_acid_structures.append(pdb_file)
            
            # Clear the structure for the next iteration
                cmd.delete("all")
            except:
                print("error:", pdb_file)
                pass

    return protein_nucleic_acid_structures

# Example usage
pdbs = sys.argv[1]
with open(pdbs) as f:
    f1 = f.readlines()
all_pdbs = f1[0].strip().split(",")
# pymol.finish_launching()  # Start PyMOL
result = find_protein_nucleic_acid_structures(all_pdbs)

# Output the results
rt = open(sys.argv[2], "w")
print("Structures containing both protein and nucleic acid:")
for pdb_file in result:
    rt.write(pdb_file+"\n")
    print(pdb_file)
rt.close()

print("PN:", result)

# Optional: Quit PyMOL after processing
# cmd.quit()

