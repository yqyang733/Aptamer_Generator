import os
import sys
from pymol import cmd

name_f = sys.argv[1]
os.makedirs(name_f)
path_f = f"../PN/{name_f}/{name_f}.cif.gz"
print(path_f)
cmd.load(path_f, "mol")
cmd.create(name_f, 'polymer.protein or polymer.nucleic')
cmd.save(f"{name_f}/{name_f}.cif", name_f)
cmd.delete("all")