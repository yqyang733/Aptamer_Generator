import os
import sys
import random
from concurrent.futures import ProcessPoolExecutor
from pymol import cmd

def get_lst(f_in):

    with open(f_in) as f:
        f1 = f.readlines()

    all_cif_lst = []

    for i in f1:
        all_cif_lst.append(i.replace("\n", "").strip())

    return all_cif_lst

def pnprocess(name_f):

    os.makedirs(name_f)
    path_f = f"../PN/{name_f}/{name_f}.cif.gz"
    print(path_f)
    cmd.load(path_f, "mol")
    cmd.create(name_f, 'polymer.protein or polymer.nucleic')
    cmd.save(f"{name_f}/{name_f}.cif.gz", name_f)
    cmd.delete("all")

def main():

    file_in = sys.argv[1]
    all_cif_lst = get_lst(file_in)
    random.shuffle(all_cif_lst)
    with ProcessPoolExecutor(max_workers=12) as executor:
        executor.map(pnprocess, all_cif_lst)  # 并行执行任务

if __name__=="__main__":
    main() 