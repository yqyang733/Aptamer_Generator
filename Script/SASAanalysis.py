import os
import sys
import shutil
from concurrent.futures import ProcessPoolExecutor

def correct_back(aa):

    with open(os.path.join(aa, "backbone_sasa.csv")) as f:
        f.readline()
        f1 = f.readlines()
    rt = open(os.path.join(aa, "backbone_sasa_1.csv"), "w")
    rt.write("Chain,Residue,Base,SASA_PRO,SASA_backbone,SASA_complex,Buried_SASA,Percent_Buried\n")
    for bb in f1:
        line = bb.replace("\n", "").split(",")
        SASA_PRO = line[3]
        # SASA_backbone = line[4]
        SASA_complex = line[5]
        Buried_SASA = line[6]
        SASA_backbone = 2*float(Buried_SASA)+float(SASA_complex)-float(SASA_PRO)
        rt.write(f"{line[0]},{line[1]},{line[2]},{line[3]},{SASA_backbone:.2f},{SASA_complex},{Buried_SASA},{line[7]}\n")
    rt.close()

def merge_back_base(aa):

    log = open("errlog", "w")

    rt = open(os.path.join(aa, "backbone_base.csv"), "w")
    rt.write("PDBID,Chain,Residue,Base,SASA_PRO,SASA_backbone,SASA_base,SASA_complex_back,SASA_complex_base,Buried_SASA_back,Buried_SASA_base,Percent_Buried_back_origin,Percent_Buried_base_origin,Percent_Buried_back,Percent_Buried_base\n")
    with open(os.path.join(aa, "backbone_sasa_1.csv")) as f:
        f.readline()
        f_back = f.readlines()
    with open(os.path.join(aa, "base_sasa.csv")) as f:
        f.readline()
        f_base = f.readlines()
    if len(f_back) == len(f_base):
        log.write(aa+": success\n")
        log.flush()
        for bb in range(len(f_back)):
            line_back = f_back[bb].replace("\n", "").split(",")
            line_base = f_base[bb].replace("\n", "").split(",")
            Chain = line_back[0]
            Residue = line_back[1]
            Base = line_back[2]
            SASA_PRO = line_back[3]
            SASA_backbone = line_back[4]
            SASA_base = line_base[4]
            SASA_complex_back = line_back[5]
            SASA_complex_base = line_base[5]
            Buried_SASA_back = line_back[6]
            Buried_SASA_base = line_base[6]
            Percent_Buried_back_origin = line_back[7]
            Percent_Buried_base_origin = line_base[7]
            if float(Percent_Buried_back_origin) < 0:
                Percent_Buried_back = "0.00"
            else:
                Percent_Buried_back = Percent_Buried_back_origin
            if float(Percent_Buried_base_origin) < 0:
                Percent_Buried_base = "0.00"
            else:
                Percent_Buried_base = Percent_Buried_base_origin
            rt.write(f"{aa},{Chain},{Residue},{Base},{SASA_PRO},{SASA_backbone},{SASA_base},{SASA_complex_back},{SASA_complex_base},{Buried_SASA_back},{Buried_SASA_base},{}\n)
    else:
        log.write(aa+": fail\n")
        log.flush()

def main():

    file_in = sys.argv[1]

    with open(file_in) as f:
        f1 = f.readlines()
    all_file = [i.replace("\n", "") for i in f1]
    # for aa in all_file:
    #     correct_back(aa)
    with ProcessPoolExecutor(max_workers=20) as executor:
        executor.map(correct_back, all_file)  # 并行执行任务

if __name__=="__main__":
    main() 