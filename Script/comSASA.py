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

def compute_SASA(PDB_FILE):

    # 定义选择
    PROTEIN_SEL = "polymer.protein"
    NUCLEIC_SEL = "polymer.nucleic"
    # NUCLEIC_BASES = "resn A+T+G+C+U"  # 核酸碱基

    OUTPUT_FILE_base = "base_sasa.csv"
    OUTPUT_FILE_backbone = "backbone_sasa.csv"
    # 设置探针半径（通常为1.4 Å）
    probe_radius = 1.4

    print(PDB_FILE)
    # 加载PDB文件
    cmd.load(os.path.join(PDB_FILE, PDB_FILE+".cif.gz"), "mol")
    
    # 获取所有核酸链信息
    nucleic_chains = cmd.get_chains(NUCLEIC_SEL)

    # 打开输出文件
    rt_base = open(os.path.join(PDB_FILE, OUTPUT_FILE_base), "w")
    rt_base.write("Chain,Residue,Base,SASA_PRO,SASA_base,SASA_complex,Buried_SASA,Percent_Buried\n")
    
    rt_backbone = open(os.path.join(PDB_FILE, OUTPUT_FILE_backbone), "w")
    rt_backbone.write("Chain,Residue,Base,SASA_PRO,SASA_backbone,SASA_complex,Buried_SASA,Percent_Buried\n")    
    
    # 计算蛋白的总SASA
    cmd.create("prot", PROTEIN_SEL)
    sasa_pro = cmd.get_area(selection="prot", load_b=1)   
    cmd.delete("prot")        

    # 遍历每条核酸链
    for chain in nucleic_chains:
        chain_sel = f"{NUCLEIC_SEL} and chain {chain}"

        # 获取该链中的所有残基编号和碱基名称
        residues = list()
        cmd.iterate(f"{chain_sel}", "residues.append((chain, resi, resn))", space=locals())
        residues = set(residues)
        residues = sorted(residues, key=lambda x: (x[0], x[1]))

        # 遍历该链上的每个碱基
        for chain, resid, base_name in residues:
            
            base_sel = f"{NUCLEIC_SEL} and chain {chain} and resi {resid} and (not backbone)"
                
            # 计算碱基的总SASA
            cmd.create("base", base_sel)
            sasa_base = cmd.get_area(selection="base", load_b=1)
            cmd.delete("base")
            
            # 计算该碱基与蛋白质结合时的SASA
            # combined_sel = f"({PROTEIN_SEL}) or ({base_sel})"
            # 将蛋白和改碱基保存为一个新的名称为alpha的object。
            pro_base_sel = f"{PROTEIN_SEL} or {base_sel}"
            cmd.create("alpha", pro_base_sel)
            combined_sel = f"alpha"
            sasa_base_complex = cmd.get_area(selection=combined_sel, load_b=1)
            
            # 计算埋藏面积
            buried_sasa = (sasa_pro + sasa_base - sasa_base_complex)/2
            percent_buried = buried_sasa / (sasa_base+0.0000001) * 100

            # 写入输出文件，增加碱基名称
            rt_base.write(f"{chain},{resid},{base_name},{sasa_pro:.2f},{sasa_base:.2f},{sasa_base_complex:.2f},{buried_sasa:.2f},{percent_buried:.2f}%\n")

            cmd.delete("alpha")

            
            backbone_sel = f"{NUCLEIC_SEL} and chain {chain} and resi {resid} and (backbone)"
            
            cmd.create("backb", backbone_sel)
            sasa_backbone = cmd.get_area(selection="backb", load_b=1)
            cmd.delete("backb")
            
            # 计算该碱基与蛋白质结合时的SASA
            # combined_sel = f"({PROTEIN_SEL}) or ({base_sel})"
            # 将蛋白和改碱基保存为一个新的名称为alpha的object。
            pro_backbone_sel = f"{PROTEIN_SEL} or {backbone_sel}"
            cmd.create("beta", pro_backbone_sel)
            combined_sel = f"beta"
            sasa_backbone_complex = cmd.get_area(selection=combined_sel, load_b=1)
            
            # 计算埋藏面积
            buried_sasa = (sasa_pro + sasa_backbone - sasa_backbone_complex)/2
            percent_buried = buried_sasa / (sasa_backbone+0.0000001) * 100

            # 写入输出文件，增加碱基名称
            rt_backbone.write(f"{chain},{resid},{base_name},{sasa_pro:.2f},{sasa_base:.2f},{sasa_backbone_complex:.2f},{buried_sasa:.2f},{percent_buried:.2f}%\n")

            cmd.delete("beta")

    cmd.delete("all")

    rt_backbone.close()
    rt_base.close()

def main():

    file_in = sys.argv[1]
    all_cif_lst = get_lst(file_in)
    random.shuffle(all_cif_lst)
    with ProcessPoolExecutor(max_workers=12) as executor:
        executor.map(compute_SASA, all_cif_lst)  # 并行执行任务

if __name__=="__main__":
    main() 