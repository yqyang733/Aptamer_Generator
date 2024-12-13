import os
import sys
import shutil
import numpy as np
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
            Percent_Buried_back_origin = line_back[7].replace("%", "")
            Percent_Buried_base_origin = line_base[7].replace("%", "")
            if float(Percent_Buried_back_origin) < 0:
                Percent_Buried_back = "0.00"
            else:
                Percent_Buried_back = Percent_Buried_back_origin
            if float(Percent_Buried_base_origin) < 0:
                Percent_Buried_base = "0.00"
            else:
                Percent_Buried_base = Percent_Buried_base_origin
            rt.write(f"{aa},{Chain},{Residue},{Base},{SASA_PRO},{SASA_backbone},{SASA_base},{SASA_complex_back},{SASA_complex_base},{Buried_SASA_back},{Buried_SASA_base},{Percent_Buried_back_origin},{Percent_Buried_base_origin},{Percent_Buried_back},{Percent_Buried_base}\n")
    else:
        log.write(aa+": fail\n")
        log.flush()

def merge_all_in1f(f_in):

    with open(f_in) as f:
        f1 = f.readlines()
    all_file = [i.replace("\n", "") for i in f1]

    rt = open("backbone_base_all.csv", "w")
    rt.write("PDBID,Chain,Residue,Base,SASA_PRO,SASA_backbone,SASA_base,SASA_complex_back,SASA_complex_base,Buried_SASA_back,Buried_SASA_base,Percent_Buried_back_origin,Percent_Buried_base_origin,Percent_Buried_back,Percent_Buried_base\n")
    for aa in all_file:
        with open(os.path.join(aa, "backbone_base.csv")) as f:
            f.readline()
            rt.write("".join(f.readlines()))
    rt.close()

def DNARNA(f_in):

    DNA = ["DA", "DG", "DC", "DT"]
    RNA = ["A", "G", "C", "U"]

    with open(f_in) as f:
        f1 = f.readlines()
    all_file = [i.replace("\n", "") for i in f1]

    rt_dna = open("DNA.csv", "w")
    rt_dna.write("PDBID,Chain,Residue,Base,Percent_Buried_back,Percent_Buried_base\n")
    rt_rna = open("RNA.csv", "w")
    rt_rna.write("PDBID,Chain,Residue,Base,Percent_Buried_back,Percent_Buried_base\n")

    for aa in all_file:
        with open(os.path.join(aa, "backbone_base.csv")) as f:
            f.readline()
            f1 = f.readlines()
        for bb in f1:
            line = bb.replace("\n", "").split(",")
            PDBID = line[0]
            Chain = line[1]
            Residue = line[2]
            Base = line[3]
            Percent_Buried_back = line[13]
            Percent_Buried_base = line[14]
            if Base in DNA:
                rt_dna.write(f"{PDBID},{Chain},{Residue},{Base},{Percent_Buried_back},{Percent_Buried_base}\n")
            elif Base in RNA:
                rt_rna.write(f"{PDBID},{Chain},{Residue},{Base},{Percent_Buried_back},{Percent_Buried_base}\n")
    rt_dna.close()
    rt_rna.close()

def compute_distribution(arr, m):
    # Step 1: 获取区间的终点数组
    min_val, max_val = np.min(arr), np.max(arr)
    bins = np.linspace(min_val, max_val, m + 1)  # m+1 生成m个区间的端点
    bin_centers = (bins[:-1] + bins[1:]) / 2     # 计算区间中点

    # Step 2: 统计每个区间中的元素数量
    hist, _ = np.histogram(arr, bins=bins)

    # Step 3: 对统计数量取以10为底的对数，注意避免对0取对数
    log_hist = np.log10(hist + 1e-10)  # 1e-10 防止分布中存在空区间
    log_hist_final = []
    for i in log_hist:
        if i < 0:
            log_hist_final.append(0)
        elif i == 0:
            log_hist_final.append(0.2)
        else:
            log_hist_final.append(i)

    return bin_centers, np.array(log_hist_final, dtype=np.float32)

def generate_range(n, interval):
    # 计算下界和上界
    lower_bound = (n // interval) * interval
    upper_bound = lower_bound + interval
    
    # 生成范围
    return list(range(0, np.int32(upper_bound + interval), interval))

def plot_scatter(f_in):

    # import seaborn as sns
    import matplotlib.pyplot as plt

    bins = 20

    back = []
    base = []

    with open(f_in) as f:
        f.readline()
        f1 = f.readlines()
    for i in f1:
        line = i.replace("\n", "").split(",")
        back.append(float(line[4]))
        base.append(float(line[5]))
    
    back = np.array(back, dtype=np.float32)
    base = np.array(base, dtype=np.float32)

    top_x, top_y = compute_distribution(back, bins)
    right_x, right_y = compute_distribution(base, bins)

    width_top = np.max(back)/bins
    print("width_top", width_top)
    heigh = np.max(base)/bins
    print("heigh", heigh)

    topfig_x = generate_range(np.max(back), 10)
    topfig_y = generate_range(np.max(top_y), 2)
    rightfig_x = generate_range(np.max(base), 10)
    rightfig_y = generate_range(np.max(right_y), 2)

    border = 0.1
    width = 0.5
    height = 0.2
    between = 0.03

    #散点图的四边
    left_o = border
    bottom_o = border
    height_o = width
    width_o = width

    # 散点图上方的图
    left_t = border
    bottom_t = border+width_o+between
    height_t = height
    width_t = width

    # 散点图右方的图
    left_f = border+width_o+between
    bottom_f = border
    height_f = width
    width_f = height

    # 使用ply.axes需要传入的位置
    ret1 = [left_o,bottom_o,width_o,height_o]
    ret2 = [left_t,bottom_t,width_t,height_t]
    ret3 = [left_f,bottom_f,width_f,height_f]

    fig = plt.figure(figsize=(12,8))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    p1 = plt.axes(ret1)
    p2 = plt.axes(ret2)
    p3 = plt.axes(ret3)

    # 去掉重叠区域的标签
    p2.set_xticks([])
    p3.set_yticks([])

    p1.scatter(back, base, s=2, c="#82B0D2", alpha=1, edgecolors=None)
    p1.set_xlabel('Backbone (Buried Percent %)', fontproperties="Arial",fontsize=18,weight="bold")
    p1.set_ylabel('Base  (Buried Percent %)', fontproperties="Arial",fontsize=18,weight="bold")
    p1.set_xticks(topfig_x)
    p1.set_xticklabels(topfig_x,fontname="Arial",rotation=0,fontsize=16,weight="bold")      # size must be after the font.
    p1.set_yticks(rightfig_x)
    p1.set_yticklabels(rightfig_x,fontname="Arial",fontsize=16,weight="bold")

    p2.bar(top_x,top_y,width=width_top,color="#BEB8DC",edgecolor="black",linewidth=0.3, alpha=0.9)  
    p2.set_ylabel('log(number)', fontproperties="Arial",fontsize=16,weight="bold")
    p2.set_yticks(topfig_y)
    p2.set_yticklabels(topfig_y,fontname="Arial",fontsize=16,weight="bold")

    p3.barh(right_x,right_y,color="#BEB8DC",height=heigh,edgecolor="black",linewidth=0.3, alpha=0.9,)  
    p3.set_xlabel('log(number)', fontproperties="Arial",fontsize=16,weight="bold")
    p3.set_xticks(rightfig_y)
    p3.set_xticklabels(rightfig_y,fontname="Arial",fontsize=16,weight="bold")

    plt.show()
    fig.savefig('huitu.png')

def main():

    file_in = sys.argv[1]

    # DNARNA(file_in)
    plot_scatter(file_in)

    # with open(file_in) as f:
    #     f1 = f.readlines()
    # all_file = [i.replace("\n", "") for i in f1]
    # for aa in all_file:
        # correct_back(aa)
        # merge_back_base(aa)
    # with ProcessPoolExecutor(max_workers=20) as executor:
        # executor.map(correct_back, all_file)  # 并行执行任务
        # executor.map(merge_back_base, all_file)
    # merge_all_in1f(file_in)

if __name__=="__main__":
    main() 