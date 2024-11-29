# Aptamer_Generator

（1）将RCSB中的所有结构都下载下来。脚本：[batchdownload.sh](./Script/batchdownload.sh)。下载命令：  
```shell
sh batchdownload.sh -f test.txt -c -p
```
（2）从下载的结构中将同时含有蛋白和核酸的结构挑选出来输出list。脚本：[judgePN.py](./Script/judgePN.py)。  
（3）根据PN的list，从结构中仅将蛋白和核酸保存为cif文件。脚本：[PNprocess.py](./Script/PNprocess.py)。
（4）计算每个cif文件中所有碱基的SASA。脚本：[comSASA.py](./Script/comSASA.py)。为了能计算快速一些，将所有数据分作10份 
