#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     fetch_gene.py
   Description :    用来根据基因列表和文件夹路径提取表达量结果和差异表达结果
   Author :       yueyao
   date：          2019/8/28
-------------------------------------------------
   Change Activity:
                   2019/8/28:
-------------------------------------------------
"""
__author__ = 'yueyao'


#
#GSE59045/4.Comparison/Obesity_vs_NASH.total_genes.xls
#
#GSE59045/2.DEG/all_exp_Obesity-vs-NAFL.xls

#GSE66676/5.Cluster/Healthy_PP_vs_NAFL_PP.fpkm.txt
#GSE66676/2.DEG/all_exp_Obesity-vs-NASH.xls

from genesearch.models import *
import os
import glob
def get_exp_list(path):
    rnaseqdeg = glob.glob(r'' + path + '/*/4.Comparison/*.total_genes.xls')
    rnaseqexp = glob.glob(r'' + path + '/*/5.Cluster/*.fpkm.txt')
    aexp = glob.glob(r'' + path + '/*/2.DEG/all_exp_*-vs-*')
    degfilelist=rnaseqdeg+aexp
    expfilelist=rnaseqexp+aexp
    degfilelist=map(lambda x:x.replace("\\","/"),degfilelist)
    expfilelist=map(lambda x:x.replace("\\","/"),expfilelist)
    return degfilelist,expfilelist

from datetime import datetime
import re
def get_result(gene,degfilelist,expfilelist,path,resultpath):
    DEGs.objects.all().delete()
    now = datetime.now()
    now=now.strftime('%Y%m%d%H')
    DEG=resultpath+"/gene_diff"+now+".xls"
    degfilename ="gene_diff"+now+".xls"
    difftable=[]
    if os.path.exists(DEG):
        os.remove(DEG)
    degf = open(DEG,'a')
    for deg in degfilelist:
        degid = deg.replace(path+"/","")
        if re.search("2\.DEG",degid):
            with open(deg,'r') as f:
                arrheader = f.readline().strip().split("\t")
                tmp=[degid,arrheader[0],arrheader[-3],arrheader[-2],arrheader[-1]]
                difftable.append(tmp)
                DEGs.objects.create(databasetype=degid, gene_id=arrheader[0],fold=arrheader[-3],pvalue=arrheader[-2],padjust=arrheader[-1])
                degf.write(degid+"\t"+"\t".join([arrheader[0],arrheader[-3],arrheader[-2],arrheader[-1]])+"\n")
                for line in f.readlines():
                    arr=line.strip().split("\t")
                    if arr[0].upper() in gene:
                        tmp = [degid, arr[0],arr[-3],arr[-2],arr[-1]]
                        difftable.append(tmp)
                        degf.write(degid+"\t"+"\t".join([arr[0],arr[-3],arr[-2],arr[-1]])+"\n")
                        DEGs.objects.create(databasetype=degid, gene_id=arr[0], fold=arr[-3],
                                            pvalue=arr[-2], padjust=arr[-1])
                difftable.append([])
                degf.write("\n\n")
        if re.search("4\.Comparison",degid):
            with open(deg,'r') as f:
                arrheader = f.readline().strip().split("\t")
                tmp=[degid,arrheader[0],arrheader[-5],arrheader[-2],arrheader[-1]]
                difftable.append(tmp)
                degf.write(degid+"\t"+"\t".join([arrheader[0],arrheader[-5],arrheader[-2],arrheader[-1]])+"\n")
                DEGs.objects.create(databasetype=degid, gene_id=arrheader[0], fold=arrheader[-5], pvalue=arrheader[-2],
                                    padjust=arrheader[-1])
                for line in f.readlines():
                    arr=line.strip().split("\t")
                    if arr[0].upper() in gene:
                        tmp = [degid, arr[0],arr[-5],arr[-2],arr[-1]]
                        difftable.append(tmp)
                        degf.write(degid+"\t"+"\t".join([arr[0],arr[-5],arr[-2],arr[-1]])+"\n")
                        DEGs.objects.create(databasetype=degid, gene_id=arr[0], fold=arr[-5],
                                            pvalue=arr[-2],
                                            padjust=arr[-1])
                difftable.append([])
                degf.write("\n\n")
    degf.close()
    EXP=resultpath+"/gene_exp"+now+".xls"
    expfilename="gene_exp"+now+".xls"
    if os.path.exists(EXP):
        os.remove(EXP)
    expf=open(EXP,'a')
    for exp in expfilelist:
        expid = exp.replace(path + "/", "")
        with open(exp,'r') as f:
            header=f.readline()
            expf.write(expid+"\t"+header)
            for line in f.readlines():
                arr=line.strip().split("\t")
                if arr[0].upper() in gene:
                    expf.write(expid+"\t"+line)
            expf.write("\n\n")
    expf.close()
    return difftable,degfilename,expfilename

if __name__ == '__main__':
    path="E:/IMA/web/imasite/data/database/Liver_cancer"
    resultpath="E:/IMA/web/imasite/result"
    degfilelist, expfilelist=get_exp_list(path)
    gene=["TSPAN6","TNMD","DPM1","SCYL3"]
    gene =map(lambda x:x.upper(),gene)
    get_result(gene,degfilelist, expfilelist,path,resultpath)