#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     fetch_general_gene.py
   Description :
   Author :       yueyao
   date：          2019/10/17
-------------------------------------------------
   Change Activity:
                   2019/10/17:
-------------------------------------------------
"""
__author__ = 'yueyao'

import os
import glob
from datetime import datetime
import re
import logging
import pandas as pd
import random

def v_code():
    ret = ""
    for i in range(6):
        num = random.randint(0, 9)
        # num = chr(random.randint(48,57))#ASCII表示数字
        letter = chr(random.randint(97, 122))#取小写字母
        Letter = chr(random.randint(65, 90))#取大写字母
        s = str(random.choice([num,letter,Letter]))
        ret += s
    return ret

def get_exp_list(path):
    rnaseqdeg = glob.glob(r'' + path + '/*/4.Comparison/*.total_genes.xls')
    rnaseqexp = glob.glob(r'' + path + '/*/5.Cluster/*.fpkm.txt')
    aexp = glob.glob(r'' + path + '/*/2.DEG/all_exp_*-vs-*')
    degfilelist=rnaseqdeg+aexp
    expfilelist=rnaseqdeg+aexp
    degfilelist=map(lambda x:x.replace("\\","/"),degfilelist)
    expfilelist=map(lambda x:x.replace("\\","/"),expfilelist)
    return degfilelist,expfilelist

def get_result(gene,degfilelist,expfilelist,path,resultpath):
    now = datetime.now()
    now=now.strftime('%Y%m%d%H')+str(v_code())
    DEG=resultpath+"/gene_diff"+now+".xls"
    degfilename ="gene_diff"+now+".xls"
    if os.path.exists(DEG):
        os.remove(DEG)
    deg_list = []
    for deg in degfilelist:
        degid = deg.replace(path+"/","")
        GSEid=degid.split('/')[0]
        if re.search("2\.DEG",degid):
            compare_tmp = degid.split('/')[-1].replace("all_exp_","").replace(".xls","")
            control=compare_tmp.split('-vs-')[0]
            treat=compare_tmp.split('-vs-')[1]
            nid=GSEid+"_"+control+"_vs_"+treat
            df = pd.read_csv(deg, header=0, sep="\t")
            df = df.dropna(axis=0, how='any')
            header = list(df.columns)
            conl = [i for i in header if re.match(control, i)]
            treatl = [i for i in header if re.match(treat, i)]
            # df = df[(df['log2FoldChange'] > 0.58496) | (df['log2FoldChange'] < -0.58496)]
            # dffold = df[df['pvalue'] < 0.05]
            dffold=df
            dffold["control_mean"] = dffold[conl].mean(axis=1)
            dffold["treat_mean"] = dffold[treatl].mean(axis=1)
            dffold["Datebase_type"] = nid
            dffold["gene_upper"] = dffold["gene"].str.upper()
            dffold["padj"] = dffold["qvalue"]
            dffold2 = dffold.loc[dffold['gene_upper'].isin(gene)]
            predate = dffold2[["Datebase_type", "gene", "control_mean", "treat_mean", "log2FoldChange","pvalue", "padj"]]
            deg_list.append(predate)
        if re.search("4\.Comparison",degid):
            compare_tmp = degid.split('/')[-1].replace(".total_genes.xls","")
            control=compare_tmp.split('_vs_')[0]
            treat=compare_tmp.split('_vs_')[1]
            nid=GSEid+"_"+control+"_vs_"+treat
            df = pd.read_csv(deg, header=0, sep="\t")
            df = df.dropna(axis=0, how='any')
            df.rename(columns={"gene_id":"gene"},inplace=True)
            header = list(df.columns)
            conl = [i for i in header if re.match(control+"-\d+\(norm\)", i)]
            treatl = [i for i in header if re.match(treat+"-\d+\(norm\)", i)]

            # df = df[(df['log2FoldChange'] > 0.58496) | (df['log2FoldChange'] < -0.58496)]
            # dffold = df[df['padj'] < 0.05]
            dffold = df
            dffold["control_mean"] = dffold[conl].mean(axis=1)
            dffold["treat_mean"] = dffold[treatl].mean(axis=1)
            dffold["Datebase_type"] = nid
            dffold["gene_upper"] = dffold["gene"].str.upper()
            dffold2 = dffold.loc[dffold['gene_upper'].isin(gene)]
            predate = dffold2[["Datebase_type", "gene", "control_mean", "treat_mean", "log2FoldChange","pvalue", "padj"]]
            deg_list.append(predate)
    raw_data = pd.concat(deg_list)
    raw_data.to_csv(DEG,header=True,sep="\t",index=False)
    EXP=resultpath+"/gene_exp"+now+".xls"
    expfilename="gene_exp"+now+".xls"
    if os.path.exists(EXP):
        os.remove(EXP)
    expf=open(EXP,'a')
    for exp in expfilelist:
        expid = exp.replace(path + "/", "")
        with open(exp,'r') as f:
            header=f.readline()
            tmp_list=[]
            tmp_list.append(expid+"\t"+header)
            # expf.write(expid+"\t"+header)
            for line in f.readlines():
                arr=line.strip().split("\t")
                if arr[0].upper() in gene:
                    tmp_list.append(expid + "\t" + line)
                    # expf.write(expid+"\t"+line)
            if len(tmp_list) ==1:
                continue
            else:
                expf.write("".join(tmp_list))
                expf.write("\n\n")
    expf.close()
    logging.warning("DEG file and Expression file successfully create!")
    return degfilename,expfilename

# if __name__ == '__main__':
    # path="E:/IMA/web/imasite/data/database/Liver_cancer"
    # resultpath="E:/IMA/web/imav2/result"
    # degfilelist, expfilelist=get_exp_list(path)
    # gene=("UBD","ZIC2","UHRF1","UBE2T")
    # get_result(gene,degfilelist, expfilelist,path,resultpath)