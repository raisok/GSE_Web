#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     get_deg_exp.py
   Description :
   Author :       yueyao
   date：          2019/9/20
-------------------------------------------------
   Change Activity:
                   2019/9/20:
-------------------------------------------------
"""
__author__ = 'yueyao'


import os
import glob
import logging
def get_exp_list(path):
    rnaseqdeg = glob.glob(r'' + path + '/*/4.Comparison/*.total_genes.xls')
    rnaseqexp = glob.glob(r'' + path + '/*/5.Cluster/*.fpkm.txt')
    aexp = glob.glob(r'' + path + '/*/2.DEG/all_exp_*-vs-*')
    degfilelist=rnaseqdeg+aexp
    expfilelist=rnaseqexp+aexp
    degfilelist=map(lambda x:x.replace("\\","/"),degfilelist)
    expfilelist=map(lambda x:x.replace("\\","/"),expfilelist)
    logging.warning("get expression file done!")
    return degfilelist,expfilelist

from datetime import datetime
import re
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

def get_result(gene,degfilelist,expfilelist,path,resultpath):
    now = datetime.now()
    now=now.strftime('%Y%m%d%H')+str(v_code())
    DEG=resultpath+"/gene_diff"+now+".xls"
    degfilename ="gene_diff"+now+".xls"
    if os.path.exists(DEG):
        os.remove(DEG)
    degf = open(DEG,'a')
    for deg in degfilelist:
        degid = deg.replace(path+"/","")
        if re.search("2\.DEG",degid):
            with open(deg,'r') as f:
                arrheader = f.readline().strip().split("\t")
                tmp=[degid,arrheader[0],arrheader[-3],arrheader[-2],arrheader[-1]]
                degf.write(degid+"\t"+"\t".join([arrheader[0],arrheader[-3],arrheader[-2],arrheader[-1]])+"\n")
                for line in f.readlines():
                    arr=line.strip().split("\t")
                    if arr[0].upper() in gene:
                        tmp = [degid, arr[0],arr[-3],arr[-2],arr[-1]]
                        degf.write(degid+"\t"+"\t".join([arr[0],arr[-3],arr[-2],arr[-1]])+"\n")
                degf.write("\n\n")
        if re.search("4\.Comparison",degid):
            with open(deg,'r') as f:
                arrheader = f.readline().strip().split("\t")
                tmp=[degid,arrheader[0],arrheader[-5],arrheader[-2],arrheader[-1]]
                degf.write(degid+"\t"+"\t".join([arrheader[0],arrheader[-5],arrheader[-2],arrheader[-1]])+"\n")
                for line in f.readlines():
                    arr=line.strip().split("\t")
                    if arr[0].upper() in gene:
                        tmp = [degid, arr[0],arr[-5],arr[-2],arr[-1]]
                        degf.write(degid+"\t"+"\t".join([arr[0],arr[-5],arr[-2],arr[-1]])+"\n")
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
    logging.warning("DEG file and Expression file successfully create!")
    return degfilename,expfilename