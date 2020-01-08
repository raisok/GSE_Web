#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     draw_bar.py
   Description :
   Author :       yueyao
   date：          2019/9/19
-------------------------------------------------
   Change Activity:
                   2019/9/19:
-------------------------------------------------
"""
__author__ = 'yueyao'

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from datetime import datetime
import os

def gettablelist(gene,degfilelist,path,resultpath):
    difftable = []
    now = datetime.now()
    now=now.strftime('%Y%m%d%H')
    for deg in degfilelist:
        degid = deg.replace(path+"/","")
        if re.search("2\.DEG",degid):
            with open(deg,'r') as f:
                for line in f.readlines()[1:]:
                    arr=line.strip().split("\t")
                    #database,name,log2fold,pvalue,qvalue
                    if arr[0].upper() in gene:
                        tmp = [degid,arr[-3],arr[-2],arr[-1]]
                        difftable.append(tmp)
        if re.search("4\.Comparison", degid):
            with open(deg,'r') as f:
                for line in f.readlines()[1:]:
                    arr=line.strip().split("\t")
                    if arr[0].upper() in gene:
                        #database,name,log2fold,pvalue,qvalue
                        tmp = [degid,arr[-5],arr[-2],arr[-1]]
                        difftable.append(tmp)
    picpathname=resultpath+'/'+now+'.output'
    if os.path.exists(picpathname+".output.pdf"):
        os.remove(picpathname+".output.pdf")
    if os.path.exists(picpathname+".output.png"):
        os.remove(picpathname+".output.png")
    genename=gene[0]
    drawpic(difftable,picpathname,genename)
    picfile=os.path.basename(picpathname)
    return picfile


def change_format(ll):
    nn = []
    for i in ll:
        GSEnum = i.split('/')[0]
        gg = i.split('/')[-1].split('.')[0]
        naa = GSEnum + "_" + gg
        nn.append(naa)
    return nn

def drawpic(tablelist,picpathname,genename):
    df=pd.DataFrame(tablelist)
    df.columns=["database","log2fold","pvalue","qvalue"]
    draw = df[['database', 'log2fold']]
    # group = map(lambda x: x.split('/')[0], draw['database'])
    # xnum = range(len(draw['database']))
    logfold = map(lambda x: float(x), draw['log2fold'].tolist())
    ngroup = change_format(draw['database'])
    # fig.size = (25.0, 10.0)
    ax = sns.barplot(x=logfold, y=ngroup, color='blue', orient='h')
    ax.set_title('Gene '+str(genename)+' expression in different database with different group', fontsize=15)
    ax.set_ylabel('Database name', fontsize=15)
    ax.set_xlabel('Log2(Fold change)', fontsize=15)
    # plt.show()
    plt.savefig(picpathname+'.png', format='png', dpi=50, bbox_inches='tight')
    plt.savefig(picpathname+'.pdf', format='pdf', dpi=100, bbox_inches='tight')