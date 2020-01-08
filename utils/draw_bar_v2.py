#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     draw_bar_v2.py
   Description :
   Author :       yueyao
   date：          2019/10/18
-------------------------------------------------
   Change Activity:
                   2019/10/18:
-------------------------------------------------
"""
__author__ = 'yueyao'


import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from datetime import datetime
import os
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

def gettablelist(gene,degfilerfile,resultpath):
    now = datetime.now()
    now=now.strftime('%Y%m%d%H')+str(v_code())
    picpathname=resultpath+'/'+now+'.output'
    if os.path.exists(picpathname+".output.pdf"):
        os.remove(picpathname+".output.pdf")
    if os.path.exists(picpathname+".output.png"):
        os.remove(picpathname+".output.png")
    genename=gene[0]
    # drawpic(degfilerfile,picpathname,genename)
    drawpic_v2(degfilerfile,picpathname,genename)
    picfile=os.path.basename(picpathname)
    return picfile


def drawpic(degfilerfile,picpathname,genename):
    df=pd.read_csv(degfilerfile,header=0,sep="\t")
    # 根据绝对值大小从大到小进行排序
    df = df.reindex(df['log2FoldChange'].abs().sort_values(ascending=False).index)
    # 取前20行
    df = df.iloc[:20, :]
    # 获取最大的fold绝对值作为x轴的上下限
    txlim = np.abs(df.loc[:, 'log2FoldChange'].tolist()[0])
    draw = df[["Datebase_type", "log2FoldChange"]]
    draw = draw[np.abs(draw['log2FoldChange']) > 0.58496]
    df['colors'] = ['blue' if x < 0 else 'red' for x in df['log2FoldChange']]
    logfold = [float(n) for n in draw['log2FoldChange'].tolist()]
    ngroup = [str(n) for n in draw['Datebase_type'].tolist()]
    ax = sns.barplot(x=logfold, y=ngroup, color='blue', orient='h')
    ax.set_title('20 Top FoldChange of Gene '+str(genename)+' expression in different database with different group', fontsize=15)
    ax.set_ylabel('Datebase name', fontsize=15)
    ax.set_xlabel('Log2(Fold change)', fontsize=15)
    ax.set_xlim([-txlim - 0.25, txlim + 0.25])
    plt.savefig(picpathname+'.png', format='png', dpi=50, bbox_inches='tight')
    plt.savefig(picpathname+'.pdf', format='pdf', dpi=100, bbox_inches='tight')
    plt.close()

def drawpic_v1(degfilerfile,picpathname,genename):

    df=pd.read_csv(degfilerfile,header=0,sep="\t")
    # 选择绝对值大于0.58496的数据
    #df = df[np.abs(df['log2FoldChange']) > 0.58496]
    # 根据绝对值大小从大到小进行排序
    df = df.reindex(df['log2FoldChange'].abs().sort_values(ascending=False).index)
    # 取前20行
    df = df.iloc[:20, :]
    # 获取最大的fold绝对值作为x轴的上下限
    txlim = np.abs(df.loc[:, 'log2FoldChange'].tolist()[0])
    #获取作图的两列
    draw = df[["Datebase_type", "log2FoldChange"]]
    print(draw.head())
    draw['colors'] = ['blue' if x < 0 else 'red' for x in draw['log2FoldChange']]
    plt.hlines(y=draw.Datebase_type, xmin=0, xmax=draw.log2FoldChange,
               color=draw.colors, linewidth=10)
    plt.gca().set(ylabel='Datebase name', xlabel='Log2(Fold change)')
    plt.title('20 Top FoldChange of Gene ' + str(genename) + ' expression in different dataset')
    plt.yticks(fontsize=8)
    plt.xlim(-txlim - 0.25, txlim + 0.25)

    plt.savefig(picpathname+'.png', format='png', dpi=100, bbox_inches='tight')
    plt.savefig(picpathname+'.pdf', format='pdf', dpi=100, bbox_inches='tight')
    plt.close()

#update at 20191213 by yueyao
def drawpic_v2(degfilerfile,picpathname,genename):

    df=pd.read_csv(degfilerfile,header=0,sep="\t")
    # 根据绝对值大小从大到小进行排序
    df = df.reindex(df['log2FoldChange'].abs().sort_values(ascending=False).index)
    # 取前20行
    df = df.iloc[:20, :]
    # 获取最大的fold绝对值作为x轴的上下限
    txlim = np.abs(df.loc[:, 'log2FoldChange'].tolist()[0])
    #获取作图需要的列
    draw = df[["Datebase_type", "log2FoldChange", "pvalue"]]
    #对新数据框从小到大排序，方便作图展示时的顺序
    draw = draw.reindex(draw['log2FoldChange'].abs().sort_values(ascending=True).index)
    #设置每个柱形需要填充的颜色
    draw['color'] = ['b' if x < 0 else 'r' for x in draw['log2FoldChange']]
    #设置每个柱形边框的颜色，p值小于0.05设置为k，也就是黑色
    draw['edge_color'] = ['k' if x < 0.05 else 'none' for x in draw['pvalue']]
    draw['hatch_color'] = ['/\\' if x < 0.05 else 'none' for x in draw['pvalue']]

    print("image is drawing....")
    for j,i in enumerate(draw.index):
        plt.barh(j,draw["log2FoldChange"][i],color=draw["color"][i],ec=draw['edge_color'][i],ls='--',lw=1,hatch=draw['hatch_color'][i])
    # plt.barh(range(len(draw["log2FoldChange"])), draw["log2FoldChange"], color=draw["color"],
    #          tick_label=draw["Datebase_type"], ec=draw['edge_color'], ls='--', lw=1)

    plt.gca().set(ylabel='Datebase name', xlabel='Log2(Fold change)')
    plt.title('20 Top FoldChange of Gene ' + str(genename) + ' expression in different dataset')
    plt.yticks(range(len(draw["log2FoldChange"])),draw["Datebase_type"],fontsize=8)

    plt.xlim(-txlim - 0.25, txlim + 0.25)

    #单独设置图例
    red_patch = mpatches.Patch(color='red', label='Up')
    blue_patch = mpatches.Patch(color='blue',label='Down')
    hatch_path = mpatches.Patch(color='white',ec='k',ls='--',hatch='/\\',label='Significant')
    plt.legend(handles=[red_patch,blue_patch,hatch_path],loc="best")

    plt.savefig(picpathname+'.png', format='png', dpi=100, bbox_inches='tight')
    plt.savefig(picpathname+'.pdf', format='pdf', dpi=100, bbox_inches='tight')
    plt.close()


#update at 20191220 by yueyao
def draw_exp_barplot(expfile,resultpath):
    filename = os.path.basename(expfile)
    pic_file_name = filename.replace(".xls",".pdf")
    cmd = "Rscript /home/user/Web/imav2/utils/muti_barplot.r "+expfile+" "+resultpath+"/"+pic_file_name
    os.system(cmd)
    if os.path.exists(resultpath+"/"+pic_file_name):
        return pic_file_name
    else:
        return "error"



def heatmap(degfilerfile,picpathname):
    df = pd.read_csv(degfilerfile, header=0, sep="\t")
    pre_heatmap = df[["control_mean", "treat_mean"]]
    pre_heatmap.index = list(df['Datebase_type'])
    sns.clustermap(pre_heatmap, cmap='RdYlGn_r', standard_scale=1)
    plt.savefig(picpathname + '.heatmap.png', format='png', dpi=100, bbox_inches='tight')
    plt.savefig(picpathname + '.heatmap.pdf', format='pdf', dpi=100, bbox_inches='tight')
    plt.close()