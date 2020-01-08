#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     viewstest.py
   Description :
   Author :       yueyao
   date：          2019/8/30
-------------------------------------------------
   Change Activity:
                   2019/8/30:
-------------------------------------------------
"""
__author__ = 'yueyao'

#测试是否正常启动django
def hello(request):
    return HttpResponse("hello my django...<br>这是一个测试")

def detail(request):
    if request.method == "POST":
        disease_type = request.POST.get('disease_type')
        species_type = request.POST.get('species_type')
        tissue_type = request.POST.get('tissue_type')
        clear_ = request.POST.get('clear_showgene_list')
        print (request.POST)
        show_=request.POST.get('show_submit')
        # print (show_)
        genelist=request.POST.get('input_gene_list')
        # print (genelist)
        query_=request.POST.get('query_submit')
        if show_ and genelist:
            glist=getgenelist(genelist)
            for gg in glist:
                print(gg)
            ggl="\r\n".join(glist)
            print (ggl)
            if not ggl.isalpha():
                return HttpResponse("你的基因列表可能有误！")
            else:
                return render(request,'detail.html',{'showgenelist':ggl,'title':"系统生物学平台"})
        elif query_ and genelist:
            glist = getgenelist(genelist)
            #这里需要根据传的值在后台去不同的GSE库里面寻找相应的基因
            fpkm="E:/IMA/web/imasite/ima/all.genes.fpkm.xls"
            hlist,context = getfromfile(glist,fpkm)
            return render(request, 'show_table.html', {'contexts': context, 'myList': hlist})
        elif clear_:
            return render(request, 'detail.html', {'title': "系统生物学平台"})
        else:
            return render(request, 'detail.html', {'preshow': "请输入需要查询的基因列表",'title':"系统生物学平台"})
    else:
        return render(request,'detail.html',{'preshow':"请输入需要查询的基因列表",'title':"系统生物学平台",'showgenelist':"这里显示需要查询的基因列表"})

def readFile(fn, buf_size=262144):  # 大文件下载，设定缓存大小
    f = open(fn, "rb")
    while True:  # 循环读取
        c = f.read(buf_size)
        if c:
            yield c
        else:
            break

def showtable(request):
    hlist, context=readfpkm("E:/IMA/web/imasite/ima/all.genes.fpkm.xls")
    return render(request,'show_table.html',{'contexts':context,'myList':hlist})

def result(request):
    if request.method == "POST":
        # print (request.POST)
        show2=request.POST.get('submit2')
        genelist=request.POST.get('gene_list')
        chaxun=request.POST.get('submit')
        if show2 and genelist:
            glist=getgenelist(genelist)
            for gg in glist:
                print(gg)
            ggl="\r\n".join(glist)
            return render(request,'result1.html',{'showgenelist':ggl})
        elif chaxun and genelist:
            glist = getgenelist(genelist)
            #这里需要根据传的值在后台去不同的GSE库里面寻找相应的基因
            fpkm="E:/IMA/web/imasite/ima/all.genes.fpkm.xls"
            hlist,context = getfromfile(glist,fpkm)
            return render(request, 'show_table.html', {'contexts': context, 'myList': hlist})
        else:
            return render(request, 'result1.html', {'preshow': "请输入需要查询的基因列表"})
    else:
        return render(request,'result1.html',{'preshow':"请输入需要查询的基因列表"})

def getgenelist(genelist):
    genelist = genelist.replace("\r","")
    glist = genelist.split("\n")
    return glist

# import pandas as pd
def getfromfile(glist,gfpkm):
    with open(gfpkm,'r') as f:
        header=f.readline()
        hlist=header.strip().split("\t")
        context=[]
        for line in f.readlines():
            genel=line.strip().split("\t")
            if genel[0] in glist:
                context.append(genel)
    return hlist,context

def readfpkm(fpkm):
    with open(fpkm,'r') as f:
        header=f.readline()
        hlist=header.strip().split("\t")
        context=[]
        for con in f.readlines()[0:100]:
            con=con.strip().split("\t")
            context.append(con)
        return hlist,context