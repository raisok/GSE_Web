from django.shortcuts import render,redirect
from django.contrib.auth.decorators import login_required
from django.contrib import auth

from django.http import HttpResponse,JsonResponse
from django.core.paginator import Paginator, Page, EmptyPage, PageNotAnInteger
from django.views.decorators.csrf import csrf_exempt

import os
import json
from  datetime import datetime
import time
import shutil
from utils import fetch_gene,get_deg_exp,draw_bar,fetch_general_gene,draw_bar_v2
from .models import *
import pandas as pd
import logging
# Create your views here.

def index(request):
    return render(request, 'genesearch.html', {'title': "系统生物学平台"})

def getgenelist(genelist):
    genelist = genelist.replace("\r","")
    glist = genelist.split("\n")
    return glist

@csrf_exempt
def general(request):
    return render(request, 'genesearch.html', {'title': "系统生物学平台"})

def general2(request):
    return render(request,'gene_general.html')

def barplot(request):
    return render(request, 'barplot.html')

def heatmap(request):
    return render(request, 'heatmap.html')

def volmap(request):
    return render(request, 'volmap.html')

def upload(request):
    return render(request, 'upload_genes.html')

def single(request):
    return render(request, 'singlegene.html')

@csrf_exempt
def gsequery(request):
    if request.is_ajax():
        print(request.POST)
        disease,sepcies,tissue = request.POST['dstype'],request.POST['sptype'],request.POST['titype']
        print(disease,sepcies,tissue)
        try:
            gse_list = dataset[disease][sepcies][tissue]
        except KeyError:
            gsedic = {'gse': []}
        else:
            gsedic = {'gse': gse_list}
        print(gsedic)
        response = JsonResponse(gsedic)
        return response

@csrf_exempt
def general_genes(request):
    if request.is_ajax():
        print(request.POST)
        query_gene_list = request.POST.get('query')
        disease_type = request.POST.get('disease_type')
        species_type = request.POST.get('species_type')
        tissue_type = request.POST.get('tissue_type')
        try:
            gsex = dataset.get(disease_type).get(species_type).get(tissue_type)
        except AttributeError:
            return HttpResponse("error")
        #######
        if len(query_gene_list) == 0:
            return render(request, 'error01.html', {'errorlog': "您输入的基因列表为空"})
        elif query_gene_list.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
            genes = getgenelist(query_gene_list)
            genes = map(lambda x: x.upper(), genes)
            genes = tuple(genes)
        else:
            return render(request, 'error01.html', {'errorlog': "您输入的基因列表包含非法字符串，请检查你的基因列表"})
        if disease_type in dataset and len(genes) !=0:
            datapath,rpath,mediapath=createpath()
            gsepath = datapath + "/" + disease_type+"/"+species_type+"/"+tissue_type
            logging.warning("dataset path is: "+gsepath)
            degfilelist, expfilelist = fetch_general_gene.get_exp_list(gsepath)
            degfilename,expfilename = fetch_general_gene.get_result(genes, degfilelist, expfilelist, gsepath, rpath)
            downexp = rpath + "/" + expfilename
            downdiff = rpath + "/" + degfilename
            inputgenenum = len(genes)
            #Resultpath.objects.all().delete()
            #Resultpath.objects.create(resultpath=rpath, expfilepath=downexp, difffilepath=downdiff,
             #                         inputgenenum=inputgenenum,
             #                         picpath=rpath)
            #判断diff文件是否有输出
            fr = open(downdiff,'r')
            row=fr.readlines()[1:]
            fr.close()
            if os.path.exists(downdiff) and os.path.exists(downexp) and len(row):
                filedic = {'downexp': downexp, 'downdiff': downdiff}
                response = JsonResponse(filedic)
                return response
            else:
                return HttpResponse(downexp)

        else:
            return HttpResponse("badchose")
    else:
        return render(request, 'genesearch.html', {'title': "系统生物学平台"})

@csrf_exempt
def general_genes_table(request):
    if request.method == 'GET':
        resultfilepath=request.GET.get('name')
        dic = []
        print (resultfilepath+".filter.xls")
        fw = open(resultfilepath+".filter.xls",'w')
        fw.write("Database_type\tgene\tcontrol_mean\ttreat_mean\tlog2FoldChange\tpvalue\tpadj\n")
        fr = open(resultfilepath, 'r')
        for inline in fr.readlines()[1:]:
            row = inline.strip().split("\t")
            if len(row):
                dic.append({'Database_type': row[0],
                            'gene': row[1],
                            'control_mean': row[2],
                            'treat_mean': row[3],
                            'log2FoldChange': row[4],
                            'pvalue':row[5],
                            'padj': row[6]
                            })
                if (float(row[5]) <0.05) and(float(row[4]) > 0.58496 or float(row[4]) < -0.58496):
                    fw.write("\t".join(row)+"\n")
        fr.close()
        fw.close()
        if len(dic):
            print("all things be ok")
            return HttpResponse(json.dumps(dic), content_type="application/json")
        else:
            print("data is blank")
            dic=[{
                'Database_type': 'NA',
                'gene': 'NA',
                'control_mean': 'NA',
                'treat_mean':'NA',
                'log2FoldChange': 'NA',
                'pvalue':'NA',
                'padj':'NA'
            }]
            return HttpResponse(json.dumps(dic), content_type="application/json")

@csrf_exempt
def single_gene(request):
    if request.is_ajax():
        print(request.POST)
        disease_type = request.POST.get('disease_type')
        species_type = request.POST.get('species_type')
        tissue_type = request.POST.get('tissue_type')
        query_gene_symbol = request.POST.get('gene_symbol')
        try:
            gsex = dataset.get(disease_type).get(species_type).get(tissue_type)
        except AttributeError:
            return HttpResponse("error")
        if len(query_gene_symbol) == 0 or len(query_gene_symbol) >= 15:
            return HttpResponse("error")
        elif query_gene_symbol.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
            genes = getgenelist(query_gene_symbol)
            genes = map(lambda x: x.upper(), genes)
            genes = tuple(genes)
        else:
            return HttpResponse("error")
        if disease_type in dataset and len(genes) != 0 :
            datapath, rpath, mediapath = createpath()
            gsepath = datapath + "/" + disease_type+"/"+species_type+"/"+tissue_type
            logging.warning("dataset path is: "+gsepath)
            degfilelist, expfilelist = fetch_general_gene.get_exp_list(gsepath)
            degfilename,expfilename = fetch_general_gene.get_result(genes, degfilelist, expfilelist, gsepath, rpath)
            downexp = rpath + "/" + expfilename
            downdiff = rpath + "/" + degfilename
            logging.warning("exp file is: " + downexp)
            logging.warning("diff file is: " + downdiff)
            ###### 展示的图片是否需要用过滤后的 #####
            degfilter=rpath+'/'+degfilename+".filter.xls"
            df = pd.read_csv(rpath+'/'+degfilename, header=0, sep="\t")
            df = df[(df['log2FoldChange'] > 0.58496) | (df['log2FoldChange'] < -0.58496)]
            df = df[df['pvalue'] < 0.05 ]

            df.to_csv(degfilter, header=True, sep="\t", index=False)
            ######
            count = wcnum(rpath+'/'+degfilename)
            if count ==1:
                return HttpResponse("error")
            else:
                picfile = draw_bar_v2.gettablelist(genes, rpath+'/'+degfilename, mediapath)
                barplot_pic = draw_bar_v2.draw_exp_barplot(rpath+'/'+expfilename, mediapath)
                pdfpic = picfile + '.pdf'
                picps = {'pic': pdfpic,'downdiff':rpath+"/"+degfilename,'downexp':rpath+"/"+expfilename,'multplot':barplot_pic}
                response = JsonResponse(picps)
                return response
        else:
            return HttpResponse("error")

@csrf_exempt
def single_gene_table(request):
    if request.method == 'GET':
        resultfilepath=request.GET.get('name')
        dic = []
        print (resultfilepath+".filter.xls")
        fw = open(resultfilepath+".filter.xls",'w')
        fw.write("Database_type\tgene\tcontrol_mean\ttreat_mean\tlog2FoldChange\tpvalue\tpadj\n")
        fr = open(resultfilepath, 'r')
        for inline in fr.readlines()[1:]:
            row = inline.strip().split("\t")
            if len(row):
                dic.append({'Database_type': row[0],
                            'gene': row[1],
                            'control_mean': row[2],
                            'treat_mean': row[3],
                            'log2FoldChange': row[4],
                            'pvalue':row[5],
                            'padj': row[6]
                            })
                if (float(row[5]) <0.05) and(float(row[4]) > 0.58496 or float(row[4]) < -0.58496):
                    fw.write("\t".join(row)+"\n")
        fr.close()
        fw.close()
        if len(dic):
            print("all things be ok")
            return HttpResponse(json.dumps(dic), content_type="application/json")
        else:
            print("data is blank")
            dic=[{
                'Database_type': 'NA',
                'gene': 'NA',
                'control_mean': 'NA',
                'treat_mean':'NA',
                'log2FoldChange': 'NA',
                'pvalue':'NA',
                'padj':'NA'
            }]
            return HttpResponse(json.dumps(dic), content_type="application/json")

@csrf_exempt
def upload_genes(request):
    if request.is_ajax():
        file_obj = request.FILES.get('file')
        print(request.POST)
        querygenelist = request.POST.get('querygenelist')
        disease_type = request.POST.get('disease_type')
        species_type = request.POST.get('species_type')
        tissue_type = request.POST.get('tissue_type')
        query_gene_symbol = request.POST.get('gene_symbol')
        try:
            gsex = dataset.get(disease_type).get(species_type).get(tissue_type)
        except AttributeError:
            return HttpResponse("error")
        if file_obj:
            if os.path.exists('media/'+file_obj.name):
                os.remove('media/'+file_obj.name)
            f = open(os.path.join('media', file_obj.name), 'wb')
            for line in file_obj.chunks():
                f.write(line)
            f.close()
            genelistpath='media/'+file_obj.name
            return HttpResponse(genelistpath)
        elif querygenelist:
            genelist=[]
            with open(querygenelist,'r') as f:
                for line in f.readlines():
                    line=line.strip()
                    if line.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
                        genelist.append(line)
                    else:
                        return render(request, 'error01.html', {'errorlog': "您输入的基因列表有误，请检查你的基因列表"})
            genes = map(lambda x: x.upper(), genelist)
            genes = tuple(genes)
            logging.warning(genes)
            if disease_type in dataset and genes is not None:
                datapath,rpath,mediapath=createpath()
                gsepath = datapath + "/" + disease_type + "/" + species_type + "/" + tissue_type
                logging.warning("dataset path is: " + gsepath)
                degfilelist, expfilelist = get_deg_exp.get_exp_list(gsepath)
                degfilename,expfilename = fetch_general_gene.get_result(genes, degfilelist, expfilelist, gsepath, rpath)
                downexp = rpath + "/" + expfilename
                downdiff = rpath + "/" + degfilename
                ###### 展示的图片是否需要用过滤后的 #####
                degfilter = rpath + '/' + degfilename + ".filter.xls"
                df = pd.read_csv(rpath + '/' + degfilename, header=0, sep="\t")
                df = df[(df['log2FoldChange'] > 0.58496) | (df['log2FoldChange'] < -0.58496)]
                df = df[df['pvalue'] < 0.05]
                df.to_csv(degfilter, header=True, sep="\t", index=False)
                if os.path.exists(downdiff) and os.path.exists(downexp) and os.path.exists(degfilter):
                    filedic={'downexp':downexp,'downdiff':downdiff,'downfilterdiff':degfilter}
                    response = JsonResponse(filedic)
                    return response
                else:
                    return HttpResponse("False")
            else:
                return  HttpResponse("bad things happen")

def wcnum(file):
    count = 0
    thefile = open(file, 'rb')
    while 1:
        buffer = thefile.read(65536)
        if not buffer: break
        count += buffer.count(b'\n')  # 通过读取换行符计算
    return count

@csrf_exempt
def showtable(request):
    gse_direct = ["Liver_cancer", "Myocardial_Infarction", "Macrophage", "VIH", "Renal_fibrosis"]
    if request.method == 'POST':
        # print (request.POST)
        #获取从表单以POST方式提交过来的数据
        query_gene_list = request.POST.get('input_gene_list')
        disease_type = request.POST.get('disease_type')
        species_type = request.POST.get('species_type')
        tissue_type = request.POST.get('tissue_type')
        #######
        if len(query_gene_list) == 0:
            return render(request, 'error01.html', {'errorlog': "您输入的基因列表为空"})
        elif query_gene_list.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
            genes = getgenelist(query_gene_list)
            genes = map(lambda x: x.upper(), genes)
            genes = tuple(genes)
        else:
            return render(request, 'error01.html', {'errorlog': "您输入的基因列表包含非法字符串，请检查你的基因列表"})
        if disease_type in gse_direct and len(genes) !=0:
            datapath,rpath,mediapath=createpath()
            gsepath = datapath + "/" + disease_type
            degfilelist, expfilelist = fetch_gene.get_exp_list(gsepath)
            difftable,degfilename,expfilename = fetch_gene.get_result(genes, degfilelist, expfilelist, gsepath, rpath)
            # return HttpResponse("成功生成了你想要的文件")
            downexp = rpath + "/" + degfilename
            downdiff = rpath + "/" + expfilename
            inputgenenum = len(genes)
            Resultpath.objects.all().delete()
            Resultpath.objects.create(resultpath=rpath, expfilepath=downexp, difffilepath=downdiff,
                                      inputgenenum=inputgenenum,
                                      picpath=rpath)
            tablepre = DEGs.objects.all()
            resultdata = Resultpath.objects.all()

            #获取表头
            hlist = ["Database_type", "gene_id", "log2FoldChange", "pvalue", "padj"]
            tablenum = int(resultdata.values('inputgenenum')[0]['inputgenenum'])+1
            resultpath = resultdata.values('resultpath')[0]['resultpath']
            downexp = resultdata.values('expfilepath')[0]['expfilepath']
            downdiff = resultdata.values('difffilepath')[0]['difffilepath']
            picpath = resultdata.values('picpath')[0]['picpath']
            paginator = Paginator(tablepre, tablenum)
            page = request.GET.get('page')
            if page:  # 判断：获取当前页码的数据集，这样在模版就可以针对当前的数据集进行展示
                data_list = paginator.page(page).object_list
            else:
                data_list = paginator.page(1).object_list
            try:  # 实现分页对象，分别判断当页码存在/不存在的情况，返回当前页码对象
                page_object = paginator.page(page)
            except PageNotAnInteger:
                page_object = paginator.page(1)
            except EmptyPage:
                page_object = paginator.page(paginator.num_pages)
            return render(request, 'show_table.html', {
                'page_object': page_object,
                'data_list': data_list,
                'myList': hlist,
                'resultpath': resultpath,
                'expfile': downexp,
                'difffile': downdiff,
                'img': "output.pdf"
            })
        else:
            return HttpResponse("其他的数据库还未开放，目前只能查询Liver_cancer数据库." )
    else:
        # 获取表格数据
        tablepre = DEGs.objects.all()
        resultdata = Resultpath.objects.all()
        #获取表头
        hlist = ["Database_type", "gene_id", "log2FoldChange", "pvalue", "padj"]
        print (Resultpath.objects.values('inputgenenum')[0]['inputgenenum'])
        tablenum = int(resultdata.values('inputgenenum')[0]['inputgenenum'])+1
        resultpath = resultdata.values('resultpath')[0]['resultpath']
        downexp = resultdata.values('expfilepath')[0]['expfilepath']
        downdiff = resultdata.values('difffilepath')[0]['difffilepath']
        print (downdiff)
        print (downexp)
        picpath = resultdata.values('picpath')[0]['picpath']
        paginator = Paginator(tablepre, tablenum)
        page = request.GET.get('page')
        if page:  # 判断：获取当前页码的数据集，这样在模版就可以针对当前的数据集进行展示
            data_list = paginator.page(page).object_list
        else:
            data_list = paginator.page(1).object_list
        try:  # 实现分页对象，分别判断当页码存在/不存在的情况，返回当前页码对象
            page_object = paginator.page(page)
        except PageNotAnInteger:
            page_object = paginator.page(1)
        except EmptyPage:
            page_object = paginator.page(paginator.num_pages)
        # 分页功能
        return render(request, 'show_table.html', {
            'page_object': page_object,
            'data_list': data_list,
            'myList': hlist,
            'resultpath': resultpath,
            'expfile': downexp,
            'difffile': downdiff,
            'img': "output.pdf"
        })


@csrf_exempt
def detail2(request):
    gse_direct = ["Liver_cancer", "Myocardial_Infarction", "Macrophage", "VIH", "Renal_fibrosis"]
    if request.is_ajax():
        file_obj = request.FILES.get('file')
        print(request.POST)
        querygenelist = request.POST.get('querygenelist')
        disease_type = request.POST.get('disease_type')
        species_type = request.POST.get('species_type')
        tissue_type = request.POST.get('tissue_type')
        query_gene_symbol = request.POST.get('gene_symbol')
        if file_obj:
            if os.path.exists('media/'+file_obj.name):
                os.remove('media/'+file_obj.name)
            f = open(os.path.join('media', file_obj.name), 'wb')
            for line in file_obj.chunks():
                f.write(line)
            f.close()
            genelistpath='media/'+file_obj.name
            return HttpResponse(genelistpath)
        elif querygenelist:
            genelist=[]
            with open(querygenelist,'r') as f:
                for line in f.readlines():
                    line=line.strip()
                    if line.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
                        genelist.append(line)
                    else:
                        return render(request, 'error01.html', {'errorlog': "您输入的基因列表有误，请检查你的基因列表"})
            genes = map(lambda x: x.upper(), genelist)
            genes = tuple(genes)
            gse_direct = ["Liver_cancer", "Myocardial_Infarction", "Macrophage", "VIH", "Renal_fibrosis"]
            print(genes)
            if disease_type in gse_direct and genes is not None:
                datapath,rpath,mediapath=createpath()
                gsepath = datapath + "/" + disease_type
                degfilelist, expfilelist = get_deg_exp.get_exp_list(gsepath)
                degfilename,expfilename = get_deg_exp.get_result(genes, degfilelist, expfilelist, gsepath, rpath)
                downexp = rpath + "/" + expfilename
                downdiff = rpath + "/" + degfilename
                if os.path.exists(downdiff) and os.path.exists(downexp):
                    filedic={'downexp':downexp,'downdiff':downdiff}
                    response = JsonResponse(filedic)
                    return response
                else:
                    return HttpResponse("False")
            else:
                return  HttpResponse("bad things happen")
        elif query_gene_symbol:
            if len(query_gene_symbol) == 0 or len(query_gene_symbol) >= 15:
                return render(request, 'error01.html', {'errorlog': "您输入的基因名称有误"})
            elif query_gene_symbol.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
                genes = getgenelist(query_gene_symbol)
                genes = map(lambda x: x.upper(), genes)
                genes = tuple(genes)
            else:
                return render(request, 'error01.html', {'errorlog': "您输入的基因列表包含非法字符串，请检查你的基因列表"})
            if disease_type in gse_direct and len(genes) != 0:
                datapath, rpath,mediapath = createpath()
                gsepath = datapath + "/" + disease_type
                degfilelist, expfilelist = fetch_gene.get_exp_list(gsepath)
                picfile = draw_bar.gettablelist(genes, degfilelist, gsepath, mediapath)
                pdfpic = picfile + '.pdf'
                picps={'pic':pdfpic}
                response = JsonResponse(picps)
                return response
            else:
                return HttpResponse("error")
        else:
            liver = {
                'Liver_cancer': ["GSE104310", "GSE112221", "GSE115193", "GSE121248", "GSE14520_U133A", "GSE14520_U133A2.0",
                                 "GSE1946", "GSE364", "GSE45050", "GSE5093", "GSE57957", "GSE60502", "GSE62232", "GSE84598",
                                 "GSE99010"], }
            response = JsonResponse(liver)
            return response
    else:
        return render(request, 'genesearch.html', {'title': "系统生物学平台"})

def createpath():
    #datapath="/media/sdb/data_mining"
    #datapath = "E:/IMA/web/imasite/data/database"
    #datapath = "/mnt/e/IMA/web/imasite/data/database"
    datapath="/home/user/data_mining"
    now = datetime.now()
    now = now.strftime('%Y%m%d%H')
    #rpath = "E:/IMA/web/imav2/result/" + now
    #rpath = "/mnt/e/IMA/web/imav2/result/" + now
    #rpath = "/media/sdc/yueyao/930/Web/test2/imav2/result/"+now
    rpath = "/home/user/Web/imav2/result/"+now
    #rpath = "/media/sdc/yueyao/Web/imav2/result/"+now
    mediapath = os.path.dirname(rpath)+"/pic"
    # if os.path.exists(mediapath):
    #     shutil.rmtree(mediapath)
    #     os.makedirs(mediapath)
    # else:
    #     os.makedirs(mediapath)
    # if os.path.exists(rpath):
    #     shutil.rmtree(rpath)
    #     os.makedirs(rpath)
    # else:
    #     os.makedirs(rpath)
    return datapath,rpath,mediapath

def readFile(fn, buf_size=262144):  # 大文件下载，设定缓存大小
    f = open(fn, "rb")
    while True:  # 循环读取
        c = f.read(buf_size)
        if c:
            yield c
        else:
            break

def download_file(request):
    if request.method == 'GET':
        # for v in request.GET.keys():
        filepath_ = request.GET.get('id')
        filename_ = os.path.basename(filepath_)
        response = HttpResponse(readFile(filepath_), content_type='APPLICATION/OCTET-STREAM')
        response['Content-Disposition'] = 'attachment; filename=' + filename_
        response['Content-Length'] = os.path.getsize(filepath_)  # 传输给客户端的文件大小
        return response
    else:
        HttpResponse("提交数据的方式存在问题？")


dataset={
    "Glucose_metabolism":{
        "Mus_musculus":{
            "Liver_Glucagon":[
                "GSE113526",
                "GSE110673"
            ],
            "Adipose_Fasted":[
                "GSE46495_adipose"
            ],
            "Liver_dbdb":[
                "GSE30140",
                "GSE59930"
            ],
            "Liver_Fasted":[
                "GSE46495_liver",
                "GSE73299",
                "GSE100313",
                "GSE10653",
                "GSE51712"
            ],
            "Liver_Insulin":[
                "GSE45694"
            ],
            "Muscle_Fasted":[
                "GSE46495_muscle"
            ]
        },
        "Homo_sapiens":{
            "Adipose_Insulin":[
                "GSE15773_subcutaneous",
                "GSE15773_omental"
            ],
            "Liver_Glucagon":[
                "GSE68144"
            ],
            "Muscle_Insulin":[
                "GSE7146_U133A"
            ]
        }
    },
    "Liver_cancer":{
        "Homo_sapiens":{
            "Liver":[
                "GSE84598",
                "GSE121248",
                "GSE57957",
                "GSE60502",
                "GSE112221",
                "GSE115193",
                "GSE62232",
                "GSE14520_U133A",
                "GSE45050",
                "GSE99010",
                "GSE104310",
                "GSE14520_U133A2.0"
            ]
        }
    },
    "Macrophage":{
        "Homo_sapiens":{
            "Macrophage":[
                "GSE117040",
                "GSE5099B",
                "GSE36537",
                "GSE86298",
                "GSE55536",
                "GSE5099A",
                "GSE32164",
                "GSE57614",
                "GSE18686"
            ]
        },
        "Mus_musculus":{
            "Macrophage":[
                "GSE84517",
                "GSE53053",
                "GSE51466",
                "GSE81922",
                "GSE69607"
            ]
        }
    },
    "Myocardial_Infarction":{
        "Mus_musculus":{
            "Heart":[
                "GSE123092",
                "GSE46395",
                "GSE76387",
                "GSE54132",
                "GSE19322",
                "GSE71906",
                "GSE62973"
            ]
        },
        "Rattus_norvegicus":{
            "Heart":[
                "GSE52313"
            ]
        }
    },
    "NASH":{
        "Homo_sapiens":{
            "Liver":[
                "GSE105127",
                "GSE48452",
                "GSE126848",
                "GSE37031",
                "GSE115193",
                "GSE59045",
                "GSE61260",
                "GSE66676",
                "GSE63067"
            ]
        },
        "Mus_musculus":{
            "Liver_HFD":[
                "GSE43106",
                "GSE51885",
                "GSE53834",
                "GSE83940",
                "GSE95428",
                "GSE9484",
                "GSE109345",
                "GSE32095",
                "GSE24031",
                "GSE119953",
                "GSE57425",
                "GSE95283",
                "GSE93132",
                "GSE70119",
                "GSE94754",
                "GSE97272"
            ],
            "Liver_HFHC":[
                "GSE57290",
                "GSE93819",
                "GSE51432"
            ]
        }
    },
    "VIH":{
        "Homo_sapiens":{
            "Endothelial_cells":[
                "GSE96962"
            ]
        },
        "Mus_musculus":{
            "Vein&Artery":[
                "GSE119549"
            ]
        }
    },
    "Myocardial_hypertrophy":{
        "Homo_sapiens":{
            "Heart":[
                "GSE3585"
            ]
        },
        "Mus_musculus":{
            "Heart":[
                "GSE1621",
                "GSE56348",
                "GSE76"
            ]
        }
    }
}

def Return_Species_Data(request):
    diseasetype = request.GET['Disease']
    species_list = []
    for species in dataset[diseasetype]:
        species_list.append(species)
    return HttpResponse(json.dumps(species_list))

def Return_Tissue_Data(request):
    disease,sepcies = request.GET['Disease'],request.GET['Species']
    tissue_list = []
    for tissue in dataset[disease][sepcies]:
        tissue_list.append(tissue)
    print(tissue_list)
    return HttpResponse(json.dumps(tissue_list))


