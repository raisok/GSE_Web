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
from utils import fetch_gene,get_deg_exp,draw_bar,fetch_general_gene,draw_bar_v2,venn_plot
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

def venn2(request):
    return render(request,'venn2.html')

def venn(request):
    return render(request,'venn.html')


@csrf_exempt
def venn_result(request):
    if request.is_ajax():
        now = datetime.now()
        datapath, rpath, mediapath = createpath()
        outname = now.strftime('%Y%m%d%H') + str(v_code())
        output = mediapath + "/" + outname
        print(request.POST)
        list1=request.POST.get('list1')
        list2 = request.POST.get('list2')
        list3 = request.POST.get('list3')
        list4 = request.POST.get('list4')
        list5 = request.POST.get('list5')
        gene_list1 = request.POST.get('gene_list1')
        gene_list2 = request.POST.get('gene_list2')
        gene_list3 = request.POST.get('gene_list3')
        gene_list4 = request.POST.get('gene_list4')
        gene_list5 = request.POST.get('gene_list5')

        all_list = [list1,list2,list3,list4,list5]
        all_gene_list = [gene_list1,gene_list2,gene_list3,gene_list4,gene_list5]
        venn_index_list = [i for i,j in enumerate(all_gene_list) if not j == '']
        print(venn_index_list)
        ch_all_gene_list = [getgenelist(n) for n in all_gene_list]
        if len(venn_index_list) == 0:
            return "bad"
        elif len(venn_index_list) == 1:
            return "bad"
        elif len(venn_index_list) == 2:
            labels = venn_plot.get_labels([ch_all_gene_list[venn_index_list[0]],ch_all_gene_list[venn_index_list[1]]],
                                     fill=['number', 'logic'])
            fig, ax = venn_plot.venn2(labels, names=[all_list[venn_index_list[0]],all_list[venn_index_list[1]]])
            fig.savefig(output + '.pdf', format='pdf', dpi=100, bbox_inches='tight')
            venn_dic = {'venn_pic': outname + '.pdf'}
            response = JsonResponse(venn_dic)
            return response
        elif len(venn_index_list) == 3:
            labels = venn_plot.get_labels([ch_all_gene_list[venn_index_list[0]],ch_all_gene_list[venn_index_list[1]],
            ch_all_gene_list[venn_index_list[2]]],
                                     fill=['number', 'logic'])
            fig, ax = venn_plot.venn3(labels, names=[all_list[venn_index_list[0]],all_list[venn_index_list[1]],all_list[venn_index_list[2]]])

            fig.savefig(output + '.pdf', format='pdf', dpi=100, bbox_inches='tight')
            venn_dic = {'venn_pic': outname + '.pdf'}
            response = JsonResponse(venn_dic)
            return response
        elif len(venn_index_list) == 4:
            labels = venn_plot.get_labels([ch_all_gene_list[venn_index_list[0]],
                                      ch_all_gene_list[venn_index_list[1]],
                                      ch_all_gene_list[venn_index_list[2]],
                                      ch_all_gene_list[venn_index_list[3]]
                                      ],
                                     fill=['number', 'logic'])
            fig, ax = venn_plot.venn4(labels, names=[all_list[venn_index_list[0]],
                                                all_list[venn_index_list[1]],
                                                all_list[venn_index_list[2]],
                                                all_list[venn_index_list[3]]
                                                ])
            fig.savefig(output + '.pdf', format='pdf', dpi=100, bbox_inches='tight')
            venn_dic = {'venn_pic': outname + '.pdf'}
            response = JsonResponse(venn_dic)
            return response
        elif len(venn_index_list) ==5:
            labels = venn_plot.get_labels([ch_all_gene_list[venn_index_list[0]],
                                      ch_all_gene_list[venn_index_list[1]],
                                      ch_all_gene_list[venn_index_list[2]],
                                      ch_all_gene_list[venn_index_list[3]],
                                      ch_all_gene_list[venn_index_list[4]],
                                      ],
                                     fill=['number', 'logic'])
            fig, ax = venn_plot.venn5(labels, names=[all_list[venn_index_list[0]],
                                                all_list[venn_index_list[1]],
                                                all_list[venn_index_list[2]],
                                                all_list[venn_index_list[3]],
                                                all_list[venn_index_list[4]],
                                                ])

            fig.savefig(output + '.pdf', format='pdf', dpi=100, bbox_inches='tight')
            venn_dic = {'venn_pic': outname + '.pdf'}
            response = JsonResponse(venn_dic)
            return response
        else:
            return "bad"


def draw_venn(request):
    return render(request,'venn.html')

from sklearn.decomposition import  PCA
from sklearn.preprocessing import  StandardScaler
import  seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # 空间三维画图
@csrf_exempt
def draw_pca(request):
    now = datetime.now()
    datapath, rpath, mediapath = createpath()
    now = now.strftime('%Y%m%d%H') + str(v_code())
    print(request.POST)
    fpkm_file = request.POST.get('fpkm')
    components2=request.POST.get('components2')
    components3=request.POST.get('components3')
    print(components2)
    print(components3)
    if fpkm_file == "draw":
        return HttpResponse("bad")
    else:
        outname = now
        data = pd.read_csv(fpkm_file, sep="\t",header=0)
        index_name = data.columns
        geneid=index_name[0]
        data2 = data.drop(geneid, axis=1)
        data3 = data2.T
        if int(components2) == 2 and int(components3) == 0:
            print("runing PCA 2 compnonent...")
            pca = PCA(n_components=2)
            pca.fit(data3)
            sample_pca_data = pca.transform(data3)
            df = pd.DataFrame(sample_pca_data, columns=["PC1", "PC2"], index=index_name[1:])
            # 如下命令需要使用seaborn 0.10版本
            sns.scatterplot(x=df['PC1'], y=df['PC2'], hue=df.index)
            plt.savefig(mediapath+"/"+outname + '.png', format='png', dpi=50, bbox_inches='tight')
            plt.savefig(mediapath+"/"+outname + '.pdf', format='pdf', dpi=100, bbox_inches='tight')
            plt.close()
            picps = {'pic': outname+".pdf"}
            response = JsonResponse(picps)
            return response
        elif int(components3) == 3 and int(components2) == 0:
            print("runing PCA 3 compnonent...")
            pca = PCA(n_components=3)
            pca.fit(data3)
            sample_pca_data = pca.transform(data3)
            dfa = pd.DataFrame(sample_pca_data, columns=["PC1", "PC2", "PC3"], index=index_name[1:])
            col_list = sns.hls_palette(n_colors=len(dfa.index), l=0.7, s=0.7)
            fig = plt.figure(figsize=(12, 9))
            ax = Axes3D(fig)
            for i in range(len(dfa.index)):
                x = dfa.iloc[i, :]['PC1']
                y = dfa.iloc[i, :]['PC2']
                z = dfa.iloc[i, :]['PC3']
                ax.scatter(x, y, z, c=col_list[i], label=dfa.index[i])
            # 绘制图例
            ax.legend(loc='best')
            # 添加坐标轴(顺序是Z, Y, X)
            ax.set_zlabel('PC3', fontdict={'size': 15, 'color': 'red'})
            ax.set_ylabel('PC2', fontdict={'size': 15, 'color': 'red'})
            ax.set_xlabel('PC1', fontdict={'size': 15, 'color': 'red'})
            plt.savefig(mediapath + "/" + outname + '.png', format='png', dpi=50, bbox_inches='tight')
            plt.savefig(mediapath + "/" + outname + '.pdf', format='pdf', dpi=100, bbox_inches='tight')
            plt.close()
            picps = {'pic': outname+".pdf"}
            response = JsonResponse(picps)
            return response
        else:
            return HttpResponse("bad")


def goenrich(request):
    return render(request,'goenrich.html')

def keggenrich(request):
    return render(request,'keggenrich.html')


def pca(request):
    return render(request,'pca.html')

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
            #                          inputgenenum=inputgenenum,
            #                          picpath=rpath)
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
#        allresult = Resultpath.objects.all()
#        for i in allresult:
#            resultfilepath =i.difffilepath
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
            # return render(request, 'error01.html', {'errorlog': "您输入的基因名称有误"})
        elif query_gene_symbol.replace("\r", "").replace("\n", "").encode('UTF-8').isalnum():
            genes = getgenelist(query_gene_symbol)
            genes = map(lambda x: x.upper(), genes)
            genes = tuple(genes)
        else:
            return HttpResponse("error")
            # return render(request, 'error01.html', {'errorlog': "您输入的基因列表包含非法字符串，请检查你的基因列表"})
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
            inputgenenum = len(genes)
#            Resultpath.objects.all().delete()
#            Resultpath.objects.create(resultpath=rpath, expfilepath=downexp, difffilepath=downdiff,
#                                      inputgenenum=inputgenenum,
#                                      picpath=rpath)
            ###### 展示的图片是否需要用过滤后的 #####
            degfilter=rpath+'/'+degfilename+".filter.xls"
            df = pd.read_csv(rpath+'/'+degfilename, header=0, sep="\t")
            df = df[(df['log2FoldChange'] > 0.58496) | (df['log2FoldChange'] < -0.58496)]
            df = df[df['pvalue'] < 0.05 ]
            #df = df[df['padj'] < 0.05]
            df.to_csv(degfilter, header=True, sep="\t", index=False)
            ######
            count = wcnum(rpath+'/'+degfilename)
            if count ==1:
                return HttpResponse("error")
            else:
                picfile = draw_bar_v2.gettablelist(genes, rpath+'/'+degfilename, mediapath)
                #picfile = "xxx"
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
 #       allresult = Resultpath.objects.all()
 #       for i in allresult:
 #           resultfilepath =i.difffilepath
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
            if os.path.exists('/home/user/Web/imav2/media/'+file_obj.name):
                os.remove('/home/user/Web/imav2/media/'+file_obj.name)
            f = open(os.path.join('/home/user/Web/imav2/media/', file_obj.name), 'wb')
            for line in file_obj.chunks():
                f.write(line)
            f.close()
            genelistpath='/home/user/Web/imav2/media/'+file_obj.name
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
                #df = df[df['padj'] < 0.05]
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

            # print (resultdata)
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
def upload_fpkm(request):
    if request.is_ajax():
        file_obj = request.FILES.get('file')
        if file_obj:
            if os.path.exists('/home/user/Web/imav2/media/'+file_obj.name):
                os.remove('/home/user/Web/imav2/media/'+file_obj.name)
            f = open(os.path.join('/home/user/Web/imav2/media', file_obj.name), 'wb')
            for line in file_obj.chunks():
                f.write(line)
            f.close()
            genelistpath='/home/user/Web/imav2/media/'+file_obj.name
            return HttpResponse(genelistpath)

import random
from datetime import datetime
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

@csrf_exempt
def draw_heatmap(request):
    print(request.POST)
    datapath, rpath, mediapath = createpath()
    fpkm_file = request.POST.get('fpkm')
    row_cluster=request.POST.get('row_cluster')
    col_cluster=request.POST.get('col_cluster')
    now = datetime.now()
    now=now.strftime('%Y%m%d%H')+str(v_code())
    if fpkm_file == "draw":
        return HttpResponse("bad")
    else:
        outname = now+".pdf"
        outpdf = mediapath+"/"+outname
        heatmap_script='''
args  <- commandArgs(TRUE)
library(pheatmap)
input <- "{0}"

myfunction <- function(rt,IdColName){{
  rt1 <- rt
  idColNum <- grep(IdColName,colnames(rt))
  mat <- as.matrix(rt[,-idColNum])
  rownames(mat) <- rt[,idColNum]
  rt1$mean <- apply(mat, 1, sum)
  aa <- rt1[order(rt1[,IdColName],rt1$mean,decreasing = T),]
  aa <- aa[!duplicated(aa[,IdColName]),]
  aa <- subset(aa,select=-c(mean))
  return(aa)
}}

gene = read.table(input,header = T,sep = "\\t",check.names=F)

gene2 <- myfunction(gene,names(gene)[1])

gene = gene2[,-1]
row.names(gene) = gene2[,1]



data_plot<-t(scale(t(gene)))
data_plot[data_plot>3]<-3
data_plot[data_plot<(-3)]<-(-3)

pheatmap ( data_plot, cluster_rows ={1}, cluster_cols ={2}, show_rownames = F,
           border_color = "white",
           color =colorRampPalette(c("green", "black", "red"))(100),
           filename = "{3}"
)
'''.format(
            fpkm_file,row_cluster,col_cluster,outpdf
        )
        with open("/home/user/Web/imav2/media"+"/"+now+".r",'w') as f:
            f.write(heatmap_script)
        cmd = "cat /home/user/Web/imav2/media/"+now+".r"+" | R --vanilla"
        os.system(cmd)
        os.remove("/home/user/Web/imav2/media/"+now+".r")
        picps = {'pic': outname}
        response = JsonResponse(picps)
        return response

@csrf_exempt
def do_go_enrich(request):
    print(request.POST)
    species=request.POST.get('species_type')
    genelist = request.POST.get('query')
    genelist = genelist.strip()
    genelist = genelist.replace("\r","")
    genelist = genelist.replace("\n","|")
    dbs = {
        'mmu': 'org.Mm.eg.db',
        'hsa': 'org.Hs.eg.db',
        'rno': 'org.Rn.eg.db',
        'mcf': 'org.Mmu.eg.db',
        'ssc': 'org.Ss.eg.db'
    }
    if species in dbs and len(genelist) != "":
        now = datetime.now()
        datapath, rpath, mediapath = createpath()
        outname=now.strftime('%Y%m%d%H')+str(v_code())
        output = mediapath + "/" + outname
        GO_enrichment_R = '''
    library(clusterProfiler)
    library({Orgdb})
    library(ggplot2)
    library(RColorBrewer)
    
    
    ### change GO Term name length
    shorten_names2 <- function(x, n_word=4, n_char=40){{
      if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
      {{
        if (nchar(x) > 40) x <- substr(x, 1, 40)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                         collapse=" "), "...", sep="")
        return(x)
      }} 
      else
      {{
        return(x)
      }}
    }}
    
    shorten_names <- function(x, n_word=4, n_char=40){{
      shortname = x["Description"]
      Go_id = x["ID"]
      if (length(strsplit(shortname, " ")[[1]]) > n_word || (nchar(shortname) > n_char))
      {{
        if (nchar(shortname) > n_char) shortname <- substr(shortname, 1, n_char)
        y <- paste(paste(strsplit(shortname, " ")[[1]][1:min(length(strsplit(shortname," ")[[1]]), n_word)],
                         collapse=" "), "...", sep="")
        y <- paste(Go_id,y,sep=" / ")
        return(y)
      }}
      else
      {{
        return (shortname)
      }}
    }}
    ### add by yueyao at 20191121
    
    ### 
    set_width <- function(para1){{
      #para<- para1[order(para1$p.adjust),][1:15,]
      Des <- para1$pretty_varname
      max1<- max(nchar(c(as.character(Des))))
      width<- max1*0.09
      if(nrow(para1)>10) {{
       height <- 8*nrow(para1)/60
      }} else {{
        height <- 1.6
    
      }}
      return(c(width,height))
    }}
    ### add by hxl at 20191121
    
    num<- function(data,n){{
      if(nrow(data)>n) {{
        data <- data[1:n,]
      }} else {{
        data <- data
      }}
      return(data)
    }}
    
    
    ###
    
    draw_bar_plot <- function(goinput){{
      
      bar <- ggplot(goinput,aes(y=Count,x=reorder(pretty_varname,-log10(pvalue)))) +
        geom_bar(stat="identity", width=0.7,aes(fill=-log10(pvalue))) + coord_flip() + 
        scale_fill_gradientn(colours=c(rev(brewer.pal(9, "RdYlBu")))) +
        labs(color=expression(-log10(pvalue)),x="GO term",y="Gene Number",title="Most Top 15 enrich GO Term")+
        theme_bw() +
        theme(
          title =element_text(size=6,face="plain",colour="black"),
          axis.title = element_text(size=6,face="plain",colour="black"),
          panel.grid=element_blank(),
          
          axis.line = element_line(size=0.25,colour = "black"),
          axis.text = element_text(color = "black",size = 6,face="plain"),
          legend.key.size=unit(0.3,'cm'),
          legend.title=element_text(size=6,colour = "black",face="plain"),
          legend.text = element_text(colour = 'black',face = 'plain',size=6)
          )
        return (bar)
    }}
    
    ###
    
    DEGs <- strsplit('{input_str}',"\\\\|")[[1]]
    
    ont <- c("MF","BP","CC")
    type <- c("molecular_function","biological_process","cellular_component")
    
    ego_MF <- enrichGO(gene=as.vector(DEGs),OrgDb={Orgdb},ont=ont[1],pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=1,keyType='SYMBOL')
    write.table(ego_MF,file=paste("{output_prefix}",".",type[1],"_enrichment.xls",sep=""),sep="\\t",quote=F,row.names=F)
    
    if (! nrow(ego_MF) == 0){{
    ego_MF<- ego_MF[order(ego_MF$pvalue),]
    ego_MF <- num(ego_MF,15)
    
    ego_MF$pretty_varname = apply(ego_MF,1,shorten_names)
    w <- set_width(ego_MF)
    
    ego_MF<- draw_bar_plot(ego_MF)
    
    ggsave(file=paste("{output_prefix}",".",type[1],"_enrichment.pdf",sep=""),ego_MF,width=w[1],height=w[2])
    ggsave(file=paste("{output_prefix}",".",type[1],"_enrichment.png",sep=""),ego_MF,width=w[1],height=w[2])
    
    svg_file=paste("{output_prefix}",".",type[1],"_enrichment.svg",sep="")
    svg(svg_file,width=w[1],height=w[2])
    plot(ego_MF)
    dev.off()
    
    }}
    
    
    ego_BP <- enrichGO(gene=as.vector(DEGs),OrgDb={Orgdb},ont=ont[2],pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=1,keyType='SYMBOL')
    write.table(ego_BP,file=paste("{output_prefix}",".",type[2],"_enrichment.xls",sep=""),sep="\\t",quote=F,row.names=F)
    
    if (! nrow(ego_BP) == 0){{
    ego_BP<- ego_BP[order(ego_BP$pvalue),]
    ego_BP <- num(ego_BP,15)
    
    
    ego_BP$pretty_varname = apply(ego_BP,1,shorten_names)
    w <- set_width(ego_BP)
    ego_BP<-draw_bar_plot(ego_BP)
    
    ggsave(file=paste("{output_prefix}",".",type[2],"_enrichment.pdf",sep=""),ego_BP,width=w[1],height=w[2])
    ggsave(file=paste("{output_prefix}",".",type[2],"_enrichment.png",sep=""),ego_BP,width=w[1],height=w[2])
    
    
    svg_file=paste("{output_prefix}",".",type[2],"_enrichment.svg",sep="")
    svg(svg_file,width=w[1],height=w[2])
    plot(ego_BP)
    dev.off()
    
    }}
    
    
    ego_CC <- enrichGO(gene=as.vector(DEGs),OrgDb={Orgdb},ont=ont[3],pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=1,keyType='SYMBOL')
    write.table(ego_CC,file=paste("{output_prefix}",".",type[3],"_enrichment.xls",sep=""),sep="\\t",quote=F,row.names=F)
    
    if (! nrow(ego_CC) == 0){{
    ego_CC<- ego_CC[order(ego_CC$pvalue),]
    ego_CC <- num(ego_CC,15)
    
    
    ego_CC$pretty_varname = apply(ego_CC,1,shorten_names)
    w <- set_width(ego_CC)
    ego_CC<-draw_bar_plot(ego_CC)
    
    ggsave(file=paste("{output_prefix}",".",type[3],"_enrichment.pdf",sep=""),ego_CC,width=w[1],height=w[2])
    ggsave(file=paste("{output_prefix}",".",type[3],"_enrichment.png",sep=""),ego_CC,width=w[1],height=w[2])
    
    svg_file=paste("{output_prefix}",".",type[3],"_enrichment.svg",sep="")
    svg(svg_file,width=w[1],height=w[2])
    plot(ego_CC)
    dev.off()
    
    }}
    
    '''.format(
            Orgdb=dbs.get(species),
            output_prefix=output,
            input_str=genelist
        )
        with open("/home/user/Web/imav2/media"+"/"+outname+'.R','w') as f:
            f.write(GO_enrichment_R)
        cmd="cat /home/user/Web/imav2/media/"+outname+'.R'+"|  R --vanilla"
        os.system(cmd)
        os.remove("/home/user/Web/imav2/media"+"/"+now+'.R')
        result_dic={'CC_enrich':output+".cellular_component_enrichment.xls",
                    "MF_enrich":output+".molecular_function_enrichment.xls",
                    "BP_enrich":output+".biological_process_enrichment.xls",
                    "CC_pic":outname+".cellular_component_enrichment.pdf",
                    "MF_pic":outname+".molecular_function_enrichment.pdf",
                    "BP_pic":outname+".biological_process_enrichment.pdf"
                    }
        response = JsonResponse(result_dic)
        return response
    else:
        return HttpResponse("bad things happen")

@csrf_exempt
def do_kegg_enrich(request):
    print(request.POST)
    species=request.POST.get('species_type')
    genelist = request.POST.get('query')
    genelist = genelist.strip()
    genelist = genelist.replace("\r","")
    genelist = genelist.replace("\n","|")
    dbs = {
        'mmu': '/home/user/Web/imav2/utils/mmu_kegg.ko_20200407',
        'hsa': '/home/user/Web/imav2/utils/hsa_kegg.ko_20200407',
        'rno': '/home/user/Web/imav2/utils/rno_kegg.ko_20200407',
    }
    if species in dbs and len(genelist) != "":
        now = datetime.now()
        datapath, rpath, mediapath = createpath()
        outname = now.strftime('%Y%m%d%H') + str(v_code())
        output = mediapath + "/" + outname
        KEGG_enrichment_R = '''
#读入基因列表
degfile <- strsplit('{input_str}',"\\\\|")[[1]]
#根据选择的物种输入相应的ko文件
species_ko_file <- "{species_type}"
species <- "{species}"
species_ko <- "{species}"
#输出结果，包含路径
output <- "{output_prefix}"


del_human = T
gmtout = F
outkopath = F
p = 0.05
fdr = 1
zs = 0

library(hash)
library(ggplot2)
library(RColorBrewer)

kegg_ko <- read.table(species_ko_file,header=F,sep="\\t",fill=T,quote="",stringsAsFactors=F)


### 3.基因对应通路等级表
ko <- data.frame(kegg_ko[grep("^A.+|^B.+|^C.+|^D.+",kegg_ko[,1]),1:2],stringsAsFactors=F)
ko$level1 <- NA
ko$level2 <- NA
ko$level3 <- NA

## level1
level1_list <- grep("^A",ko[,1])
n <- length(level1_list)
for(i in 1:(n-1)){{
ko$level1[level1_list[i]:(level1_list[i+1]-1)] <- ko[level1_list[i],1]
}}
ko$level1[level1_list[n]:nrow(ko)] <- ko[level1_list[n],1]

## level2
level2_list <- grep("^B",ko[,1])
n <- length(level2_list)
for(i in 1:(n-1)){{
ko$level2[level2_list[i]:(level2_list[i+1]-1)] <- ko[level2_list[i],1]
}}
ko$level2[level2_list[n]:nrow(ko)] <- ko[level2_list[n],1]

## level3
level3_list <- grep("^C",ko[,1])
n <- length(level3_list)
for(i in 1:(n-1)){{
ko$level3[level3_list[i]:(level3_list[i+1]-1)] <- ko[level3_list[i],1]
}}
ko$level3[level3_list[n]:nrow(ko)] <- ko[level3_list[n],1]

ko_path <- ko[grep("^D",ko[,1]),]
ko_path <- ko_path[grep(paste("PATH:",species_ko,sep=""),ko_path[,5]),]
ko_path[,1] <- sub(";.*","",ko_path[,1])
ko_path[,1] <- sub("^D +[0-9].*? ","",ko_path[,1])
ko_path[,2] <- sub(" +.*","",ko_path[,2])
ko_path[,3] <- sub("^A.*? ","",ko_path[,3])
ko_path[,4] <- sub("^B +[0-9].*? ","",ko_path[,4])
ko_path[,5] <- sub("^C +[0-9].*? ","",ko_path[,5])
ko_path[,5] <- sub("/"," ",ko_path[,5])
ko_path[,5] <- sub(" +"," ",ko_path[,5])
ko_path$pathway <- sub(".*:","",ko_path[,5])
ko_path$pathway <- sub("]","",ko_path$pathway)
ko_path[,5] <- sub(" \\\\[.*\\\\]","",ko_path[,5])
names(ko_path)[1:2] <- c("Gene","ko")
if(length(grep("^$",ko_path$Gene))>0){{
ko_path<-ko_path[-grep("^$",ko_path$Gene),]}}
## 是否去人类疾病及输出ko_path
if(del_human==T){{
ko_path <- ko_path[!(ko_path[,3]=="Human Diseases"),]
if(outkopath==T){{
write.table(ko_path,"ko_path_symbol_nohuman.xls",sep="\\t",quote=F,row.names=F)
}}
}} else {{
if(outkopath==T){{write.table(ko_path,"ko_path_symbol.xls",sep="\\t",quote=F,row.names=F)
}}
}}

### 4.hash构建
## 基因hash
pre_hash <- ko_path[,c(1,5)]
hash_gene <- unique(pre_hash[,1])

hash_gene_path <- hash()
for(i in 1:length(hash_gene)){{
.set(hash_gene_path,keys=hash_gene[i],values=list(pre_hash[pre_hash[,1]==hash_gene[i],2]))
}}

## 通路hash
hash_path <- unique(pre_hash[,2])
hash_path_gene <- hash()
for(i in 1:length(hash_path)){{
.set(hash_path_gene,keys=hash_path[i],values=list(unique(pre_hash[pre_hash[,2]==hash_path[i],1])))
}}

## 是否输出gmt
if(gmtout==T){{
gmt <- data.frame(keys(hash_path_gene),stringsAsFactors=F)
gmt$gene <- NA
for(i in 1:nrow(gmt)){{
gmt[i,]$gene <- paste(unlist(values(hash_path_gene,gmt[i,1])),collapse=";",sep="")
}}
gmt[,1] <- gsub("-","_",gmt[,1])
gmt$out <- paste(gmt[,1],gmt$gene,sep=";")
gmt_out <- data.frame(gmt$out,stringsAsFactors=F)
gmt_out <- gsub(";","\\t",gmt_out[,1])
write.table(gmt_out,paste(species_ko,"_",time,".gmt",sep=""),sep="\\t",quote=F,row.names=F,col.names=F)
}}

### 5.kegg富集
## 通路文件准备
gene_annotation <- unique(ko_path$Gene)

deg <- as.vector(degfile)
deg_annotation <- intersect(deg,gene_annotation)
if(length(deg_annotation)==0){{
  a="blank_file"
  write.table(a,paste(output,".kegg_enrichment.xls",sep=""),sep="\\t",col.names=F,row.names=F,quote=F)
 q()
}}
path_counts <- unique(unlist(values(hash_gene_path,deg_annotation)))
path_counts <- data.frame(path_counts,stringsAsFactors=F)
kegg_out <- ko_path[,c(3:6)]
kegg_out <- kegg_out[!duplicated(kegg_out[,3]),]
kegg_out <- merge(kegg_out,path_counts,by.x="level3",by.y="path_counts")
kegg_out <- kegg_out[,c(4,2,3,1)]
kegg_out$pathway_degs_num <- 0
degs_num_inpathway <- length(deg_annotation)
kegg_out$degs_num_inpathway <- degs_num_inpathway
kegg_out$pathway_gs_num_genome <- 0
gs_num_inpathway_genome <- length(gene_annotation)
kegg_out$gs_num_inpathway_genome <- gs_num_inpathway_genome
kegg_out$pvalue <- 1
kegg_out$FDR <- 1
kegg_out$zscore <- 0
kegg_out$pathway_degs <- NA
kegg_out$pathway_nondegs <- NA

##通路富集分析
for(i in 1:nrow(kegg_out)){{
path_gene <- unlist(values(hash_path_gene,kegg_out[i,]$level3))
deg_inpath <- intersect(deg_annotation,path_gene)
kegg_out[i,]$pathway_degs_num <- length(deg_inpath)
kegg_out[i,]$pathway_gs_num_genome <- length(path_gene)
kegg_out[i,]$pvalue <- phyper(kegg_out[i,]$pathway_degs_num-1,kegg_out[i,]$pathway_gs_num_genome,gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome,degs_num_inpathway,lower.tail=FALSE)
kegg_out[i,]$zscore <- (kegg_out[i,]$pathway_degs_num-degs_num_inpathway*kegg_out[i,]$pathway_gs_num_genome/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome)) / (sqrt(degs_num_inpathway*(kegg_out[i,]$pathway_gs_num_genome/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome))*(1-kegg_out[i,]$pathway_gs_num_genome/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome))*(1-(degs_num_inpathway-1)/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome-1))))
kegg_out[i,]$pathway_degs <- paste(deg_inpath,collapse=";")
kegg_out[i,]$pathway_nondegs <- paste(path_gene,collapse=";")
}}
kegg_out$FDR <- p.adjust(kegg_out$pvalue,method="fdr")
kegg_out <- kegg_out[order(kegg_out$pvalue),]
kegg_out_filter <- kegg_out[kegg_out$pvalue<p,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$FDR<fdr,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$zscore>zs,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$pathway_degs_num>=3,]
write.table(kegg_out_filter,paste(output,".kegg_enrichment.xls",sep=""),sep="\\t",row.names=F,quote=F)
write.table(kegg_out,paste(output,"unfiltered.kegg_enrichment.xls",sep=""),sep="\\t",row.names=F,quote=F)

### 6.绘图
## 数据准备
pathway <- kegg_out_filter
if(nrow(pathway) ==0){{
  cat(c("Warning ",output,".kegg_enrichment.xls is blank\\n"))
  q()
}}

if(nrow(pathway)>15) {{
data <- pathway[1:15,]
}} else {{
data <- pathway
}}



if(nrow(data)>10) {{
  h <- 8*nrow(data)/60
}} else if(nrow(data)<=5){{
  h <- 1.2
}}else{{
  h <- 1.6
}}

shorten_names <- function(x, n_word=4, n_char=40){{
  shortname = x["level3"]
  if (length(strsplit(shortname, " ")[[1]]) > n_word || (nchar(shortname) > n_char))
  {{
    if (nchar(shortname) > n_char) shortname <- substr(shortname, 1, n_char)
    y <- paste(paste(strsplit(shortname, " ")[[1]][1:min(length(strsplit(shortname," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(y)
  }} 
  else
  {{
    return (shortname)
  }}
}}

data$pretty_varname = apply(data,1,shorten_names)
max1<- max(nchar(c(as.character(data$pretty_varname))))
w<- max1*0.09


## 柱状图

bar <- ggplot(data,aes(y=pathway_degs_num,x=reorder(pretty_varname,-log10(pvalue)))) +
  geom_bar(aes(fill=-log10(pvalue)),stat="identity",width=0.8,position = position_dodge(0.9)) +
  scale_fill_gradientn(colours=c(rev(brewer.pal(11, "RdYlBu")[2:10]))) +
  #scale_fill_gradientn(colours=c(rev(brewer.pal(11, "Spectral")[1:5]))) +
  coord_flip() +labs(color=expression(-log10(pvalue)),y="Gene Number",x="Pathway") +
  theme_bw() +
  theme(
		title =element_text(size=6,face="plain",colour="black"),
		axis.title = element_text(size=6,face="plain",colour="black"),
		panel.grid=element_blank(),
      
		axis.line = element_line(size=0.25,colour = "black"),
		axis.text = element_text(color = "black",size = 6,face="plain"),
		
		legend.key.size=unit(0.3,'cm'),
		legend.title=element_text(size=6,colour = "black",face="plain"),
		legend.text = element_text(colour = 'black',face = 'plain',size=6)
		
)
ggsave(paste(output,".kegg_enrichment_bar.png",sep=""),bar,width=w,height=h)
ggsave(paste(output,".kegg_enrichment_bar.pdf",sep=""),bar,width=w,height=h)


## 分类图
colors <- length(data$level2[!duplicated(data$level2)])
category <-  ggplot(data,aes(x=-log10(pvalue),y=reorder(pretty_varname,-log10(pvalue)))) +
  geom_segment(aes(yend=reorder(pretty_varname,-log10(pvalue))),xend=0,colour="grey50") +
  geom_point(aes(size=pathway_degs_num,colour=level2)) +
  guides(colour=FALSE) +scale_size_area(max_size=5) +
  facet_grid(level2~.,scales="free_y",space="free_y") +
  labs(color=expression(-log10(pvalue)),size="Gene number",x="-log10(P-value)",y="Pathway") +
  theme_bw()+
  theme(strip.text.y=element_text(angle=0,size=6,colour="black"),
		strip.background = element_rect(fill=brewer.pal(11, "RdYlBu")[9]),
		legend.position="top",
        title =element_text(size=6,face="plain",colour="black"),
        axis.title = element_text(size=6,face="plain",colour="black"),
        panel.grid=element_blank(),
      
      axis.text = element_text(color = "black",size = 6,face="plain"),
      legend.title=element_text(size=6,colour = "black",face="plain"),
      legend.text = element_text(colour = 'black',face = 'plain',size=6)
  )
  
w=w*1.8
h=h*2
ggsave(paste(output,".kegg_enrichment_category.png",sep=""),category,width=w,height=h)
ggsave(paste(output,".kegg_enrichment_category.pdf",sep=""),category,width=w,height=h)

### 7.打印概况
cat("Genome annotation num of ",species,": ",gs_num_inpathway_genome,"\\n",sep="")
cat("DEG num: ",nrow(deg),"\\n",sep="")
cat("DEG annotation num of ",species,": ",degs_num_inpathway,"\\n",sep="")
cat(nrow(kegg_out)," pathways are enriched, and ",nrow(kegg_out_filter)," pathways are passed Fisher's exact test.","\\n",sep="")
###

    '''.format(
            species_type=dbs.get(species),
            species=species,
            output_prefix=output,
            input_str=genelist
        )
        with open("/home/user/Web/imav2/media/"+outname+".R",'w') as f:
            f.write(KEGG_enrichment_R)
        cmd = "cat /home/user/Web/imav2/media/"+outname+".R" + " |  R --vanilla"
        os.system(cmd)
        os.remove("/home/user/Web/imav2/media/"+outname+".R")
        result_dic = {'kegg_enrich': output + ".kegg_enrichment.xls",
                      "unfiltered_kegg_enrich": output + ".unfiltered..kegg_enrichment.xls",
                      "kegg_bar": outname + ".kegg_enrichment_bar.pdf",
                      "kegg_category": outname + ".kegg_enrichment_category.pdf",
                      }
        response = JsonResponse(result_dic)
        return response
    else:
        return HttpResponse("bad things happen")

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
            if os.path.exists('/home/user/Web/imav2/media/'+file_obj.name):
                os.remove('/home/user/Web/imav2/media/'+file_obj.name)
            f = open(os.path.join('/home/user/Web/imav2/media', file_obj.name), 'wb')
            for line in file_obj.chunks():
                f.write(line)
            f.close()
            genelistpath='/home/user/Web/imav2/media/'+file_obj.name
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
    #datapath="F:/GSE_Web-master"
    now = datetime.now()
    now = now.strftime('%Y%m%d%H')
    #rpath = "E:/IMA/web/imav2/result/" + now
    #rpath = "/mnt/e/IMA/web/imav2/result/" + now
    #rpath = "/media/sdc/yueyao/930/Web/test2/imav2/result/"+now
    #rpath ="F:/GSE_Web-master/result/"+now
    rpath = "/home/user/Web/imav2/result/"+now
    #rpath = "/media/sdc/yueyao/Web/imav2/result/"+now
    mediapath = os.path.dirname(rpath)+"/pic"
    #if os.path.exists(mediapath):
    #    shutil.rmtree(mediapath)
    #    os.makedirs(mediapath)
    #else:
    #    os.makedirs(mediapath)
    if not os.path.exists(rpath):
         os.makedirs(rpath)
    #    shutil.rmtree(rpath)
    #    os.makedirs(rpath)
    #else:
    #    os.makedirs(rpath)
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


def get_data(json_file):
    with open(json_file,'r') as f:
        dataset = json.load(f)
    return dataset

dataset = get_data("/home/user/Web/imav2/database.json")

dataset2={
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
                "GSE76",
                "GSE24242"
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


