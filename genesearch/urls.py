#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     urls.py
   Description :
   Author :       yueyao
   date：          2019/9/27
-------------------------------------------------
   Change Activity:
                   2019/9/27:
-------------------------------------------------
"""
__author__ = 'yueyao'

from django.urls import path,re_path
from . import views
#from django.conf import settings
#from django.conf.urls.static import static

urlpatterns=[
    path('',views.index),
    path('query',views.detail2),
    path('gene_general',views.general),
    path('gse_query',views.gsequery),
    path('general_gene_query',views.general_genes),
    path('general_genes_table',views.general_genes_table),
    path('single_gene',views.single_gene),
    path('single_gene_table',views.single_gene_table),
    path('upload_gene',views.upload_genes),
    path('showtable',views.showtable),
    re_path(r'download',views.download_file,name='genesearch'),
    re_path(r'^GetSpecies/$', views.Return_Species_Data),
    re_path(r'^GetTissue/$', views.Return_Tissue_Data),
    path('gene_general2',views.general2),
    path('barplot',views.barplot),
    path('heatmap',views.heatmap),
    path('volmap',views.volmap),
    path('upload',views.upload),
    path('single',views.single)
]