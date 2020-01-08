from django.db import models

# Create your models here.

from django.db import models

# Create your models here.
class Disease(models.Model):
    type=models.CharField(max_length=30)
    def __str__(self):
        return self.type

class Species(models.Model):
    disease=models.ForeignKey(Disease,on_delete=models.CASCADE)
    type=models.CharField(max_length=30)
    def __str__(self):
        return self.type

class Tissue(models.Model):
    type=models.CharField(max_length=30)
    disease=models.ForeignKey(Disease,on_delete=models.CASCADE)
    species=models.ForeignKey(Species,on_delete=models.CASCADE)
    def __str__(self):
        return self.type

class Articles(models.Model):
    title=models.CharField(max_length=32)
    content=models.TextField()
    def __str__(self):
        return self.title

class DEGs(models.Model):
    databasetype=models.CharField(max_length=32)
    gene_id=models.CharField(max_length=32)
    fold=models.CharField(max_length=32)
    pvalue=models.CharField(max_length=32)
    padjust=models.CharField(max_length=32)
    def __str__(self):
        return self.databasetype

class Resultpath(models.Model):
    resultpath=models.CharField(max_length=32)
    expfilepath=models.CharField(max_length=32)
    difffilepath=models.CharField(max_length=32)
    inputgenenum=models.CharField(max_length=32)
    picpath=models.CharField(max_length=32)
    def __str__(self):
        return self.resultpath

# from django.contrib.auth.models import AbstractUser
# class MyUser(AbstractUser):
#     username=models.CharField('username',max_length=20)
#     userpwd=models.CharField('userpwd',max_length=20)
#     def __str__(self):
#         return self.username,self.userpwd
