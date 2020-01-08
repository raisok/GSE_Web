# Generated by Django 2.2.1 on 2019-09-29 02:31

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Articles',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=32)),
                ('content', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='DEGs',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('databasetype', models.CharField(max_length=32)),
                ('gene_id', models.CharField(max_length=32)),
                ('fold', models.CharField(max_length=32)),
                ('pvalue', models.CharField(max_length=32)),
                ('padjust', models.CharField(max_length=32)),
            ],
        ),
        migrations.CreateModel(
            name='Disease',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type', models.CharField(max_length=30)),
            ],
        ),
        migrations.CreateModel(
            name='Resultpath',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('resultpath', models.CharField(max_length=32)),
                ('expfilepath', models.CharField(max_length=32)),
                ('difffilepath', models.CharField(max_length=32)),
                ('inputgenenum', models.CharField(max_length=32)),
                ('picpath', models.CharField(max_length=32)),
            ],
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type', models.CharField(max_length=30)),
                ('disease', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genesearch.Disease')),
            ],
        ),
        migrations.CreateModel(
            name='Tissue',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('type', models.CharField(max_length=30)),
                ('disease', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genesearch.Disease')),
                ('species', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genesearch.Species')),
            ],
        ),
    ]
