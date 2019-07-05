from flask import Flask, render_template, flash, request, redirect, url_for, send_from_directory, session
from app import app
from app.forms import NumberOfSamplesForm
import pandas as pd
import numpy as np
import json
import plotly
import plotly.graph_objs as go
import plotly.io as pio
import cufflinks as cf
import io
import os
import seaborn as sns
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import base64
from werkzeug.utils import secure_filename
from flask_wtf import FlaskForm 
from flask_wtf.file import FileField
from wtforms import StringField, PasswordField, BooleanField, SubmitField, IntegerField, TextField, HiddenField
from wtforms.validators import DataRequired
import urllib.parse

#Function to read and process deleterious variants file
def read_var_file(path):
    deleterious_variants = pd.read_csv(path, sep='\t')
    del_var = deleterious_variants.iloc[: ,[0,1,3,4,5,6]]
    
    proper_transcript = []
    variant_type = []
    reference_variant = []
    alternate_variant = []

    for i in del_var['Transcript Variant']:
        for j in i.replace(" ", "").split(";"):        
            #Pick first entry with c and exit
            if j[0] == "c":
                proper_transcript.append(j)
                if len(j.split(">")) > 1:
                    variant_type.append("SNV")
                    reference_variant.append(j.split(">")[0][-1])
                    alternate_variant.append(j.split(">")[1])
                elif len(j.split("del")) > 1:
                    variant_type.append("del")
                    reference_variant.append(j.split("del")[1])
                    alternate_variant.append("-")
                elif len(j.split("dup")) > 1:  
                    variant_type.append("ins")
                    reference_variant.append("-")
                    alternate_variant.append(j.split("dup")[1]) 
                elif len(j.split("ins")) > 1:
                    variant_type.append("ins")
                    reference_variant.append("-")
                    alternate_variant.append(j.split("ins")[1])
                else:
                    variant_type.append("-")
                    reference_variant.append("-")
                    alternate_variant.append("-")
                break
    
    #Add lists as new columns
    del_var = del_var.assign(Proper_Transcript=proper_transcript)
    del_var = del_var.assign(Variant_Type=variant_type)
    del_var = del_var.assign(Reference_Variant=reference_variant)
    del_var = del_var.assign(Alternate_Variant=alternate_variant)

    #Deletions need to have position -1
    del_var.loc[del_var['Variant_Type'] == 'del', 'Position'] -= 1
    return del_var
    
#Function to read vcf file    
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

#Function to read gVCDF and select variants
def select_genomic_variants(path, concat_del_vars, date, sample_name):
    genomic_vcf = read_vcf(path)
    genomic_vcf['CHROM'] = genomic_vcf['CHROM'].str.replace('chr', '')

    #Select variants from genomic VCF file
    selected_variants = pd.merge(concat_del_vars, genomic_vcf,  how='left', left_on=['Chromosome','Position'], right_on = ['CHROM','POS'])

    #Select wanted columns
    selected_variants = selected_variants.iloc[: ,[0,1,2,3,4,5,6,7,8,-3,-1]]

    #Rename sample unique column to generic name, add Unique_Symbol column, split info field
    selected_variants = selected_variants.rename(columns={ selected_variants.columns[-1]: "Sample_Info" })
    selected_variants['Unique_Symbol'] = selected_variants['Gene Symbol']+"_"+selected_variants['Proper_Transcript']
    selected_variants['VAF'] = selected_variants['Sample_Info'].str.split(":", expand = True)[3].astype('float64').round(4)
    selected_variants['INFO'] = selected_variants['INFO'].str.replace('DP=', '')
    selected_variants['Depth'] = selected_variants['INFO'].str.split(";", expand = True)[0]
    selected_variants['Date'] = date
    selected_variants['SampleName'] = sample_name
    
    #Sort by VAF in order to remove duplicates and only keep the variant with largets VAF
    selected_variants = selected_variants.sort_values(by='VAF').drop_duplicates(subset=['Chromosome', 'Position'], keep = 'last')

    return selected_variants

@app.route('/')
@app.route('/index', methods=['GET', 'POST'])
def index():

    form = NumberOfSamplesForm()

    if form.validate_on_submit():
        session['numberofsamples']=form.numberofsamples.data
        session['samplename']=form.samplename.data
        return redirect(url_for('submitinfo'))

    return render_template('index.html', form=form)


@app.route('/submitinfo', methods=['GET', 'POST'])
def submitinfo():
    
    numberofsamples = session['numberofsamples']
    tidpunkter = range(1,int(numberofsamples)+1)

    class MyStandardForm(FlaskForm):
        pass
    
    class F(MyStandardForm):
        pass

    for name in tidpunkter:
        date = 'Date_'+str(name)
        sample_name = 'SampleName_'+str(name)
        del_file = 'DelFile_'+str(name)
        gvcf_file = 'gVCFFile_'+str(name)
        setattr(F, date, TextField(date.title(), validators=[DataRequired()]))
        setattr(F, sample_name, TextField(sample_name.title(), validators=[DataRequired()]))
        setattr(F, del_file, FileField(del_file.title(), validators=[DataRequired()]))
        setattr(F, gvcf_file, FileField(gvcf_file.title(), validators=[DataRequired()]))
    setattr(F, 'submit', SubmitField('Fortsätt'))
    form = F()
    
    if request.method == 'POST':

        dates=[]
        sample_names=[]
        del_files=[]
        gvcf_files=[]

        for i in tidpunkter:
            date = 'Date_'+str(i)
            sample_name = 'SampleName_'+str(i)
            del_file = 'DelFile_'+str(i)
            gvcf_file = 'gVCFFile_'+str(i)
            
            dates.append(request.form[date])
            sample_names.append(request.form[sample_name])

            f = request.files[del_file]
            f.save('uploads/' + secure_filename(f.filename))
            del_files.append(secure_filename(f.filename))

            g = request.files[gvcf_file]
            g.save('uploads/' + secure_filename(g.filename))
            gvcf_files.append(secure_filename(g.filename))

        return redirect(url_for('results', dates=dates, sample_names=sample_names, del_files=del_files, gvcf_files=gvcf_files)) 
        
    return render_template('submitinfo.html', form=form)

@app.route('/results', methods=("POST", "GET"))
def results():	

    if request.method == 'POST':
        plot_info=request.form.getlist('use_in_plot')
        return redirect(url_for('final_results', plot_info=plot_info))
    
    numberofsamples = session['numberofsamples']
    
    del_files=request.args.getlist('del_files')
    gvcf_files=request.args.getlist('gvcf_files')
    sample_names=request.args.getlist('sample_names')
    dates=request.args.getlist('dates')

    # Set session variables
    session['dates']=request.args.getlist('dates')
    session['sample_names']=request.args.getlist('sample_names')
    session['del_files']=request.args.getlist('del_files')
    session['gvcf_files']=request.args.getlist('gvcf_files')

    dates_list = request.args.getlist('dates')
    sample_name_list = request.args.getlist('sample_names')
    
    del_files_list=['uploads/' + x for x in request.args.getlist('del_files')]
    gvcf_files_list=['uploads/' + x for x in request.args.getlist('gvcf_files')]

    setup_file = pd.DataFrame(
    {'Date': dates_list,
     'SampleName': sample_name_list,
     'Del_File': del_files_list,
     'gVCF_File': gvcf_files_list
    })
    #Put all del_var panda dfs into dict and concatenate
    del_files2={}
    for x in range(0,len(setup_file['Del_File'])):
        del_files2["del_var_{0}".format(x+1)]=read_var_file(setup_file.loc[x,'Del_File'])
    
    concat_del_var = pd.concat(del_files2, ignore_index=True)

    #Convert Chromosome column to str to enable duplication detection
    concat_del_var['Chromosome'] = concat_del_var['Chromosome'].astype(str)
    concat_del_var = concat_del_var.drop_duplicates(subset=['Chromosome', 'Position'], keep = 'first')

    gvcf_variants={}
    for x in range(0,len(setup_file['gVCF_File'])):
        gvcf_variants["gvcf_var_{0}".format(x+1)]=select_genomic_variants(setup_file.loc[x,'gVCF_File'],concat_del_var, setup_file.loc[x,'Date'], setup_file.loc[x,'SampleName'])
    concat_plot_frame = pd.concat(gvcf_variants, ignore_index=True)

    concat_plot_frame['VAF'] = pd.to_numeric(concat_plot_frame['VAF'])*100
    concat_plot_frame['Date'] = pd.to_datetime(concat_plot_frame['Date'])
    concat_plot_frame = concat_plot_frame.sort_values(by=['Date'])

    #Group by date to get separate dfs for each occasion
    df_rn = concat_plot_frame.groupby(['Date'])

    #Extract VAF and Depth for each timepoint
    final_dict={}
    for x in range(0,len([df_rn.get_group(x) for x in df_rn.groups])):
        final_dict["VAF_{0}".format(x+1)]=([df_rn.get_group(x) for x in df_rn.groups][x].sort_values(by=['Chromosome', 'Position']).reset_index()['VAF'])
        final_dict["Depth_{0}".format(x+1)]=([df_rn.get_group(x) for x in df_rn.groups][x].sort_values(by=['Chromosome', 'Position']).reset_index()['Depth'])

    #Use the concated deleterious variants as report base, sort by chrom and pos to get same order
    concat_del_var_sorted = concat_del_var.sort_values(by=['Chromosome', 'Position']).reset_index()

    #Add VAF and Depth columns
    for i in final_dict:
        concat_del_var_sorted[i] = pd.Series(final_dict[i])

    #Long bit for renaming the headers. Find a better solution for this junk.
    corrected_headers = {}
    date_list = concat_plot_frame.Date.unique()
    sample_name_list = concat_plot_frame.SampleName.unique()
    n=0
    
    for i in range(numberofsamples*2*-1,-1,2):
        new_col = []

        VAFS = concat_del_var_sorted.iloc[:, i].values.round(2)
        depths = concat_del_var_sorted.iloc[:, i+1].values
        
        for x in range(0,len(VAFS)):
            new_col.append(str(VAFS[x])+' ('+str(depths[x])+')')

        corrected_headers[np.datetime_as_string(date_list[n], unit='D')+' - '+sample_name_list[n]]=new_col
        n+=1

    for i in corrected_headers:
        concat_del_var_sorted[i] = pd.Series(corrected_headers[i])

    # Print the variants to file 
    output_variants = concat_del_var_sorted.iloc[:, 1:10].copy()
    output_variants= pd.concat([output_variants, concat_del_var_sorted.iloc[:, numberofsamples*-1:]], axis=1, join='inner')


    output_variants.to_csv('app/static/Output/'+session['samplename']+'_temp.txt',index=False, sep="\t", decimal=',')
    csv_url=session['samplename']+'_temp.txt'

    #Select a subset of the fields to display on results page
    print_variants = concat_del_var_sorted.iloc[:, [3,5,6]].copy()
    print_variants= pd.concat([print_variants, concat_del_var_sorted.iloc[:, numberofsamples*-1:]], axis=1, join='inner')
    print_variants['Plot_Variant'] = ""

    for i in range(0, print_variants.shape[0]):
        print_variants['Plot_Variant'][i] = '<input type="checkbox" name="use_in_plot" value="'+str(i)+'" checked>'

    pd.set_option('display.max_colwidth', -1)


    ### PLOTLY ###
    
    unique_variants = concat_plot_frame.Unique_Symbol.unique()
    unique_dates = concat_plot_frame.Date.unique()
        
    fig = {
    'data': [
        {
            'x': concat_plot_frame[concat_plot_frame['Unique_Symbol']==variant]['Date'],
            'y': concat_plot_frame[concat_plot_frame['Unique_Symbol']==variant]['VAF'],
            'name': variant, 'mode': 'lines+markers',
        } for variant in unique_variants
    ],
    'layout': {
        'title': session['samplename'],
        'xaxis': {'type': 'date'},
        'yaxis': {'title': "Variantallelfrekvens (VAF)"},
        'legend': {'orientation': 'h'}
    }
    }
    cf.set_config_file(offline=True, world_readable=False)

    div = plotly.offline.plot(fig, show_link=False, output_type="div", include_plotlyjs=False)

    return render_template('results.html', tables=[print_variants.to_html(classes='u-max-full-width', index=False, escape=False)], csv_url=csv_url, plotly_graph=div, dates=dates, del_files=del_files, gvcf_files=gvcf_files, title=session['samplename'])


@app.route('/final_results', methods=("POST", "GET"))
def final_results():	
    
    numberofsamples = session['numberofsamples']
    
    dates=session['dates']
    sample_names=session['sample_names']
    del_files=session['del_files']
    gvcf_files=session['gvcf_files']
    plot_info=request.args.getlist('plot_info')
    #convert str to int
    plot_info = [int(i) for i in plot_info]
    
    dates_list = dates
    sample_name_list = sample_names
    del_files_list=['uploads/' + x for x in del_files]
    gvcf_files_list=['uploads/' + x for x in gvcf_files]

    setup_file = pd.DataFrame(
    {'Date': dates_list,
     'SampleName': sample_name_list,
     'Del_File': del_files_list,
     'gVCF_File': gvcf_files_list
    })
    #Put all del_var panda dfs into dict and concatenate
    del_files2={}
    for x in range(0,len(setup_file['Del_File'])):
        del_files2["del_var_{0}".format(x+1)]=read_var_file(setup_file.loc[x,'Del_File'])
    
    concat_del_var = pd.concat(del_files2, ignore_index=True)

    #Convert Chromosome column to str to enable duplication detection
    concat_del_var['Chromosome'] = concat_del_var['Chromosome'].astype(str)
    concat_del_var = concat_del_var.drop_duplicates(subset=['Chromosome', 'Position'], keep = 'first')

    gvcf_variants={}
    for x in range(0,len(setup_file['gVCF_File'])):
        gvcf_variants["gvcf_var_{0}".format(x+1)]=select_genomic_variants(setup_file.loc[x,'gVCF_File'],concat_del_var, setup_file.loc[x,'Date'], setup_file.loc[x,'SampleName'])
    concat_plot_frame = pd.concat(gvcf_variants, ignore_index=True)

    concat_plot_frame['VAF'] = pd.to_numeric(concat_plot_frame['VAF'])*100
    concat_plot_frame['Date'] = pd.to_datetime(concat_plot_frame['Date'])
    concat_plot_frame = concat_plot_frame.sort_values(by=['Date'])

    #Group by date to get separate dfs for each occasion
    df_rn = concat_plot_frame.groupby(['Date'])

    #Extract VAF and Depth for each timepoint
    final_dict={}
    for x in range(0,len([df_rn.get_group(x) for x in df_rn.groups])):
        final_dict["VAF_{0}".format(x+1)]=([df_rn.get_group(x) for x in df_rn.groups][x].sort_values(by=['Chromosome', 'Position']).reset_index()['VAF'])
        final_dict["Depth_{0}".format(x+1)]=([df_rn.get_group(x) for x in df_rn.groups][x].sort_values(by=['Chromosome', 'Position']).reset_index()['Depth'])

    #Use the concated deleterious variants as report base, sort by chrom and pos to get same order
    concat_del_var_sorted = concat_del_var.sort_values(by=['Chromosome', 'Position']).reset_index()

    #Add VAF and Depth columns
    for i in final_dict:
        concat_del_var_sorted[i] = pd.Series(final_dict[i])

    #Long bit for renaming the headers. Find a better solution for this junk.
    corrected_headers = {}
    date_list = concat_plot_frame.Date.unique()
    sample_name_list = concat_plot_frame.SampleName.unique()
    n=0
    
    for i in range(numberofsamples*2*-1,-1,2):
        new_col = []

        VAFS = concat_del_var_sorted.iloc[:, i].values.round(2)
        depths = concat_del_var_sorted.iloc[:, i+1].values
        
        for x in range(0,len(VAFS)):
            new_col.append(str(VAFS[x])+' ('+str(depths[x])+')')

        corrected_headers[np.datetime_as_string(date_list[n], unit='D')+' - '+sample_name_list[n]]=new_col
        n+=1

    for i in corrected_headers:
        concat_del_var_sorted[i] = pd.Series(corrected_headers[i])


    ### CORRECT FOR SELECTED VARIANTS ###
    concat_del_var_sorted = concat_del_var_sorted.iloc[plot_info]

    # Print the variants to file 
    output_variants = concat_del_var_sorted.iloc[:, 1:10].copy()
    output_variants= pd.concat([output_variants, concat_del_var_sorted.iloc[:, numberofsamples*-1:]], axis=1, join='inner')


    output_variants.to_excel('app/static/Output/'+session['samplename']+'.xlsx',index=False)
    csv_url=session['samplename']+'.txt'

    #Select a subset of the fields to display on results page
    print_variants = concat_del_var_sorted.iloc[:, [3,5,6]].copy()
    print_variants= pd.concat([print_variants, concat_del_var_sorted.iloc[:, numberofsamples*-1:]], axis=1, join='inner')

    ### PLOTLY ###
    
    concat_plot_frame = concat_plot_frame.loc[concat_plot_frame['Transcript Variant'].isin(concat_del_var_sorted['Transcript Variant'])]

    unique_variants = concat_plot_frame.Unique_Symbol.unique()
    unique_dates = concat_plot_frame.Date.unique()
        
    fig = {
    'data': [
        {
            'x': concat_plot_frame[concat_plot_frame['Unique_Symbol']==variant]['Date'],
            'y': concat_plot_frame[concat_plot_frame['Unique_Symbol']==variant]['VAF'],
            'name': variant, 'mode': 'lines+markers',
        } for variant in unique_variants
    ],
    'layout': {
        'title': session['samplename'],
        'xaxis': {'type': 'date'},
        'yaxis': {'title': "Variantallelfrekvens (VAF)", 'rangemode' : 'tozero'},
        'legend': {'orientation': 'h'}
    }
    }
    cf.set_config_file(offline=True, world_readable=False)

    div = plotly.offline.plot(fig, show_link=False, output_type="div", include_plotlyjs=False)
    pio.write_image(fig, 'app/static/Output/'+session['samplename']+'.pdf')
    #########

    return render_template('final_results.html', tables=[print_variants.to_html(classes='u-max-full-width', index=False, escape=False)], csv_url=csv_url, plotly_graph=div, dates=dates, del_files=del_files, gvcf_files=gvcf_files, title=session['samplename'])


if __name__ == '__main__':
    app.run(debug=True)