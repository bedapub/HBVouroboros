import plotly as py
import pandas as pd
import numpy as np
import os
import sys
import glob
import argparse
from io import StringIO
import plotly.tools as plotly_tools
#import plotly.graph_objs as go
import plotly.graph_objects as go
import plotly.figure_factory as ff
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
import plotly.express as px 
from IPython.display import HTML





##################################################inferred genotype####################################################################
def infGen(infref_path, perSamp, inptRef):

    Inferred_genome=''' '''
    
    
    if perSamp == True:
        Inferred_genotype=[]
        matching = glob.glob(infref_path + "/**/infref_strain.fasta", recursive = True)
        matching.sort()
        matching=[i.split('/')[-2]+ "/" + i.split('/')[-1] for i in matching]
        for i in matching:
            infrefFasta=infref_path + "/" + i
            with open(infrefFasta) as f:
                first_line = f.readline()
            Inferred_genotype.append(first_line.split(" ")[0].replace(">",""))
        
        
        subfig = go.Figure(data=[go.Table(
                header=dict(values=['Sample', 'Genotype'],
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values=[[i.split('/')[-2] for i in matching], Inferred_genotype],
                            fill_color='lavender',
                            align='center', height=33))])       
        subfig.update_layout(title_text="Inferred genomes per sample", title_x=0.5)
        #subfig.layout.height=15*len(df)
        genotype_tab=py.offline.plot(subfig, include_plotlyjs=False, output_type='div') 
        Inferred_genome=Inferred_genome + "<h2>2 Inferred genomes</h3>\n"
        Inferred_genome=Inferred_genome + "<div class='scroll';align='center'>\n"
        Inferred_genome=Inferred_genome + "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n"
        Inferred_genome=Inferred_genome+genotype_tab+ "\n"
        Inferred_genome=Inferred_genome + "</div>\n"        
        
        
        
    else:
        if inptRef == True:
            infrefFasta=infref_path + "/inpt_strain.fasta"
            Inferred_genome=Inferred_genome + "<h2>2 Input genome</h2>\n"
        else:
            infrefFasta=infref_path + "/infref_strain.fasta"
            Inferred_genome=Inferred_genome + "<h2>2 Inferred genome</h2>\n"
        with open(infrefFasta) as f:
            first_line = f.readline()
        Inferred_genotype=first_line.split(" ")[0].replace(">","")
        Inferred_genome=Inferred_genome + "<ul>\n"    
        Inferred_genome=Inferred_genome+"<li> <p style='font-size: 15px; color: #0000FC'> <strong>" + Inferred_genotype +" </strong></p></li>\n"
        Inferred_genome=Inferred_genome + "</ul>\n"
    

   
    
    
    return(Inferred_genome)



##########################################################coverage plot of all samples################################################
def covplots(coverage_path, perSamp):
    
    
    
    layout = go.Layout(
    autosize=False,
    width=1200,
    height=600,

    xaxis= go.layout.XAxis(linecolor = 'black',
                          linewidth = 1,
                          mirror = True),

    yaxis= go.layout.YAxis(linecolor = 'black',
                          linewidth = 1,
                          mirror = True)
)

    
    if perSamp == True:
        data=[]
        matching = glob.glob(coverage_path + "/perSamp/*genome_depth.tsv", recursive = True)
        matching.sort()
        matching=[i.split('/')[-2]+ "/" + i.split('/')[-1] for i in matching]
        for i in matching:
            covFile=coverage_path + "/" + i
            with open(covFile) as f:
                df = pd.read_csv(f, sep="\t", header=0, index_col=0)
                covLine = go.Scatter( x=df.iloc[:,0], y=df.iloc[:,1], name=df.columns[1]) 

                data.append(covLine)

        cov_plot = go.Figure(data = data,layout= layout)
        cov_plot.update_layout(legend_title_text='Sample')
        cov_plot.update_layout( width=1000, height=600)
        cov_plot.update_layout(xaxis_title="Position", yaxis_title="Read count")
        cov_plot=py.offline.plot(cov_plot, include_plotlyjs=False, output_type='div')  

        coverage_str=''' '''
        coverage_str=coverage_str + "<h2>4 Per-nucleotide count</h2>\n"
        coverage_str=coverage_str + "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n"
        coverage_str=coverage_str + cov_plot + "\n"    
    
    else:
        cov=pd.read_csv(coverage_path, delimiter='\t')
        cov.columns=cov.columns.str.replace('infref_','')
        fig = px.line(cov, x=cov['POS'], y=cov.columns[2:])

        fig.update_layout(legend_title_text='Sample')
        fig.update_layout( width=1000, height=600, xaxis= go.layout.XAxis(linecolor = 'black', linewidth = 1, mirror = True), yaxis= go.layout.YAxis(linecolor = 'black', linewidth = 1, mirror = True), autosize=False)
        fig.update_layout(xaxis_title="Position", yaxis_title="Read count")
    
        cov_plot=py.offline.plot(fig, include_plotlyjs=False, output_type='div')    
        coverage_str=''' '''
        coverage_str=coverage_str + "<h2>4 Per-nucleotide count</h2>\n"
        coverage_str=coverage_str + "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n"
        coverage_str=coverage_str + cov_plot + "\n"
    return(coverage_str)





#####################################################count table#######################################################################
def covTable(totalCount_path, perSamp):

         
    if perSamp == True:
        counts=[]
        samples=[]
        matching = glob.glob(totalCount_path + "/perSamp/*genome_count.tsv", recursive = True)
        matching.sort()
        matching=[i.split('/')[-2]+ "/" + i.split('/')[-1] for i in matching]
        for i in matching:
            infrefFasta=totalCount_path + "/" + i
            with open(infrefFasta) as f:
                df = pd.read_csv(f, sep="\t", header=0, index_col=None)
                counts.append(df.iloc[0,1])
                samples.append(i.split('/')[-2])
            #Inferred_genotype.append(first_line.split(" ")[0].replace(">",""))
        
        subfig = go.Figure(data=[go.Table(
                    header=dict(values=['Sample', 'Genotype'],
                                fill_color='paleturquoise',
                                align='left'),
                    cells=dict(values=[samples, counts],
                                fill_color='lavender',
                                align='center', height=33))])       
        subfig.update_layout(title_text="Count per sample", title_x=0.5)
        #subfig.layout.height=15*len(df)
        count_table=py.offline.plot(subfig, include_plotlyjs=False, output_type='div')
            
            
    
    
    else:
        df = pd.read_csv(totalCount_path, sep="\t", names=['ID', 'Coverage'], index_col=None)

        subfig = go.Figure(data=[go.Table(
            header=dict(values=['ID', 'Coverage'],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[df.ID[1:].str.split('infref_').str[-1], df.Coverage[1:]],
                        fill_color='lavender',
                        align='center', height=33))])
        subfig.update_layout(title_text="Count per sample", title_x=0.5)
        #subfig.layout.height=15*len(df)
        count_table=py.offline.plot(subfig, include_plotlyjs=False, output_type='div')
    Count_table=''' '''
    Count_table=Count_table + "<h2>3 Count tables</h2>\n"

    Count_table=Count_table + "<div class='scroll';align='center'>\n"
    Count_table=Count_table + "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n"
    Count_table=Count_table+count_table+ "\n"
    Count_table=Count_table + "</div>\n"
    return(Count_table)


##########################################################SNP plots per sample#########################################################

def SNP_Plots(vc_path, perSamp):


    data=[] 
    SNP_plots=[]
    SNP_plot=''''''
   
    if perSamp == True:
        matching = glob.glob(vc_path + "/**/*primitive*", recursive = True)
        matching.sort()
        matching=[i.split('/')[-1].replace('/','') for i in matching]
        matching=[i.split('_cleaned')[0]+'/'+i for i in matching]
    else:        
        mm=os.listdir(vc_path)
        matching = [s for s in mm if "primitive" in s]
        matching.sort()
        
    for j in range(0,len(matching)):
        vcFile=open(vc_path+'/'+matching[j],'r')
        vcfLines=vcFile.readlines()
        aline=[]
        for i in vcfLines:
            #possibly fragile
            if (i[0] != '#'):
                aline.append(i)
        
        lines="".join(aline)
        df = pd.read_csv(StringIO(lines), sep="\t", names=['Genome', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','unknown'], index_col=None)
    
    
        subfig = go.Figure(data=[go.Table(
            header=dict(values=['Sample','POS','REF','ALT', 'QUAL'],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[df.Genome, df.POS, df.REF, df.ALT, df.QUAL],
                       fill_color='lavender',
                       align='center', height=33))])
        #df=df.iloc[:,0:6]            
        #subfig =  ff.create_table(df,height_constant=20)
        #thrFig=py.offline.plot(fig,  output_type='div')
        SNP_plots.append(subfig)
                
                
                
              
        subfig.data[0]['columnwidth'] = [10,2,2,2,2];
        subfig.update_layout(width=600)
        subfig.update_layout(title_text=matching[j].split('_cleaned')[0].replace('infref_',''), title_x=0.5)
        subfig.layout.width=800
        #subfig.show()
        if len(df) > 0:
            SNP_figs=py.offline.plot(subfig, include_plotlyjs=False, output_type='div')
            data.append(SNP_figs)
   

    SNP_plot="<h2>5 SNP tables</h2>\n" + SNP_plot
    SNP_plot=SNP_plot + "<div class='scroll';>\n" 
    for i in range(0,len(data)):
        SNP_plot=SNP_plot + "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n"
        SNP_plot=SNP_plot+data[i]+ "\n"
    SNP_plot=SNP_plot + "</div>\n" 
    return(SNP_plot)



######################################################multiQC path#####################################################################

def multiQC_link(multiqc_path) :  

    multiqc_str=''' '''
    multiqc_str=multiqc_str + "<ul>\n"
    multiqc_str =  multiqc_str +"<li class='item'> <a href=file://" + str(multiqc_path) + " style='font-size: 15px; color: color: #0000FC'><strong>Multi_QC report<strong></a> </li>\n"
    multiqc_str=multiqc_str+"</ul>\n"
    return(multiqc_str)
       
    
    
    
    
    
    

    

def main(args):

    parser = argparse.ArgumentParser(
        description='HBVouroboros html summary report(s) generator')
    parser.add_argument('mode',
        help = 'persamp, inpt or infref')
    parser.add_argument('outputFile',
        help = 'output file to be saved')
    args = parser.parse_args(args)


    run_mode = args.mode
    outputFile = str(args.outputFile)


    #print(run_mode)
    if (run_mode != 'infref' and run_mode != 'inpt' and run_mode != 'persamp'):
        print("Bad input: choose one of the valid modes: infref, inpt, persamp")
        return(False)
        sys.exit()

    #multiqc_path='results/multiqc/infref/infref_multiqc_report.html'
    #multiqc_input_path='results/multiqc/inpt/inpt_multiqc_report.html'
    #multiqc_persamp_path='results/multiqc/perSamp/perSamp_multiqc_report.html'

    multiqc_path=os.path.realpath(os.path.join(os.path.dirname("__file__"), 'results', 'multiqc', 'infref','infref_multiqc_report.html'))
    multiqc_input_path=os.path.realpath(os.path.join(os.path.dirname("__file__"), 'results', 'multiqc', 'inpt','inpt_multiqc_report.html'))
    multiqc_persamp_path=os.path.realpath(os.path.join(os.path.dirname("__file__"), 'results', 'multiqc', 'perSamp','perSamp_multiqc_report.html'))

    coverage_path='results/coverage/infref/infref_genome_depth.tsv'
    coverage_persamp_path='results/coverage'
    coverage_inpt_path='results/coverage/inpt/inpt_genome_depth.tsv'

    totalCount_path='results/coverage/infref/infref_genome_count.tsv'
    totalCount_persamp_path='results/coverage'
    totalCount_inpt_path='results/coverage/inpt/inpt_genome_count.tsv'

    vc_path='results/variant-calling/infref'
    vc_persamp_path='results/variant-calling/perSamp'
    vc_inpt_path='results/variant-calling/inpt'

    infref_path='results/infref'
    infref_persamp_path='results/perSamp'
    infref_inpt_path='results/inpt'

    
    
    
    
    theReport = '''
    <!DOCTYPE html>
    <html>
    <head>
    <style>

    .flex-container {
        overflow-x: auto;
        display: flex;
        flex-wrap: nowrap;
    }
    
    .flex-container > div {

        background-color: #f1f1f1;
        margin: 10px;
        text-align: center;
        line-height: 75px;
        font-size: 30px;
    }



    div.scroll {
        padding:4px;
        width: 900px;
         height: 500px;
         overflow-x: hidden;
         overflow-y: auto;
        text-align:justify;
        border: #000000 1px solid;
        align: center;
    }

    a {
        text-decoration: none;
    }
    a:link, a:visited {
        color: #0000FC;
    }
    a:hover {
    color: #0000FC;
    }

    </style>
    </head>
    <body>
        <head>
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
'''

    if (run_mode == "infref"):
    
        theReport=theReport+"<h1 style='font-size:35px;'>HBVouroboros pipeline summary report using dataset-level-inferred reference genotype</h1>\n"
    if (run_mode == "inpt"):   
        theReport=theReport+"<h1 style='font-size:35px;'>HBVouroboros pipeline summary report using user-input reference genotype</h1>\n" 
    if (run_mode == "persamp"):
        theReport=theReport+"<h1 style='font-size:35px;'>HBVouroboros pipeline summary report using per-sample-inferred reference genotype</h1>\n"
    theReport=theReport+'''

                <hr>
                <br> <br>
                <h2>1 QC reports</h2>

'''
 
       
    if (run_mode == "infref"):
        theReport = theReport + multiQC_link(multiqc_path) + infGen(infref_path, False, False) + covTable(totalCount_path, False) + covplots(coverage_path, False) + SNP_Plots(vc_path, False)

    if (run_mode == "persamp"):
        theReport = theReport +multiQC_link(multiqc_persamp_path) + infGen(infref_persamp_path, True, False) + covTable(totalCount_persamp_path, True) + covplots(coverage_persamp_path, True) + SNP_Plots(vc_persamp_path, True)

    if (run_mode == "inpt"):
        theReport = theReport + multiQC_link(multiqc_input_path) + infGen(infref_inpt_path, False, True) + covTable(totalCount_inpt_path, False) + covplots(coverage_inpt_path, False) + SNP_Plots(vc_inpt_path, False)




    f = open(outputFile,'w')
    f.write(theReport)
    f.close()
    return True
if __name__ == '__main__':
    main(sys.argv[1:])
