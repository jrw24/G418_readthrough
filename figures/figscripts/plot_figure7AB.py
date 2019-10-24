
### plotting TE's

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import math
from scipy import stats
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
pd.set_option('display.max_columns', 40)
import argparse
import importlib

parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)
import rphelper as rph

rootDir = args.rootDir
threadNumb = str(args.threadNumb)

### set colors 

# outfileOrder = ['7A_right', '7A_left', '7S1A']

inputDir = "%s/Data/RNA/FPassignment/hg38_protCode/HEK_G418/analysis/Xtail" %  rootDir

colorList = ['#000000', '#ffb000', '#63cfff', '#eb4300', '#00c48f', '#eb68c0', '#fff71c', '#006eb9']
colorList= colorList
colorDict = {
   'black':colorList[0],
    'orange':colorList[1],
    'cyan':colorList[2],
    'red':colorList[3],
    'green':colorList[4],
    'pink':colorList[5],
    'yellow':colorList[6],
    'blue':colorList[7] 
}

def log_trans_b10(x):
    try:
        return math.log(x, 10)
    except:
#         return float(-6.00)
        return float("NaN")
def log_trans_b2(x):
    try:
        return math.log(x, 2)
    except:
#         return float("NaN")
        return float(-15.00) # set arbitrarily low value



def load_Xtail_results():
	d1 = pd.read_csv('%s/Untr_vs_g4182K_XtailResTab.csv' % inputDir)
	d1 = d1.rename(columns = {'Unnamed: 0':'index'})
	d1.set_index('index', inplace=True)
	d1.head()

	d2 = pd.read_csv('%s/Untr_vs_g41810min_XtailResTab.csv' % inputDir)
	d2 = d2.rename(columns = {'Unnamed: 0':'index'})
	d2.set_index('index', inplace=True)
	d2.head()

	d3 = pd.read_csv('%s/Untr_vs_g418_500_XtailResTab.csv' % inputDir)
	d3 = d3.rename(columns = {'Unnamed: 0':'index'})
	d3.set_index('index', inplace=True)
	d3.head()

	return d1, d2, d3


def plot_RNA_RPF_fc(dinput, figPos):
	
	### set annotated selenoproteins and chaperones
	chapList = ['HSPA5', 'HSPA1A', 'HSPA8', 'HSPH1', 'HSP90B1', 'DNAJA1', 'CALR', 'ERP44', 
	            'PDIA3', 'PDIA4', 'PDIA6', 'PFDN5', 'FKBP4', 'ST13', 'STIP1', 'HSPA1B']    

	selenList = ['GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX6', 'DIO1', 'DIO2', 'DIO3', 'SEPHS2',
	            'SELENOH', 'SELENOI', 'SELENOK', 'SELENOM', 'SELENON', 'SELENOO', 
	            'SELENOP', 'SELENOF', 'SELENOS', 'SELENOT', 'SELENOV', 'SELENOW',
	            'MSRB1', 'TXNRD1', 'TXNRD2', 'TXNRD3']

	dfp = dinput.copy()

	xVar = 'mRNA_log2FC'
	yVar = 'RPF_log2FC'

	g = sns.JointGrid(x=xVar, y=yVar, data=dfp, xlim=(-1.5,1.5), ylim=(-1.5,1.5), height = 6)
	g = g.plot_marginals(sns.distplot, kde=False)
	plt.setp(g.ax_marg_x.patches, color='black')
	plt.setp(g.ax_marg_y.patches, color='black')
	g = g.plot_joint(plt.scatter, color = "black", s=1)

	ax = plt.gca()

	for row in dfp.index:
	    geneName = dfp.loc[row, "gene_name"]
	    if "HIST" in geneName:
	        ax.scatter(x=dfp.loc[row,xVar],
	                       y=dfp.loc[row,yVar],
	                       color = colorDict['red'],
	                       edgecolor = 'black',
	                       linewidth = 0.5,
	                       s = 30)
	    if geneName in selenList:
	        ax.scatter(x=dfp.loc[row,xVar],
	                       y=dfp.loc[row,yVar],
	                       color = colorDict['green'],
	                       edgecolor = 'black',
	                       linewidth = 0.5,
	                       s = 30)
	    if geneName  in chapList:
	        ax.scatter(x=dfp.loc[row,xVar],
	                       y=dfp.loc[row,yVar],
	                       color= colorDict['pink'],
	                       edgecolor = 'black',
	                       linewidth = 0.5,
	                       s = 30)
	    if geneName == "AMD1":
	        ax.scatter(x=dfp.loc[row,xVar],
	                       y=dfp.loc[row,yVar],
	                       color = colorDict['blue'],
	                       edgecolor = 'black',
	                       linewidth = 0.5,
	                       s = 30)
	    if geneName == "ATF4":
	        ax.scatter(x=dfp.loc[row,xVar],
	                       y=dfp.loc[row,yVar],
	                       color = colorDict['orange'],
	                       edgecolor = 'black',
	                       linewidth = 0.5,
	                       s = 30)

	ax.axhline(y=0, linestyle="--", color='gray') 
	ax.axvline(x=0, linestyle="--", color='gray')       

	outfile = "%s/figures/Fig%s.pdf" % (rootDir, figPos)
	plt.savefig(outfile, format='pdf', bbox_inches="tight")

def main():
	d1, d2, d3, = load_Xtail_results()
	plot_RNA_RPF_fc(d1, "7B")
	plot_RNA_RPF_fc(d2, "7A")
	plot_RNA_RPF_fc(d3, "7S1A")



if __name__ == "__main__":
	main()






