### generate plots of luciferase data:


### Import dependencies
import matplotlib
matplotlib.use('Agg') ### set backend
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines
plt.tight_layout()

import sys
import math
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats
# sys.path.append('/home/jamie/.local/lib/python2.7/site-packages')
import seaborn as sns
from pylab import *
import argparse
import importlib

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/figures" % args.rootDir)
sys.path.append("%s/figures/figscripts" % args.rootDir)


rootDir = args.rootDir
threadNumb = str(args.threadNumb)

### set custom colors
colorList = ['#000000', '#ffb000', '#63cfff', '#eb4300', '#00c48f', '#eb68c0', '#fff71c', '#006eb9']
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

#### functions:

def input_luc_data():
	df = pd.read_csv("/home/jamie/G418paper/RawData/LucAssays/JRW-R-41_Luc_tidy.csv")


	df['Ratio_Untr'] = df['Nluc_Untr']/df['Fluc_Untr']
	df['Ratio_G418'] = df['Nluc_G418']/df['Fluc_G418']

	for i in df.index:
	    df.loc[i,'StopCodon'] = df.loc[i,'PTC']+"_"+df.loc[i,'Codon']

	return df

def calc_luc_ratios(df):

	index_all = range(len(df))
	tr_index = range(0, len(df), 3)
	repN = 3

	treatments = []

	for i in tr_index:
	    
	    trOut = df.loc[i,'StopCodon']
	    treatments.append(trOut)

	fM_Untr = []
	fSD_Untr = []
	fSEM_Untr = []

	fM_G418 = []
	fSD_G418 = []
	fSEM_G418 = []

	nM_Untr = []
	nSD_Untr = []
	nSEM_Untr = []

	nM_G418 = []
	nSD_G418 = []
	nSEM_G418 = []

	rM_Untr =[]
	rSD_Untr =[]
	rSEM_Untr =[]

	rM_G418 =[]
	rSD_G418 =[]
	rSEM_G418 =[]

	for i in tr_index:
	    
	    fluc_vals_Untr = df.loc[i:i+repN-1,'Fluc_Untr']
	    nluc_vals_Untr = df.loc[i:i+repN-1,'Nluc_Untr']
	    ratio_vals_Untr = df.loc[i:i+repN-1, 'Ratio_Untr']
	    
	    fluc_vals_G418 = df.loc[i:i+repN-1,'Fluc_G418']
	    nluc_vals_G418 = df.loc[i:i+repN-1,'Nluc_G418']
	    ratio_vals_G418 = df.loc[i:i+repN-1, 'Ratio_G418']
	    
	    
	    ### Fluc
	    fluc_mean_Untr = np.mean(fluc_vals_Untr)
	    fluc_sd_Untr = np.std(fluc_vals_Untr, ddof=1)
	    fluc_sem_Untr = stats.sem(fluc_vals_Untr)
	    
	    fluc_mean_G418 = np.mean(fluc_vals_G418)
	    fluc_sd_G418 = np.std(fluc_vals_G418, ddof=1)
	    fluc_sem_G418 = stats.sem(fluc_vals_G418)
	    
	    fM_Untr.append(fluc_mean_Untr)
	    fSD_Untr.append(fluc_sd_Untr)
	    fSEM_Untr.append(fluc_sem_Untr)
	    
	    fM_G418.append(fluc_mean_G418)
	    fSD_G418.append(fluc_sd_G418)
	    fSEM_G418.append(fluc_sem_G418)

	    
	    ### Nluc
	    nluc_mean_Untr = np.mean(nluc_vals_Untr)
	    nluc_sd_Untr = np.std(nluc_vals_Untr, ddof=1)
	    nluc_sem_Untr = stats.sem(nluc_vals_Untr)
	    
	    nluc_mean_G418 = np.mean(nluc_vals_G418)
	    nluc_sd_G418 = np.std(nluc_vals_G418, ddof=1)
	    nluc_sem_G418 = stats.sem(nluc_vals_G418)
	    
	    nM_Untr.append(nluc_mean_Untr)
	    nSD_Untr.append(nluc_sd_Untr)
	    nSEM_Untr.append(nluc_sem_Untr)
	    
	    nM_G418.append(nluc_mean_G418)
	    nSD_G418.append(nluc_sd_G418)
	    nSEM_G418.append(nluc_sem_G418)
	    
	    
	    ### ratio 
	    ratio_mean_Untr = np.mean(ratio_vals_Untr)
	    ratio_sd_Untr = np.std(ratio_vals_Untr, ddof=1)
	    ratio_sem_Untr = stats.sem(ratio_vals_Untr)
	    
	    ratio_mean_G418 = np.mean(ratio_vals_G418)
	    ratio_sd_G418 = np.std(ratio_vals_G418, ddof=1)
	    ratio_sem_G418 = stats.sem(ratio_vals_G418)
	    
	    rM_Untr.append(ratio_mean_Untr)
	    rSD_Untr.append(ratio_sd_Untr)
	    rSEM_Untr.append(ratio_sem_Untr)
	    
	    rM_G418.append(ratio_mean_G418)
	    rSD_G418.append(ratio_sd_G418)
	    rSEM_G418.append(ratio_sem_G418)

	new_ind = range(len(df)/repN) 


	df_dict = {'StopCodon':treatments,
	           'fluc_mean_Untr':fM_Untr,
	           'fluc_sd_Untr':fSD_Untr,
	           'fluc_sem_Untr':fSEM_Untr,
	           'nluc_mean_Untr':nM_Untr,
	           'nluc_sd_Untr':nSD_Untr,
	           'nluc_sem_Untr':nSEM_Untr,
	           'ratio_mean_Untr':rM_Untr,
	           'ratio_sd_Untr':rSD_Untr,
	           'ratio_sem_Untr':rSEM_Untr,
	           'fluc_mean_G418':fM_G418,
	           'fluc_sd_G418':fSD_G418,
	           'fluc_sem_G418':fSEM_G418,
	           'nluc_mean_G418':nM_G418,
	           'nluc_sd_G418':nSD_G418,
	           'nluc_sem_G418':nSEM_G418,
	           'ratio_mean_G418':rM_G418,
	           'ratio_sd_G418':rSD_G418,
	           'ratio_sem_G418':rSEM_G418
	          }

	dfout = pd.DataFrame.from_dict(df_dict)


	### normalization

	minVal = dfout['ratio_mean_Untr'].min()

	dfout['ratio_mean_Untr_Norm'] = dfout['ratio_mean_Untr']*(1/minVal)
	dfout['ratio_mean_G418_Norm'] = dfout['ratio_mean_G418']*(1/minVal)
	dfout['ratio_sd_Untr_Norm'] = dfout['ratio_sd_Untr']*(1/minVal)
	dfout['ratio_sd_G418_Norm'] = dfout['ratio_sd_G418']*(1/minVal)

	return dfout, minVal

def adjust_ratio(df, minVal):

	df['Ratio_Untr_Norm'] = df['Ratio_Untr']*(1/minVal)
	df['Ratio_G418_Norm'] = df['Ratio_G418']*(1/minVal)

	return df

def plot_luc_data(df, dfout):
	fig, ax = plt.subplots(figsize=(6,3))

	xVar = 'StopCodon'
	yVar = 'Ratio_Untr_Norm'
	yVar2 = 'ratio_mean_Untr_Norm'
	yErr = 'ratio_sd_Untr_Norm'

	sns.swarmplot(data=df, x=xVar, y=yVar, size=5, color= colorDict['black'], edgecolor='black', linewidth=0.5)
	plt.errorbar(dfout[xVar], dfout[yVar2], yerr=dfout[yErr], capsize = 3, ls='none', color = "gray")
	plt.bar(x=dfout[xVar], height=dfout[yVar2], edgecolor='gray', fill=False)
	
	ax=plt.gca()
	for item in ax.get_xticklabels():
	    item.set_rotation(90)

	ax.set_ylim(0, 110)

	plt.savefig("%s/figures/Fig4S1A.pdf" % rootDir, format="pdf", bbox_inches = "tight")



def main():
	df = input_luc_data()
	dfout, minVal = calc_luc_ratios(df)
	df = adjust_ratio(df, minVal)
	plot_luc_data(df, dfout)

if __name__ == '__main__':
	main()


