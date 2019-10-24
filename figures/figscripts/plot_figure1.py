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
import seaborn as sns
from pylab import *
import argparse

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/figures" % args.rootDir)
sys.path.append("%s/figures/figscripts" % args.rootDir)


rootDir = args.rootDir
threadNumb = str(args.threadNumb)


#### functions:

def input_luc_data():
	df = pd.read_csv("%s/Data/Luc/Figure1_LucData.csv" % rootDir)
	exclusion_list = ['Untrans', 'H2O', 'Tris']
	df = df.loc[~df['Treatment'].isin(exclusion_list)]
	df['Ratio'] = df['Nluc']/df['Fluc']


	for i in df.index:
		df.loc[i,'Test'] = str(df.loc[i,'Treatment']) + "_" + str(df.loc[i,'Conc'])

	df.reset_index(drop=True, inplace=True)
	return df



def calc_luc_ratios(df):
	### for data collected in triplicate:

	index_all = range(len(df))

	# repN = 6
	tr_index = [0, 3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81]


	# tr_index = index_all[::repN]

	treatments = []
	means = []

	for i in tr_index:
		tr = str(df.loc[i,'Treatment'])
		conc = str(df.loc[i,'Conc'])
		
		trOut = tr+"_"+conc
		treatments.append(trOut)

	print treatments    

	fM = []
	fSD = []
	fSEM = []

	nM = []
	nSD = []
	nSEM = []

	rM =[]
	rSD =[]
	rSEM =[]


	for i in tr_index:
		
		if i == 0:
			repN = 3
		else:
			repN = 6

		fluc_vals = df.loc[i:i+repN-1,'Fluc']
		nluc_vals = df.loc[i:i+repN-1,'Nluc']
		ratio_vals = df.loc[i:i+repN-1, 'Ratio']
		
		
		### Fluc
		fluc_mean = np.mean(fluc_vals)
		fluc_sd = np.std(fluc_vals, ddof=1)
		fluc_sem = stats.sem(fluc_vals)
		
		fM.append(fluc_mean)
		fSD.append(fluc_sd)
		fSEM.append(fluc_sem)

		
		### Nluc
		nluc_mean = np.mean(nluc_vals)
		nluc_sd = np.std(nluc_vals, ddof=1)
		nluc_sem = stats.sem(nluc_vals)
	   
		nM.append(nluc_mean)
		nSD.append(nluc_sd)
		nSEM.append(nluc_sem)
		
		### ratio 
		ratio_mean = np.mean(ratio_vals)
		ratio_sd = np.std(ratio_vals, ddof=1)
		ratio_sem = stats.sem(ratio_vals)
		
		rM.append(ratio_mean)
		rSD.append(ratio_sd)
		rSEM.append(ratio_sem)

	new_ind = range(len(df)/repN)    

	df_dict = {'Test':treatments,
			   'fluc_mean':fM,
			   'fluc_sd':fSD,
			   'fluc_sem':fSEM,
			   'nluc_mean':nM,
			   'nluc_sd':nSD,
			   'nluc_sem':nSEM,
			   'ratio_mean':rM,
			   'ratio_sd':rSD,
			   'ratio_sem':rSEM
			  }

	dfout = pd.DataFrame.from_dict(df_dict)

	fluc_ctrl_mean = dfout.loc[0, 'fluc_mean']
	nluc_ctrl_mean = dfout.loc[0, 'nluc_mean']

	ctrl_mean = dfout.loc[0,'ratio_mean']
	# ctrl_sd = dfout.loc[1,'ratio_sd']
	# ctrl_sem = dfout.loc[1,'ratio_sem']

	dfout['rM_FC'] = dfout['ratio_mean']/ctrl_mean
	dfout['rSD_FC'] = dfout['ratio_sd']/ctrl_mean
	dfout['rSEM_FC'] = dfout['ratio_sem']/ctrl_mean


	return dfout, fluc_ctrl_mean, ctrl_mean



def plot_1B(df, dfout, fluc_ctrl_mean):
	### plot Fluc

	fig, ax = plt.subplots(figsize=(6,2.5))

	sns.swarmplot(data=df, x='Test', y='Fluc', size=3)
	sns.boxplot(data=dfout, x='Test', y='fluc_mean', showbox=False , width = 0.5, 
		showcaps=False, color = "black")
	plt.errorbar(dfout['Test'], dfout['fluc_mean'], yerr=dfout['fluc_sd'], capsize = 3, ls='none', color = "black")
	ax=plt.gca()

	for item in ax.get_xticklabels():
		item.set_rotation(90)

	axhline(y=fluc_ctrl_mean/2, linestyle="--", color="#999999")
	ax.set_ylim(0, 1000000)

	plt.savefig("%s/figures/Fig1B.pdf" % rootDir, format="pdf", bbox_inches = "tight")


def plot_1C(df, dfout):
	### plot Nluc
	# sns.catplot(x='Test', y='Nluc', data=df, kind='swarm', height = 6, aspect=12/6)

	fig, ax = plt.subplots(figsize=(6,2.5))

	sns.swarmplot(data=df, x='Test', y='Nluc', size=3)
	sns.boxplot(data=dfout, x='Test', y='nluc_mean', showbox=False , width = 0.5, 
		showcaps=False, color = "black")
	plt.errorbar(dfout['Test'], dfout['nluc_mean'], yerr=dfout['nluc_sd'], capsize = 3, ls='none', color = "black")
	ax=plt.gca()

	for item in ax.get_xticklabels():
		item.set_rotation(90)


	ax.set_ylim(0, 1000000)

	plt.savefig("%s/figures/Fig1C.pdf" % rootDir, format="pdf", bbox_inches = "tight")


def plot_1D_low(df, dfout, ctrl_mean):

	dfp1 = df.copy()
	dfp2 = dfout.copy()

	dfp1['RatioScaled'] = dfp1['Ratio']/ctrl_mean

	fig, ax = plt.subplots(figsize=(6,1.25))

	sns.swarmplot(data=dfp1, x='Test', y='RatioScaled', size=3)
	sns.boxplot(data=dfp2, x='Test', y='rM_FC', showbox=False , width = 0.5, 
					 showcaps=False, color = "black")
	plt.errorbar(dfp2['Test'], dfout['rM_FC'], yerr=dfout['rSD_FC'], capsize = 3, ls='none', color = "black")
	ax=plt.gca()

	for item in ax.get_xticklabels():
		item.set_rotation(90)

	ax.set_ylim(0, 40)
	ax.set_yticks([0, 10, 20, 30, 40])

	plt.savefig("%s/figures/Fig1D_lower.pdf" % rootDir, format="pdf", bbox_inches = "tight")



def plot_1D_high(df, dfout, ctrl_mean):
	### plot Ratio

	dfp1 = df.copy()
	dfp2 = dfout.copy()

	dfp1['RatioScaled'] = dfp1['Ratio']/ctrl_mean

	fig, ax = plt.subplots(figsize=(6,1.25))

	sns.swarmplot(data=dfp1, x='Test', y='RatioScaled', size=3)
	sns.boxplot(data=dfp2, x='Test', y='rM_FC', showbox=False , width = 0.5, 
					 showcaps=False, color = "black")
	plt.errorbar(dfp2['Test'], dfout['rM_FC'], yerr=dfout['rSD_FC'], capsize = 3, ls='none', color = "black")
	ax=plt.gca()

	for item in ax.get_xticklabels():
		item.set_rotation(90)

	ax.set_ylim(50, 1000)
	ax.set_yticks([50, 250, 500, 750, 1000])

	plt.savefig("%s/figures/Fig1D_high.pdf" % rootDir, format="pdf", bbox_inches = "tight")



def main():
	df = input_luc_data()
	dfout, fluc_ctrl_mean, ctrl_mean = calc_luc_ratios(df)
	plot_1B(df, dfout, fluc_ctrl_mean)
	plot_1C(df, dfout)
	plot_1D_low(df, dfout, ctrl_mean)
	plot_1D_high(df, dfout, ctrl_mean)

if __name__ == '__main__':
	main()

