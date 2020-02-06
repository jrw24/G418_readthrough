"""
This script plots codon occupancies relative to untreated cells

"""
## plotting
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
import pandas as pd
import pysam
import numpy as np
import math
import argparse
import importlib
from scipy import stats

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)


rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)

threadNumb = str(args.threadNumb)
import rphelper as rph


### inputs:

readlength = "28to35"
shift = "A"

## indexes
ctrl_1_index = 0
ctrl_2_index = 1
treat_index_list = [2,4,6,8,10,12,14,16]
treatList = ['G418-500', 'G418-2000', 'Genta', 'Parom', 'Neo', 'Tobra', 'Amik', 'G418-500-10min']

## colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'

colorDict = {
	'D':cyan,
	'G':orange,
	'I':green
}

### AA to codon dict
codonDict = {
	'A':['GCT', 'GCC', 'GCA', 'GCG'],
	'C':['TGT', 'TGC'],
	'D':['GAT', 'GAC'],
	'E':['GAA', 'GAG'],
	'F':['TTT', 'TTC'],
	'G':['GGT', 'GGC', 'GGA', 'GGG'],
	'H':['CAT', 'CAC'],
	'I':['ATT', 'ATC', 'ATA'],
	'K':['AAA', 'AAG'],
	'L':['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
	'M':['ATG'],
	'N':['AAT', 'AAC'],
	'P':['CCT', 'CCC', 'CCA', 'CCG'],
	'Q':['CAA', 'CAG'],
	'R':['AGA','AGG','CGT','CGC','CGA','CGG'],
	'S':['TCT', 'TCC', 'TCA', 'TCG'],
	'T':['ACT', 'ACC', 'ACA', 'ACG'],
	'V':['GTT', 'GTC', 'GTA', 'GTG'],
	'W':['TGG'],
	'Y':['TAT', 'TAC']
}




def log_trans_b10(x):
	try:
		return math.log(x, 10)
	except:
		return float(-6.00)
#         return float("NaN")
def log_trans_b2(x):
	try:
		return math.log(x, 2)
	except:
#         return float("NaN")
		return float(-15.00) # set arbitrarily low value

def corrfunc(x, y, **kws):
	r, _ = stats.pearsonr(x, y)
	rho, _ = stats.spearmanr(x, y)
	ax = plt.gca()
	ax.annotate("R = {:.3f}".format(r),
				xy=(.1, .9), xycoords=ax.transAxes)






def plot_codons():

	### AG4 datasets
	
	counter = 0
	for i in treat_index_list:

		samp1 = samplelist[ctrl_1_index]
		samp2 = samplelist[ctrl_2_index]

		samp3 = samplelist[i]
		samp4 = samplelist[i+1]

		codon_dir = "%s/FPassignment/%s/%s" % (rootpath, genome_name, experiment)

		df1 = pd.read_csv("%s/%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % 
						 (codon_dir, samp1, samp1, shift, readlength), index_col=0)
		df2 = pd.read_csv("%s/%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % 
						 (codon_dir, samp2, samp2, shift, readlength), index_col=0)
		df3 = pd.read_csv("%s/%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % 
						 (codon_dir, samp3, samp3, shift, readlength), index_col=0)
		df4 = pd.read_csv("%s/%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % 
						 (codon_dir, samp4, samp4, shift, readlength), index_col=0)

		dflist = [df1, df2, df3, df4]

		dfout = pd.concat(dflist, axis=1)
		dfout.drop(index=['TAA', 'TAG', 'TGA'],axis=0,inplace=True)
		dfout = dfout[[samp1,samp2,samp3,samp4]].apply(pd.to_numeric)

		for col in dfout.columns:
			dfout[col+'_log2'] = dfout[col].apply(log_trans_b2)


		dfout['ctrl_mean'] = dfout[[samp1+'_log2',samp2+'_log2']].mean(axis=1)
		dfout['ctrl_std'] = dfout[[samp1+'_log2',samp2+'_log2']].std(axis=1)
		dfout['ctrl_sem'] = dfout[[samp1+'_log2',samp2+'_log2']].sem(axis=1)

		dfout['treat_mean'] = dfout[[samp3+'_log2',samp4+'_log2']].mean(axis=1)
		dfout['treat_std'] = dfout[[samp3+'_log2',samp4+'_log2']].std(axis=1)
		dfout['treat_sem'] = dfout[[samp3+'_log2',samp4+'_log2']].sem(axis=1)



		### plotting:
		fig, ax = plt.subplots(figsize=(6,6))

		for aa in codonDict:
			for cdn in codonDict[aa]:

				plt.scatter(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							color="black", 
							edgecolor="black", 
							linewidth = 0.5, 
							s=20
						   )
				plt.errorbar(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							xerr = dfout.loc[cdn,'ctrl_std'],
							yerr = dfout.loc[cdn,'treat_std'],
							color="black", 
							ecolor="black", 
							linewidth = 0.75
							)	

		for aa in colorDict:
			for cdn in codonDict[aa]:
	#           
				plt.scatter(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							color=colorDict[aa], 
							edgecolor="black", 
							linewidth = 0.5, 
							s=20
						   )

				plt.text(x=dfout.loc[cdn,'ctrl_mean']+0.00, 
						 y=dfout.loc[cdn,'treat_mean']+0.00, 
						 color=colorDict[aa], 
						 s=aa, 
						 fontsize = 16)
				plt.text(x=dfout.loc[cdn,'ctrl_mean']+0.00, 
						 y=dfout.loc[cdn,'treat_mean']-0.00, 
						 color=colorDict[aa], #"black", 
						 s=cdn, 
						 fontsize = 10, 
						 va="top", 
						 ha="left")

		plt.title('Untr'+' vs '+treatList[counter], color='gray')

		lims = [
			np.min([-1.0, -1.0]),
			np.max([1.5, 1.5]),
		]

		ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
		ax.set_aspect('equal')
		ax.set_xlim(lims)
		ax.set_ylim(lims)
		corrfunc(dfout['ctrl_mean'], dfout['treat_mean'])

		outfile = "%s/figures/Fig2S3B_Untr_vs_%s.pdf" % (rootDir, treatList[counter])
		plt.savefig(outfile, format="pdf")

		counter +=1


def main():
	plot_codons()

if __name__ == '__main__':
	main()






