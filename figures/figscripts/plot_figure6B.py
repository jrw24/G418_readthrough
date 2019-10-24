#####


import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import math
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import importlib
pd.set_option('display.max_columns', 40)


### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/RNAseq" % args.rootDir)
sys.path.append("%s/RNAseq/libsettings" % args.rootDir)


rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)
sys.path.append('%s/riboseq' % rootDir)


def load_samples():

	FPassignpath = "%s/FPassignment/%s/%s" % (rootpath, genome_name, experiment)
	namelist = []
	dflist= []

	for samp in samplelist:

	    dftemp = pd.read_csv('%s/%s/countTables/%s_5p_rnaseq_rpm_countTable.csv' % (FPassignpath, samp, samp))
	    dftemp['sampname'] = samp
	    dflist.append(dftemp)
	    namelist.append(samp)

	df = pd.concat(dflist, axis=0, ignore_index = True)

	return df, dflist

def plot_rpkm(df):
	dfp = df.loc[df['#transcript']== 'ENST00000445888.6'] ## TP53
	dfp = dfp.iloc[0:4,]
	dfp['mrnaRPKM'] = dfp['mrnaDensity'] * 1000
	dfp['treat'] = ['untr', 'untr', 'g418', 'g418']
	outfile = "%s/figures/Fig6B.pdf" % rootDir

	sns.catplot(data=dfp, x='treat', y='mrnaRPKM', hue='treat', kind='swarm', height = 3, aspect = 0.75)
	sns.boxplot(data=dfp, x='treat', y='mrnaRPKM', showbox=False , width = 0.25, showcaps=False,
	           whiskerprops= dict(color="white"))
	plt.savefig(outfile, format='pdf', bbox_inches="tight")

def main():

	df, dflist = load_samples()
	plot_rpkm(df)

if __name__ == '__main__':
	main()
