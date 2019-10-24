

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
parser.add_argument('--libSetFile', help = 'libset file')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)
import rphelper as rph

ootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)
### set colors 

# outfileOrder = ['7A_right', '7A_left', '7S1A']



### set colors
# black = '#000000'
# orange = '#ffb000'
# cyan = '#63cfff'
# red = '#eb4300'
# green = '#00c48f'
# pink = '#eb68c0'
# yellow = '#fff71c'
# blue = '#006eb9'

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

def load_samples():
	### 
	ThreshRAW = 128

	FPassignpath = "%s/FPassignment/%s/%s" % (rootpath, genome_name, experiment)
	namelist = []
	dflist= []
	# dflu3 = []
	for samp in samplelist:
	#     dftemp = pd.read_csv('%s/%s/countTables/%s_fl_rpm_28to35_countTable_rpkm.csv' % (FPassignpath, samp, samp))
	#     dftemp['sampname'] = samp
	# #     print dftemp.head()
	#     dflist.append(dftemp)
		
		###
		
		dftemp = pd.read_csv('%s/%s/countTables/%s_fl_rpm_28to35_countTable_rpkm_utr3adj.csv' % (FPassignpath, samp, samp))
		dftemp['sampname'] = samp
		dftemp['cdsCountsLog2'] = dftemp['cdsCounts'].apply(log_trans_b2)
		dftemp['RAWcdsCountsLog2'] = dftemp['RAWcdsCounts'].apply(log_trans_b2)
		dftemp['utr3CountsLog2'] = dftemp['utr3Counts'].apply(log_trans_b2)
		dftemp['utr3AdjCountsLog2'] = dftemp['utr3AdjCounts'].apply(log_trans_b2)
		dftemp['utr3OccLog2'] = dftemp['utr3_occupancy'].apply(log_trans_b2)
		dftemp['RRTSlog2'] = dftemp['RRTS'].apply(log_trans_b2)
		dftemp['cdsDenstilyRPKMlog2'] = dftemp['cdsDensity_rpkm'].apply(log_trans_b2)
		dftemp['cdslenLog2'] = dftemp['cdslen'].apply(log_trans_b2)
		dftemp['utr3LenAdjLog2'] = dftemp['utr3LenAdj'].apply(log_trans_b2)
	#     print dftemp.head()

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > ThreshRAW]
		print samp
		print len(dftemp)

		dflist.append(dftemp)
		
		
		namelist.append(samp)
		
	# ### combine into one master dataframe
	# df = pd.concat(dflist, axis=0, ignore_index = True)
	df = pd.concat(dflist, axis=0, ignore_index = True)

	return df, dflist, namelist

def plot_corr(dflist):

	samp_ind = [18, 20]

	treatList = ['Untr', 'G418_500']
	colorlist = [colorDict['black'],
				 colorDict['orange'],
				 colorDict['red'],
				 colorDict['green'],
				 colorDict['blue'],
				 colorDict['yellow'],
				 colorDict['pink'],
				 colorDict['cyan'],
				 colorDict['orange']]


	counter = 0
	for i in samp_ind:
		s1 = i
		s2 = i+1

		d1 = dflist[s1]
		d2 = dflist[s2]

		treat = treatList[counter]
		d = d1.merge(d2, on="#transcript")

		g = sns.JointGrid(x='RRTSlog2_x', y='RRTSlog2_y', data=d, xlim=(-16,4), ylim=(-16,4))
		g = g.plot_joint(plt.scatter, color = "black", s=1)
		g = g.plot_marginals(sns.distplot, kde=False, bins=20)# colorlist[counter])
		plt.setp(g.ax_marg_x.patches, alpha=0.8, color=colorlist[counter], edgecolor="black")
		plt.setp(g.ax_marg_y.patches, alpha=0.8, color=colorlist[counter], edgecolor="black")

		figout = "%s/figures/Fig3S1A_%s_%s.pdf" % (rootDir, counter, treat)
		plt.savefig(figout, format='pdf', bbox_inches="tight")
		counter +=1

def main():
	df, dflist, namelist = load_samples()
	plot_corr(dflist)

if __name__ == '__main__':
	main()




















