"""
script is based off of riboseq_readingFrame_utr3meta
"""
## plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
import math
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scipy.stats as stats
pd.set_option('display.max_columns', 50)
import seaborn as sns
from pylab import *
import argparse
import importlib

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
sys.path.append('%s/riboseq' % rootDir)


### set colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'


# colorList = [black, orange, red, green, blue, yellow, pink]
colorList = ['lightgray', 'gray', 'black', 'gold', 'darkorange','orangred']

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
		# return float("NaN")
		return float(-15.00) # set arbitrarily low value


def load_countTables():
	FPassignpath = "%s/FPassignment/%s/%s" % (rootpath, genome_name, experiment)
	namelist = []
	dflist= []
	# dflu3 = []
	for samp in samplelist:
	
		###
		dftemp = pd.read_csv('%s/%s/countTables/%s_fl_rpm_28to35_countTable_rpkm_utr3adj.csv' % (FPassignpath, samp, samp))
		dftemp['sampname'] = samp
		dftemp['cdsCountsLog2'] = dftemp['cdsCounts'].apply(log_trans_b2)
		dftemp['RAWcdsCountsLog2'] = dftemp['RAWcdsCounts'].apply(log_trans_b2)
		dftemp['utr3CountsLog2'] = dftemp['utr3Counts'].apply(log_trans_b2)
		dftemp['utr3AdjCountsLog2'] = dftemp['utr3AdjCounts'].apply(log_trans_b2)
		dftemp['utr3OccLog2'] = dftemp['utr3_occupancy'].apply(log_trans_b2)
		dftemp['RRTSlog2'] = dftemp['RRTS'].apply(log_trans_b2)

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > 128]
		dflist.append(dftemp)
		namelist.append(samp)
		
	# ### combine into one master dataframe
	df = pd.concat(dflist, axis=0, ignore_index = True)
	return df, dflist, namelist

def plot_RRTS_boxplot_4nt(dflist, namelist):

	# sampInd = 20
	# dfp = dflist[sampInd]

	dfp = dflist[20]

	figout = '%s/figures/Fig4D.pdf' % rootDir
	dfp = dfp.sort_values("stop4nt", ascending=True)

	scO = dfp['stopcodon'].unique()
	stopCodeOrd = dfp['stop4nt'].unique()

	scl = []
	for i in stopCodeOrd:
	    scl.append(i)

	fig, ax = subplots(figsize=(12, 8))
	ax = sns.boxplot(x=dfp['stop4nt'],  y=dfp['RRTS'], 
	                 orient = 'v', fliersize = 0, order=scl, notch=True, width=0.8)

	ax.set_ylim((0,1.4))
	plt.savefig(figout, format='pdf', bbox_inches="tight")

def main():
	df, dflist, namelist = load_countTables()
	plot_RRTS_boxplot_4nt(dflist, namelist)

if __name__ == '__main__':
	main()

