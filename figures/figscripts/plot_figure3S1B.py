


"""
Script for plotting Venn Diagrams comparing total number of transcripts with detectable readthrough

"""


import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

from matplotlib_venn import venn2
import sys, os
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
import csv
from datetime import datetime
import argparse
from pathos.multiprocessing import ProcessingPool as Pool
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

rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)


samplelist = samplelist[18:22]



### set colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'

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

def plot_venn(dflist):
	"""
	create ven diagrams for total number of transcripts with detectable readthorugh
	"""

	### Untreated
	figout_untr = "%s/figures/Fig3S1B_untr.pdf" % (rootDir)

	untr_1 = dflist[0]
	untr_2 = dflist[1]

	print "plotting Venn Diagrams of %s vs %s" % (samplelist[0], samplelist[1])

	df_untr = untr_1.merge(untr_2, on="#transcript")
	u1v = df_untr.loc[df_untr['RRTS_x']>0]
	u2v = df_untr.loc[df_untr['RRTS_y']>0]

	fig, ax = plt.subplots(figsize=(4,4))
	v = venn2 ([set(u1v['#transcript']), set(u2v['#transcript'])])
	plt.savefig(figout_untr, format='pdf', bbox_inches="tight")


	### G418
	figout_g418 = "%s/figures/Fig3S1B_g418.pdf" % (rootDir)

	g418_1 = dflist[2]
	g418_2 = dflist[3]

	print "plotting Venn Diagrams of %s vs %s" % (samplelist[2], samplelist[3])

	df_g418 = g418_1.merge(g418_2, on="#transcript")
	g1v = df_g418.loc[df_g418['RRTS_x']>0]
	g2v = df_g418.loc[df_g418['RRTS_y']>0]

	fig, ax = plt.subplots(figsize=(4,4))
	v = venn2 ([set(g1v['#transcript']), set(g2v['#transcript'])])
	plt.savefig(figout_g418, format='pdf', bbox_inches="tight")



def main():
	df, dflist, namelist = load_samples()
	plot_venn(dflist)

if __name__ == '__main__':
	main()
