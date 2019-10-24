



import matplotlib
matplotlib.use('TkAgg') ## set backend here
import logomaker as lm
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
pd.set_option('display.max_columns', 100)
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
import statsmodels.stats.multitest as mtest
import argparse
import importlib

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from statsmodels.regression.linear_model import OLS

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


mrnadf = pd.read_csv("%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_mRNAseqs.csv" % rootDir)
scdf = pd.read_csv(stopcodons)


upstreamNTs = 15
downstreamNTs = 12

samplelist = [samplelist[18]]+[samplelist[20]]
treatlist = ['Untr', 'G418']

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

def build_master_df():
	### 
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
		dftemp['cdslenLog2'] = dftemp['cdslen'].apply(log_trans_b2)
		dftemp['utr3LenAdjLog2'] = dftemp['utr3LenAdj'].apply(log_trans_b2)
	#     print dftemp.head()

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > 128]
		print len(dftemp)

		dflist.append(dftemp)
		
		
		namelist.append(samp)
		
	# ### combine into one master dataframe
	# df = pd.concat(dflist, axis=0, ignore_index = True)
	df = pd.concat(dflist, axis=0, ignore_index = True)

	return df, dflist


# def choose_single_sample(dflist, samplenumb):

# 	dfsingle = dflist[samplenumb]
# 	print len(dfsingle)
# 	dfsingle = dfsingle.merge(mrnadf, on = '#transcript')
# 	print len(dfsingle)
# 	dfsingle = dfsingle.merge(scdf, on='#transcript')
# 	dfsingle = find_stop_context(dfsingle)

# 	return dfsingle


def find_stop_context(d, upstreamNTs, downstreamNTs):
	
	for i in range(len(d)):
		sequence = d.loc[i,'mRNAseqs']
		pos = d.loc[i,'5utr_len'] + d.loc[i,'cds_len'] 
		
		d.loc[i,'stopContext'] = sequence[pos-3-upstreamNTs:pos+downstreamNTs].upper()

	return d


def scale_RRTS(d, scalor):

	df = d.copy()
	df['RRTS_scaled'] = df['RRTS']*scalor
	df['RRTS_scaled'] = df['RRTS_scaled'].astype('int')

	print len(df)
	return df 


def winsorize_RRTS_scaled(d, winsor_thresh):

	d['RRTS_scaled'] = stats.mstats.winsorize(d['RRTS_scaled'], limits = (0, winsor_thresh))
	print "new winsor max == ", max(d['RRTS_scaled'])
	d.reset_index(drop=True, inplace=True)
	return d


def build_nt_matrix(df, upstreamNTs, downstreamNTs):
	
	seqz = df['stopContext']

	dfoutlist = []


	### first build dataframe to hold these data
	possible_nts = ['A', 'C', 'G', 'T']

	# sequence_positions = range(-30, 33)
	sequence_positions = range(-upstreamNTs, downstreamNTs+3)
	sequence_positions

	index_list = []

	for i in sequence_positions:
	#     print i
		for j in possible_nts:
			index_list.append(str(i)+"_"+j)
	#         print str(i)+"_"+j

	index_list

	dfreg = pd.DataFrame(columns=index_list)
	print dfreg

	geneCounter = 0
	for seq in seqz:
		### build a single dictionary for a given sequence
		
	#     print len(seq)
		if len(seq) < len(sequence_positions):
			print "too short"
			geneCounter+=1
			continue
		
		ntdict = OrderedDict()

	#     counter = -30
		counter = -upstreamNTs
		
		for c in seq:
			ntdict[counter] = c
			counter+=1

	#     print ntdict
			
		dftemplate = dfreg.copy()
		val_list = [0 for char in dfreg.columns]

		cnter = 0
		for i in dfreg.columns:
			query = i.split("_")
			qval = int(query[0])

			if ntdict[qval] == query[1]:
				val_list[cnter] += 1

			cnter += 1

			
		gene = df.loc[geneCounter, '#transcript']
		geneCounter +=1
		
		tempdf = pd.DataFrame([val_list], columns = dfreg.columns)
		
		tempdf.index = [gene]

		dfoutlist.append(tempdf)

	dfout= pd.concat(dfoutlist)
	return dfout

def build_nt_matrix_collapseStop(df, upstreamNTs, downstreamNTs):
	"""
	collapse the stop codon into a single position
	"""
	
	seqz = df['stopContext']

	dfoutlist = []


	### first build dataframe to hold these data
	possible_nts = ['A', 'C', 'G', 'T']
	possible_stops = ['TAA', 'TAG', 'TGA']
	stop_positions = [0, 1, 2]


	# sequence_positions = range(-30, 33)
	sequence_positions = range(-upstreamNTs, downstreamNTs+3)
	sequence_positions

	index_list = []

	for i in sequence_positions:
		if i == 1:
			continue
		elif i == 2:
			continue
		elif i in stop_positions:
			for j in possible_stops:
				index_list.append(str(i)+"_"+j)
		else:
			for j in possible_nts:
				index_list.append(str(i)+"_"+j)
	#         print str(i)+"_"+j

	index_list

	dfreg = pd.DataFrame(columns=index_list)
	print dfreg

	geneCounter = 0
	for seq in seqz:
		### build a single dictionary for a given sequence
		
	#     print len(seq)
		if len(seq) < len(sequence_positions)-2:
			print "too short"
			geneCounter+=1
			continue
		
		ntdict = OrderedDict()

	#     counter = -30
		counter = -upstreamNTs
		
		for c in seq:
#         print c
			if counter == 1:
				counter+=1
				continue
			elif counter == 2:
				counter +=1
				continue
			elif counter == 0:
				# stopCodon = seq[6:6+3]
				stopCodon = seq[upstreamNTs:upstreamNTs+3]
				ntdict[0] = stopCodon
				counter+=1
			else:
				ntdict[counter] = c
				counter+=1      

		dftemplate = dfreg.copy()
		val_list = [0 for char in dfreg.columns]

		cnter = 0
		for i in dfreg.columns:
			query = i.split("_")
			qval = int(query[0])

			if ntdict[qval] == query[1]:
				val_list[cnter] += 1

			cnter += 1

			
		gene = df.loc[geneCounter, '#transcript']
		geneCounter +=1
		
		tempdf = pd.DataFrame([val_list], columns = dfreg.columns)
		
		tempdf.index = [gene]

		dfoutlist.append(tempdf)

	dfout= pd.concat(dfoutlist)
	return dfout


def get_RRTS_vals(df, dfout):
	rtvals = df['RRTS']
	# rtvals = df['RRTS_scaled']
	rtvals.index = df['#transcript']
	rtvals = rtvals.loc[rtvals.index.isin(dfout.index)]
	return rtvals


def plot_Ridge_Reg(dfout, rtvals, alphA, samp):
	X = dfout.values 
	y = rtvals 

	clf = Ridge(alpha=alphA)
	clf.fit(X, y)

	outcoeffs = clf.coef_
	outcoeffs

	outcoef = pd.DataFrame(outcoeffs, index = dfout.columns)


	import matplotlib.lines as mlines

	fig, ax = plt.subplots(figsize=(6,6))

	ntCols = {"A": "blue",
			  "C" : "orange",
			  "G" : "green",
			  "T" : "red"}


	for pos in outcoef.index:
	#     print pos, coeffs_lm.loc[pos]
		
		ntpos = int(pos.split("_")[0])
		nt = pos.split("_")[1]
		
		val = outcoef.loc[pos]
		
		plt.plot(ntpos, val, marker='o', color = ntCols[nt], markersize = 6)
		
		
	patchA = mpatches.Patch(color="blue", label = "A")
	patchC = mpatches.Patch(color="orange", label = "C")
	patchG = mpatches.Patch(color="green", label = "G")
	patchT = mpatches.Patch(color="red", label = "T")

	plt.legend(handles = [patchA, patchC, patchG, patchT])

	plt.axhline(0, linestyle="--", color = "gray")

	plotoutfile = "/home/jamie/misc/rt/ridge/%s_RRTS_RidgeReg_%s_to_%s.pdf" % (samp, upstreamNTs, downstreamNTs)

	plt.savefig(plotoutfile, format='pdf')


def plot_Ridge_Reg_StopCollapse(dfout, rtvals, alphA, samp):

	X = dfout.values 
	y = rtvals 

	clf = Ridge(alpha=alphA)
	clf.fit(X, y)

	outcoeffs = clf.coef_
	outcoeffs

	outcoef = pd.DataFrame(outcoeffs, index = dfout.columns)

	print outcoef.head()


	import matplotlib.lines as mlines

	fig, ax = plt.subplots(figsize=(6,6))

	# ntCols = {"A": "blue",
	# 		  "C" : "orange",
	# 		  "G" : "green",
	# 		  "T" : "red",
	# 		  "TAA": "purple",
	# 		  "TAG": "brown",
	# 		  "TGA": "pink"}

	ntCols = {"A": sns.color_palette()[0],
		  "C" : sns.color_palette()[1],
		  "G" : sns.color_palette()[2],
		  "T" : sns.color_palette()[3],
		  "TAA": sns.color_palette()[4],
		  "TAG": sns.color_palette()[5],
		  "TGA": sns.color_palette()[6]}

	markerDict = {"A": "v",
		  "C" : "s",
		  "G" : "P",
		  "T" : "D",
		  "TAA": "o",
		  "TAG": "p",
		  "TGA": "h"}

	for pos in outcoef.index:
	#     print pos, coeffs_lm.loc[pos]
		
		ntpos = int(pos.split("_")[0])
		nt = pos.split("_")[1]
		
		if ntpos > 0:
			ntpos = ntpos-2
			
		val = outcoef.loc[pos]
		
		plt.plot(ntpos, val, marker=markerDict[nt], color = ntCols[nt], markersize = 6)
		
		
	# patchA = mpatches.Patch(color="blue", label = "A")
	# patchC = mpatches.Patch(color="orange", label = "C")
	# patchG = mpatches.Patch(color="green", label = "G")
	# patchT = mpatches.Patch(color="red", label = "T")
	# patchUAA = mpatches.Patch(color="purple", label = "UAA")
	# patchUAG = mpatches.Patch(color="brown", label = "UAG")
	# patchUGA = mpatches.Patch(color="pink", label = "UGA")

	patchA = mpatches.Patch(color=sns.color_palette()[0], label = "A")
	patchC = mpatches.Patch(color=sns.color_palette()[1], label = "C")
	patchG = mpatches.Patch(color=sns.color_palette()[2], label = "G")
	patchT = mpatches.Patch(color=sns.color_palette()[3], label = "T")
	patchUAA = mpatches.Patch(color=sns.color_palette()[4], label = "UAA")
	patchUAG = mpatches.Patch(color=sns.color_palette()[5], label = "UAG")
	patchUGA = mpatches.Patch(color=sns.color_palette()[6], label = "UGA")


	plt.legend(handles = [patchA, patchC, patchG, patchT, patchUAA, patchUAG, patchUGA], loc=1, prop={'size':6})

	# plt.axhline(0, linestyle="--", color = "gray")
	ax.set_xticks(range(-upstreamNTs,downstreamNTs+1))
	ax.axhline(y=0, color="gray", linestyle="--")
	# ax.set_ylim(-0.07, 0.07)

	plotoutfile = "/home/jamie/misc/rt/ridge/%s_RRTS_RidgeReg_StopCollapse_%s_to_%s.pdf" % (samp, upstreamNTs, downstreamNTs)

	plt.savefig(plotoutfile, format='pdf')



def plot_coef_abs(dfout, rtvals, alphA, samp):

	X = dfout.values 
	y = rtvals 

	clf = Ridge(alpha=alphA)
	clf.fit(X, y)

	outcoeffs = clf.coef_
	outcoeffs

	outcoef = pd.DataFrame(outcoeffs, index = dfout.columns, columns=['coef'])
	outcoef['abs'] = outcoef['coef'].abs()


	import matplotlib.lines as mlines

	fig, ax = plt.subplots(figsize=(6,6))

	posCol = ['-6','-5','-4','-3','-2','-1',"stop",'+4', '+5', '+6', '+7', '+8', '+9']
	print outcoef

	coef_list = [
		outcoef.iloc[0:4,1].sum(),
		outcoef.iloc[4:8,1].sum(),
		outcoef.iloc[8:12,1].sum(),
		outcoef.iloc[12:16,1].sum(),
		outcoef.iloc[16:20,1].sum(),
		outcoef.iloc[20:24,1].sum(), #up
		outcoef.iloc[24:27,1].sum(), #stop
		outcoef.iloc[27:31,1].sum(),
		outcoef.iloc[31:35,1].sum(),
		outcoef.iloc[35:39,1].sum(),
		outcoef.iloc[39:43,1].sum(),
		outcoef.iloc[43:47,1].sum(),
		outcoef.iloc[47:51,1].sum(),

	]
	print coef_list


	fig, ax = plt.subplots(figsize=(12,6))

	plt.bar(x=posCol, height=coef_list)

	plt.savefig('/home/jamie/misc/rt/abs_test_reg.pdf', format='pdf')

# def plot_logomaker():

def plot_Ridge_Reg_StopCollapse_logmaker(dfout, rtvals, alphA, samp, treat):

	X = dfout.values 
	y = rtvals 

	clf = Ridge(alpha=alphA)
	clf.fit(X, y)

	outcoeffs = clf.coef_
	outcoeffs

	outcoef = pd.DataFrame(outcoeffs, index = dfout.columns)

	print outcoef



	dlm = pd.DataFrame(np.zeros((4,upstreamNTs+3+downstreamNTs)))
	dlm.index = ['A', 'C', 'G', 'U']
	dlm.columns = range(-upstreamNTs, downstreamNTs+3)

	print dlm 

	for i in outcoef.index:
		pos, nt = i.split('_')
		pos = int(pos)
		if nt == 'T':
			nt = 'U'
		coefVal = outcoef.loc[i].item()

	### handle the stop codon: split accross 

		if pos == 0: # this should be the stop codon
			if nt == 'TAA':
				nt = 'A'
			if nt == 'TAG':
				nt = 'C'
			if nt == 'TGA':
				nt = 'G'

			dlm.loc[nt, 0] = coefVal
			dlm.loc[nt, 1] = coefVal
			dlm.loc[nt, 2] = coefVal


		else:

			print pos, nt
			print coefVal
			dlm.loc[nt, pos] = coefVal

	print dlm # now have complete position weight matrix

	dlm = dlm.T
	dlm

	fig, ax = plt.subplots(figsize=(12,6))

	logo = lm.Logo(dlm, 
				   font_name = 'Arial',
				   ax=ax
				  )

	logo.draw_baseline(linewidth = 1, color ='black', linestyle="-")
	logo.highlight_position_range(0,2,alpha=0.5,color='lightgray')

	# ax.set_ylim(-7, 7)

	# outpath = '%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_sccLogo.pdf' % (
	# 	rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	outpath = '%s/figures/Fig5S1B_%s.pdf' % (rootDir, treat)

	plt.savefig(outpath, format='pdf')

	###


	# import matplotlib.lines as mlines

	# fig, ax = plt.subplots(figsize=(6,6))

	# # ntCols = {"A": "blue",
	# # 		  "C" : "orange",
	# # 		  "G" : "green",
	# # 		  "T" : "red",
	# # 		  "TAA": "purple",
	# # 		  "TAG": "brown",
	# # 		  "TGA": "pink"}

	# ntCols = {"A": sns.color_palette()[0],
	# 	  "C" : sns.color_palette()[1],
	# 	  "G" : sns.color_palette()[2],
	# 	  "T" : sns.color_palette()[3],
	# 	  "TAA": sns.color_palette()[4],
	# 	  "TAG": sns.color_palette()[5],
	# 	  "TGA": sns.color_palette()[6]}

	# markerDict = {"A": "v",
	# 	  "C" : "s",
	# 	  "G" : "P",
	# 	  "T" : "D",
	# 	  "TAA": "o",
	# 	  "TAG": "p",
	# 	  "TGA": "h"}

	# for pos in outcoef.index:
	# #     print pos, coeffs_lm.loc[pos]
		
	# 	ntpos = int(pos.split("_")[0])
	# 	nt = pos.split("_")[1]
		
	# 	if ntpos > 0:
	# 		ntpos = ntpos-2
			
	# 	val = outcoef.loc[pos]
		
	# 	plt.plot(ntpos, val, marker=markerDict[nt], color = ntCols[nt], markersize = 6)
		
		
	# # patchA = mpatches.Patch(color="blue", label = "A")
	# # patchC = mpatches.Patch(color="orange", label = "C")
	# # patchG = mpatches.Patch(color="green", label = "G")
	# # patchT = mpatches.Patch(color="red", label = "T")
	# # patchUAA = mpatches.Patch(color="purple", label = "UAA")
	# # patchUAG = mpatches.Patch(color="brown", label = "UAG")
	# # patchUGA = mpatches.Patch(color="pink", label = "UGA")

	# patchA = mpatches.Patch(color=sns.color_palette()[0], label = "A")
	# patchC = mpatches.Patch(color=sns.color_palette()[1], label = "C")
	# patchG = mpatches.Patch(color=sns.color_palette()[2], label = "G")
	# patchT = mpatches.Patch(color=sns.color_palette()[3], label = "T")
	# patchUAA = mpatches.Patch(color=sns.color_palette()[4], label = "UAA")
	# patchUAG = mpatches.Patch(color=sns.color_palette()[5], label = "UAG")
	# patchUGA = mpatches.Patch(color=sns.color_palette()[6], label = "UGA")


	# plt.legend(handles = [patchA, patchC, patchG, patchT, patchUAA, patchUAG, patchUGA], loc=1, prop={'size':6})

	# # plt.axhline(0, linestyle="--", color = "gray")
	# ax.set_xticks(range(-upstreamNTs,downstreamNTs+1))
	# ax.axhline(y=0, color="gray", linestyle="--")
	# # ax.set_ylim(-0.07, 0.07)

	# plotoutfile = "/home/jamie/misc/rt/ridge/%s_RRTS_RidgeReg_StopCollapse_%s_to_%s.pdf" % (samp, upstreamNTs, downstreamNTs)

	# plt.savefig(plotoutfile, format='pdf')


def main():
	df, dflist = build_master_df()
	print df.head()
	
	for i in range(len(samplelist)):

		# if i != 6:
		# 	continue

		samp = samplelist[i]
		treat = treatlist[i]
		d = dflist[i]
		d = d.merge(mrnadf, on="#transcript")
		d = d.merge(scdf, on="#transcript")

		# print d.head()

		d = find_stop_context(d, upstreamNTs, downstreamNTs)
		d = scale_RRTS(d, 1000)
		d = winsorize_RRTS_scaled(d, 0.01)

		print d.head()

		# dfout = build_nt_matrix(d, upstreamNTs, downstreamNTs)
		dfoutStopCollapse = build_nt_matrix_collapseStop(d, upstreamNTs, downstreamNTs)


		print dfoutStopCollapse.head(20)
		# rtvals = get_RRTS_vals(d, dfout)
		rtvals = get_RRTS_vals(d, dfoutStopCollapse)

		print rtvals.head()

		plot_Ridge_Reg_StopCollapse_logmaker(dfoutStopCollapse, rtvals, 1000.0, samp, treat)

if __name__ == "__main__":
	main()










