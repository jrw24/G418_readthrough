#####

#####

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


### inputs

samplelist = [samplelist[18]]
treatList = ['Untr']

upstreamNts = 12
downstreamNts = 12

mrnadf = pd.read_csv("%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_mRNAseqs.csv" % rootDir)
scdf = pd.read_csv(stopcodons)

RAW_reads_thresh = 128

pop = "fl"
assignment = "5"
ribosome_shift = "A"
norm = "rpm"

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
else:
	print "oops"
	

### define functions:

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
		# print samp 
	#     dftemp = pd.read_csv('%s/%s/countTables/%s_fl_rpm_28to35_countTable_rpkm.csv' % (FPassignpath, samp, samp))
	#     dftemp['sampname'] = samp
	# #     print dftemp.head()
	#     dflist.append(dftemp)
		
		###
		
		dftemp = pd.read_csv('%s/%s/countTables/%s_%s_%s_%sto%s_countTable_rpkm_utr3adj.csv' % (
			FPassignpath, samp, samp, pop, norm, minlen, maxlen))
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

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > RAW_reads_thresh]
		print len(dftemp)

		dflist.append(dftemp)
		
		
		namelist.append(samp)
		
	# ### combine into one master dataframe
	# df = pd.concat(dflist, axis=0, ignore_index = True)
	df = pd.concat(dflist, axis=0, ignore_index = True)

	return df, dflist


def find_stop_context(d, upstreamNts, downstreamNts):

	print len(d)
	
	for i in range(len(d)):
		sequence = d.loc[i,'mRNAseqs']
		pos = d.loc[i,'5utr_len'] + d.loc[i,'cds_len'] 

		d.loc[i,'stopContext'] = sequence[pos-upstreamNts-3:pos+downstreamNts].upper() # set to 6 and 6 by default

		# print d.loc[i,'stopContext']

		if len(d.loc[i, 'stopContext']) < upstreamNts+3+downstreamNts: ### filter out sequences that are too short
			# print d.loc[i,'stopContext']
			d.drop(i, axis=0, inplace=True)

		# d.loc[i,'stopContext'] = sequence[pos-9:pos+6].upper()    
#         d.loc[i,'nt4_1stInFrame'] = sequence[pos:pos+4].upper()

	print len(d)
	return d


def find_stop_context_3tc(d, upstreamNts, downstreamNts):

	print len(d)

	print d.head()

	scdf = pd.read_csv(stopcodons, index_col=0)
	# print scdf.head()
	
	for i in range(len(d)):
		tr = d.loc[i, '#transcript']
		sccLen = upstreamNts+3+downstreamNts
		sequence = d.loc[i,'mRNAseqs']
		utr3LenAdj = scdf.loc[tr, 'frameZeroUtr3LenAdj']
		pos = d.loc[i,'5utr_len'] + d.loc[i,'cds_len']+utr3LenAdj+3 ## need to add +3 to inlucde stop codon

		# if tr == 'ENST00000222803.9': ### this checks out
		# 	print sequence[pos-upstreamNts-3:pos+downstreamNts].upper()

		if upstreamNts > utr3LenAdj: ## set too short sequences to no context
			d.loc[i,'stopContext3tc'] == ""
		else:
			d.loc[i,'stopContext3tc'] = sequence[pos-upstreamNts-3:pos+downstreamNts].upper() # set to 6 and 6 by default

	for i in d.index:
		if len(d.loc[i, 'stopContext3tc']) < sccLen:
			d.drop(i, axis=0, inplace=True)

	d.reset_index(drop=True, inplace=True)

	print len(d)
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




def build_stopseq_RRTS_kpLogo(dfinput, samp):

	df = dfinput.copy()
	seqdf = df[['stopContext3tc', 'RRTS']]

	print seqdf.head()
	print len(seqdf)


	# t = 'AAGCGATGAGGGCCC'

	# print len(t)


	seqlist = seqdf['stopContext3tc']

	outseqs = []

	for seq in seqlist:
		outlist = []

		counter = 0

		for i in seq:
			# print i, counter

			if counter <upstreamNts:
				outlist.append(i)
				counter+=1
			elif counter == upstreamNts:
				stopCodon = seq[counter:counter+3]
				# print stopCodon
				stopDict = {"TAA":"A",
							"TAG":"C",
							"TGA":"G"}
				outlist.append(stopDict[stopCodon])
				counter+=1
			elif counter < upstreamNts+3:
				counter +=1 # wait for stop codon to catch up
			else:
		#         print i
				outlist.append(i)
				counter +=1

		outcontext = "".join(outlist)
		outseqs.append(outcontext)

	# outdf = pd.DataFrame({'stopOut':outseqs,
	# 					  'RRTS':seqdf['RRTS']})

	print len(outseqs)


	seqdf2 = seqdf.copy()
	seqdf2['stopCollapse'] = outseqs

	print seqdf2.head()

	outdf = seqdf2[['stopCollapse', 'RRTS']]

	print outdf.head()


	outdir = "%s/FPassignment/%s/%s/analysis/kpLogo" % (rootpath, genome_name, experiment)
	if not os.path.exists(outdir): os.makedirs(outdir)

	sampledir = "%s/%s" % (outdir, samp)
	if not os.path.exists(sampledir): os.makedirs(sampledir)

	outdf.to_csv("%s/%s_%supTO%sdown_kpLogoTest_stopCollapse.txt" % (sampledir, samp, upstreamNts, downstreamNts), sep='\t', header=False, index=False)

	startPosit = upstreamNts+1
	kpLogCmnd = "/home/jamie/tools/kpLogo/kpLogo-1.1/bin/kpLogo %s/%s_%supTO%sdown_kpLogoTest_stopCollapse.txt -k 1 -startPos %s -o %s/%s_%supTO%sdown_kpLogo -weighted -gapped -pc 0.05" % (
		sampledir, samp, upstreamNts, downstreamNts, startPosit, sampledir, samp, upstreamNts, downstreamNts)


	os.system(kpLogCmnd)
	print kpLogCmnd


def plot_logomaker(samp):
	inpath = "%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_kpLogo.pass.p.cutoff.txt" % (
		rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	df = pd.read_csv(inpath, sep = '\t')
	print df 
	df = df.sort_values(by=['position', '#kmer'])
	df = df.reset_index(drop=True)

	### get raw p-values from -log10 vals

	pvals = []

	for i in df['p-value']:
		pvals.append(10 ** (-i))
	df['pvals'] = pvals

	### bonferroni correction
	reject, pvals_cor, sidak, bonfMeth = mtest.multipletests(pvals, alpha=0.1, method='bonferroni')
	df['bonf_pval'] = pvals_cor

	### FRD, Benjamini Hochberg adjusted pvalues
	reject, pvals_cor, sidak, bonfMeth = mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
	
	pvals_cor_log10 = []
	for i in pvals_cor:
		p10 = -math.log10(i)
		pvals_cor_log10.append(p10)

	df['bh_pval'] = pvals_cor
	df['bh_pval_neglog10'] = pvals_cor_log10

	### extract desired values
	dp = df.copy()
	dp = dp[['#kmer', 'position', 'shift','statistics', 'pvals', 'p-value', 'bh_pval', 'bh_pval_neglog10']]
	dp['difflog10'] = dp['p-value']*dp['bh_pval_neglog10']
	dp['log10p'] = dp['p-value']*np.sign(dp['statistics'])
	dp['log10pBH'] = dp['bh_pval_neglog10']*np.sign(dp['statistics'])

	### set up matrix for lm
	dl = dp.copy()
	dl = dl[['#kmer', 'position', 'log10pBH']]
	dl.loc[dl['position']>0,'position'] -=1

	dl.loc[dl['#kmer'] == 'T','#kmer'] = 'U'


	# dlm = pd.DataFrame(np.zeros((4,upstreamNts+1+downstreamNts)))
	# dlm.index = ['A', 'C', 'G', 'T']
	# dlm.columns = range(-upstreamNts, downstreamNts+1)
	# # dlm.drop(0, axis=1,inplace=True)

	# dlm


	# for i in dl.index:
	# 	nt = dl.loc[i,'#kmer']
	# 	pos = dl.loc[i, 'position']
	# 	val = dl.loc[i, 'log10pBH']
		
	# 	dlm.loc[nt, pos] = val

	# dlm = dlm.T
	# dlm

	# fig, ax = plt.subplots(figsize=(6,6))

	# logo = lm.Logo(dlm, 
	# 			   font_name = 'Arial',
	# 			   ax=ax
	# 			  )

	# logo.draw_baseline(linewidth = 1, color ='black', linestyle="-")

	# outpath = '%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_sccLogo.pdf' % (
	# 	rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	# plt.savefig("/home/jamie/misc/rt/logoTest.pdf", format='pdf')

	### rep sc 3 positions

	print dl 
	print ""

	dlm = pd.DataFrame(np.zeros((4,upstreamNts+3+downstreamNts)))
	dlm.index = ['A', 'C', 'G', 'U']
	dlm.columns = range(-upstreamNts, downstreamNts+3)
	# dlm.drop(0, axis=1,inplace=True)

	print dlm


	for i in dl.index:
		nt = dl.loc[i,'#kmer']
		pos = dl.loc[i, 'position']
		val = dl.loc[i, 'log10pBH']
		
		if pos < 0:
			dlm.loc[nt, pos] = val
		elif pos == 0:
			if nt == 'A':
				dlm.loc['A',0] = val
				dlm.loc['A',1] = val
				dlm.loc['A',2] = val
			elif nt == 'C':
				dlm.loc['C',0] = val
				dlm.loc['C',1] = val
				dlm.loc['C',2] = val
			elif nt == 'G':
				dlm.loc['G',0] = val
				dlm.loc['G',1] = val
				dlm.loc['G',2] = val
			else:
				print "stop not set"
		elif pos >0:
			newpos = pos +2
			dlm.loc[nt, newpos] = val

	

	print dlm


	dlm = dlm.T
	dlm

	fig, ax = plt.subplots(figsize=(12,6))

	logo = lm.Logo(dlm, 
				   font_name = 'Arial',
				   ax=ax
				  )

	logo.draw_baseline(linewidth = 1, color ='black', linestyle="-")
	logo.highlight_position_range(0,2,alpha=0.5,color='lightgray')

	ax.set_ylim(-7, 7)

	outpath = '%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_sccLogo_3tc.pdf' % (
		rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	plt.savefig(outpath, format='pdf')


def plot_logomaker_abs(samp):
	inpath = "%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_kpLogo.pass.p.cutoff.txt" % (
		rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	df = pd.read_csv(inpath, sep = '\t')
	print df 
	df = df.sort_values(by=['position', '#kmer'])
	df = df.reset_index(drop=True)

	### get raw p-values from -log10 vals

	pvals = []

	for i in df['p-value']:
		pvals.append(10 ** (-i))
	df['pvals'] = pvals

	### bonferroni correction
	reject, pvals_cor, sidak, bonfMeth = mtest.multipletests(pvals, alpha=0.1, method='bonferroni')
	df['bonf_pval'] = pvals_cor

	### FRD, Benjamini Hochberg adjusted pvalues
	reject, pvals_cor, sidak, bonfMeth = mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
	
	pvals_cor_log10 = []
	for i in pvals_cor:
		p10 = -math.log10(i)
		pvals_cor_log10.append(p10)

	df['bh_pval'] = pvals_cor
	df['bh_pval_neglog10'] = pvals_cor_log10

	### extract desired values
	dp = df.copy()
	dp = dp[['#kmer', 'position', 'shift','statistics', 'pvals', 'p-value', 'bh_pval', 'bh_pval_neglog10']]
	dp['difflog10'] = dp['p-value']*dp['bh_pval_neglog10']
	dp['log10p'] = dp['p-value']*np.sign(dp['statistics'])
	dp['log10pBH'] = dp['bh_pval_neglog10']*np.sign(dp['statistics'])

	print dp['bh_pval_neglog10']

	### set up matrix for lm
	dl = dp.copy()
	dl = dl[['#kmer', 'position', 'bh_pval_neglog10', 'log10pBH']]
	dl.loc[dl['position']>0,'position'] -=1

	dl.loc[dl['#kmer'] == 'T','#kmer'] = 'U'


	# dlm = pd.DataFrame(np.zeros((4,upstreamNts+1+downstreamNts)))
	# dlm.index = ['A', 'C', 'G', 'T']
	# dlm.columns = range(-upstreamNts, downstreamNts+1)
	# # dlm.drop(0, axis=1,inplace=True)

	# dlm


	# for i in dl.index:
	# 	nt = dl.loc[i,'#kmer']
	# 	pos = dl.loc[i, 'position']
	# 	val = dl.loc[i, 'log10pBH']
		
	# 	dlm.loc[nt, pos] = val

	# dlm = dlm.T
	# dlm

	# fig, ax = plt.subplots(figsize=(6,6))

	# logo = lm.Logo(dlm, 
	# 			   font_name = 'Arial',
	# 			   ax=ax
	# 			  )

	# logo.draw_baseline(linewidth = 1, color ='black', linestyle="-")

	# outpath = '%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_sccLogo.pdf' % (
	# 	rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	# plt.savefig("/home/jamie/misc/rt/logoTest.pdf", format='pdf')

	### rep sc 3 positions

	print dl 
	print ""

	dlm = pd.DataFrame(np.zeros((4,upstreamNts+3+downstreamNts)))
	dlm.index = ['A', 'C', 'G', 'U']
	dlm.columns = range(-upstreamNts, downstreamNts+3)
	# dlm.drop(0, axis=1,inplace=True)

	print dlm


	for i in dl.index:
		nt = dl.loc[i,'#kmer']
		pos = dl.loc[i, 'position']
		val = dl.loc[i, 'bh_pval_neglog10']
		
		if pos < 0:
			dlm.loc[nt, pos] = val
		elif pos == 0:
			if nt == 'A':
				dlm.loc['A',0] = val
				dlm.loc['A',1] = val
				dlm.loc['A',2] = val
			elif nt == 'C':
				dlm.loc['C',0] = val
				dlm.loc['C',1] = val
				dlm.loc['C',2] = val
			elif nt == 'G':
				dlm.loc['G',0] = val
				dlm.loc['G',1] = val
				dlm.loc['G',2] = val
			else:
				print "stop not set"
		elif pos >0:
			newpos = pos +2
			dlm.loc[nt, newpos] = val

	

	print dlm


	dlm = dlm.T
	dlm

	fig, ax = plt.subplots(figsize=(12,6))

	logo = lm.Logo(dlm, 
				   font_name = 'Arial',
				   ax=ax
				  )

	# for i in dl.index: # this is not really working...
	# 	if np.sign(dl.loc[i,'log10pBH']) == -1.0:
	# 		# print 'neg', dl.loc[i,['log10pBH']]
	# 		nt = dl.loc[i,'#kmer']
	# 		pos = dl.loc[i, 'position']

	# 		logo.style_single_glyph(pos, nt, flip=True, edgecolor='black', color='white')



	logo.draw_baseline(linewidth = 1, color ='black', linestyle="-")
	logo.highlight_position_range(0,2,alpha=0.5,color='lightgray')

	ax.set_ylim(0, 15)

	outpath = '%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_sccLogo_ABS_3tc.pdf' % (
		rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	plt.savefig(outpath, format='pdf')


def plot_frequency(samp):
	inpath = "%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_kpLogo.pass.p.cutoff.txt" % (
		rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	df = pd.read_csv(inpath, sep = '\t')
	print df 
	df = df.sort_values(by=['position', '#kmer'])
	df = df.reset_index(drop=True)

	### get raw p-values from -log10 vals

	pvals = []

	freqs = []

	for i in df.index:

		count1 = float(df.loc[i,'n1'])
		countOther = float(df.loc[i, 'n2'])
		freqOne = (count1/(count1+countOther))
		freqs.append(freqOne)
	
	df['freq'] = freqs
	### extract desired values

	### set up matrix for lm
	dl = df.copy()
	dl = dl[['#kmer', 'position', 'freq']]
	dl.loc[dl['position']>0,'position'] -=1

	dl.loc[dl['#kmer'] == 'T','#kmer'] = 'U'

	print dl 
	print ""

	dlm = pd.DataFrame(np.zeros((4,upstreamNts+3+downstreamNts)))
	dlm.index = ['A', 'C', 'G', 'U']
	dlm.columns = range(-upstreamNts, downstreamNts+3)
	# dlm.drop(0, axis=1,inplace=True)

	print dlm


	for i in dl.index:
		nt = dl.loc[i,'#kmer']
		pos = dl.loc[i, 'position']
		val = dl.loc[i, 'freq']
		
		if pos < 0:
			dlm.loc[nt, pos] = val
		elif pos == 0:
			if nt == 'A':
				dlm.loc['A',0] = val
				dlm.loc['A',1] = val
				dlm.loc['A',2] = val
			elif nt == 'C':
				dlm.loc['C',0] = val
				dlm.loc['C',1] = val
				dlm.loc['C',2] = val
			elif nt == 'G':
				dlm.loc['G',0] = val
				dlm.loc['G',1] = val
				dlm.loc['G',2] = val
			else:
				print "stop not set"
		elif pos >0:
			newpos = pos +2
			dlm.loc[nt, newpos] = val

	
	print dlm


	dlm = dlm.T
	dlm

	fig, ax = plt.subplots(figsize=(30,6))

	logo = lm.Logo(dlm, 
				   font_name = 'Arial',
				   ax=ax,
				   vpad=0.05
				  )

	logo.draw_baseline(linewidth = 1, color ='black', linestyle="-")
	# logo.highlight_position_range(0,2,alpha=0.5,color='lightgray')

	ax.set_ylim(0,1)

	outpath = '%s/FPassignment/%s/%s/analysis/kpLogo/%s/%s_%supTO%sdown_FREQ_3tc.pdf' % (
		rootpath, genome_name, experiment, samp, samp, upstreamNts, downstreamNts)

	outpath = '%s/figures/Fig5C.pdf' % rootDir
	plt.savefig(outpath, format='pdf', bbox_inches = 'tight')


def main():
	df, dflist = build_master_df()
	print df.head()
	

	for i in range(len(samplelist)):

		samp = samplelist[i]
		treat = treatList[i]

		d = dflist[i]
		d = d.merge(mrnadf, on="#transcript")
		d = d.merge(scdf, on="#transcript")

		d = find_stop_context_3tc(d, upstreamNts, downstreamNts)
		print d[['stopContext3tc', 'RRTS']].head(50)
		d = scale_RRTS(d, 1000)
		d = winsorize_RRTS_scaled(d, 0.01)

		build_stopseq_RRTS_kpLogo(d, samp)
		plot_frequency(samp)

if __name__ == '__main__':
	main()



