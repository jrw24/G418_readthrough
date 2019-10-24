#####

import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import math
from scipy import stats
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import twobitreader
import ast
import argparse
import importlib
pd.set_option('display.max_columns', 80)

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
import rphelper as rph


### inputs

samps = [20]

# samp = samplelist[sampInd]

mRNAdf = pd.read_csv('%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_mRNAseqs.csv' % rootDir,sep=',', index_col=0)
firstStopsCSV = stopcodons
stopPositions = pd.read_csv('%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_utr3StopPositions.csv' % rootDir ,sep=',', index_col=0)

### nucleotides to include in the window
upstreamNTs = 30
downstreamNTs = 30
threshold = 10

pop = "fl"
# pop = "custom"
# customSize = 32

assignment = "5"
ribosome_shift = "A"
norm_type = 'rpm'
densitystring = 'Density_rpm'
inset_choice = 'default'
norm = 'eq'

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "read lengths not set"


def load_genomes(UTRfilestring, firstStopsCSV, twobitfile):
	"""
	make this a separate function so that these only need to be loaded a single time
	"""
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	utr3adj = pd.read_csv(firstStopsCSV, index_col=0)
	genome= twobitreader.TwoBitFile(twobitfile) # do we actually need to load this in here?
	return UTRdict, utr3adj, genome

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
	FPassignpath = "%s/FPassignment/%s/%s" % (rootpath, genome_name, experiment)
	namelist = []
	dflist= []
	for samp in samplelist:
		###
		dftemp = pd.read_csv('%s/%s/countTables/%s_%s_%s_%sto%s_countTable_rpkm_utr3adj.csv' % (
			FPassignpath, samp, samp, pop, norm_type, minlen, maxlen))
		dftemp['sampname'] = samp
		dftemp['cdsCountsLog2'] = dftemp['cdsCounts'].apply(log_trans_b2)
		dftemp['RAWcdsCountsLog2'] = dftemp['RAWcdsCounts'].apply(log_trans_b2)
		dftemp['utr3CountsLog2'] = dftemp['utr3Counts'].apply(log_trans_b2)
		dftemp['utr3AdjCountsLog2'] = dftemp['utr3AdjCounts'].apply(log_trans_b2)
		dftemp['utr3OccLog2'] = dftemp['utr3_occupancy'].apply(log_trans_b2)
		dftemp['RRTSlog2'] = dftemp['RRTS'].apply(log_trans_b2)
	#     print dftemp.head()

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > 128]
#         print len(dftemp)

		dflist.append(dftemp)


		namelist.append(samp)

	# ### combine into one master dataframe
	df = pd.concat(dflist, axis=0, ignore_index = True)
	
	return df, dflist, namelist

def filter_transcript_features(stops):
	### require that all transcripts have:
	
	dropList = [] ## add index of transcripts to be removed
	
	zero3TC = 0
	ntcTooClose = 0
	tc3TooClose = 0
	utr3TooShortDownstream = 0
	secondStopTooClose = 0
	
	
	
	## 1) At least one in-frame stop
	## 2) 3'TC is greater than upstreamNTs downstream of from NTC
	## 3) 3'TC is greater than downstreamNTs downstream of the NTC
	## 4) 3'UTR has at least downstreamNTs from end of 3'UTR
	## 5) no other in-frame 3'TCs are present in downstreamNTs region
	
	for tr in stops.index:
		
		frameZeroStops = stops.loc[tr, 'frameZeroStopCount']
		utr3AdjustedLen = stops.loc[tr, 'frameZeroUtr3LenAdj'] # zero based position of 'U' of stop codon, also 3'UTR adjsuted len
		inFrameStopList = ast.literal_eval(stops.loc[tr, 'frameZeroStopPosUtr3'])
		utr3RawLen = stops.loc[tr, '3utr_len']
		
		if frameZeroStops < 1: 
			dropList.append(tr)
			zero3TC +=1
			
		if utr3AdjustedLen < upstreamNTs:
			dropList.append(tr)
#             print tr
#             print stops.loc[tr, 'frameZeroUtr3LenAdj']
			ntcTooClose +=1
			
		if utr3AdjustedLen < downstreamNTs:
			dropList.append(tr)
			tc3TooClose +=1
			
		if utr3AdjustedLen+3+downstreamNTs > utr3RawLen:
			dropList.append(tr)
			utr3TooShortDownstream +=1
#             print frameZeroStops, utr3RawLen, inFrameStopList
			
			
		if frameZeroStops > 1: # if at least 2 inframe stops in 3'UTR
			minInterStopDist = inFrameStopList[0]+3+upstreamNTs
			if inFrameStopList[1] < minInterStopDist:
				dropList.append(tr)
				secondStopTooClose +=1
#                 print inFrameStopList, tr
				
			
			
			
		## len utr3lenadj is exact length of the 3'UTR without including the stopcodon or 3'TC
		
			
	### drop filtered transcripts from the data frame         
	filtStops = stops.drop(dropList, axis=0)
#     print dropList
	print len(stops) ## total
	print len(set(dropList)) ## unique dropped transcripts 
	print len(filtStops) ## remaining transcripts
	
	return filtStops
	
def retrieve_stop_density(din, samp):
	
	d = din.copy()
	
	dropList = []
	averagegeneNTC= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)]
	averagegeneTC3= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)]
	
	d['NTC_up_counts'] = np.zeros(len(d))
	d['NTC_up_dense'] = np.zeros(len(d))
	d['NTC_down_counts'] = np.zeros(len(d))
	d['NTC_down_dense'] = np.zeros(len(d))
#     d['NTC_ratio'] = np.zeros(len(d))
	d['TC3_up_counts'] = np.zeros(len(d))
	d['TC3_up_dense'] = np.zeros(len(d))
	d['TC3_down_counts'] = np.zeros(len(d))
	d['TC3_down_dense'] = np.zeros(len(d))
#     d['TC3_ratio'] = np.zeros(len(d))
	
	lowCdsdense = 0
	totalCountedTranscripts = 0
	
	
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, samp)
	trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
		fp_assign_path, densitystring, assignment, ribosome_shift, 
		minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, samp)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	
	defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
	zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
	customInsets  = { 'utr5Inset3' : 15, 'cdsInset5' : 24, 'cdsInset3' : 15, 'utr3Inset5' : 15 }

	if inset_choice == "default":
		insets = defaultInsets
	elif inset_choice == "zero":
		insets = zeroInsets
	elif inset_choice == "custom":
		insets = customInsets
	else:
		print "Insets were not set"
		sys.exit()

	### This is where all the count files are loaded into a dictionary
	trspdict= rph.readcountsf(trspfilestring) ### this takes a minute
	
	for tr in d.index:
#         print tr
		
		mrnalen = int(d.loc[tr, 'mrna_len'])
		cdslen = int(d.loc[tr, 'cds_len'])
		utr5len = int(d.loc[tr, '5utr_len'])
		utr3len = int(d.loc[tr, '3utr_len'])
		assert mrnalen == cdslen + utr5len + utr3len
		
		exonsplicedcounts = trspdict[tr]
		utr3LenAdj = int(d.loc[tr, 'frameZeroUtr3LenAdj'])+3 ## adding +3 here to now include length of stop codon
		
		### calculate normalized counts to cds
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len # cdsend == first position of utr3

		cdslenMod = cdslen-insets['cdsInset5']-insets['cdsInset3']
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		cdsDensity = (cdsCounts/cdslenMod) 
		
#         print cdsDensity*float(1000)
		
		if cdsDensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
			lowCdsdense += 1
			dropList.append(tr)
			continue
			
		exonNormCounts = [rpf/cdsDensity for rpf in exonsplicedcounts] ## counts normalized by cds density, using only region within insets
#         print sum(exonNormCounts)
#         print sum(exonNormCounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		utr3CountList = exonNormCounts[cdsend:]
		
		
		ntcCounts = exonNormCounts[(cdsend-upstreamNTs-3):(cdsend+downstreamNTs)]
		tc3Counts = utr3CountList[(utr3LenAdj-upstreamNTs-3):(utr3LenAdj+downstreamNTs)]
		
#         print len(ntcCounts)
#         print ntcCounts
#         print len(tc3Counts)
#         print tc3Counts

		### define regions: upstream counts == [0:27], downstream counts == [36:0]
		### more generally: upstream counts == [:upstreamNTs-3], downstream counts == [downstreamNTs+6:]

		d.at[tr, 'NTC_up_counts'] = sum(ntcCounts[:upstreamNTs-3])
		d.at[tr, 'NTC_up_dense'] = sum(ntcCounts[:upstreamNTs-3])/(upstreamNTs-3) ## avoiding codon before stop
		d.at[tr, 'NTC_down_counts'] = sum(ntcCounts[downstreamNTs+6:])
		d.at[tr, 'NTC_down_dense'] = sum(ntcCounts[downstreamNTs+6:])/(downstreamNTs-3) ## avoiding codon after stop
#         d.at[tr, 'NTC_ratio'] = d.loc[tr, 'NTC_down_dense']/d.loc[tr, 'NTC_up_dense']
		d.at[tr, 'TC3_up_counts'] = sum(tc3Counts[:upstreamNTs-3])
		d.at[tr, 'TC3_up_dense'] = sum(tc3Counts[:upstreamNTs-3])/(upstreamNTs-3)
		d.at[tr, 'TC3_down_counts'] = sum(tc3Counts[downstreamNTs+6:])
		d.at[tr, 'TC3_down_dense'] = sum(tc3Counts[downstreamNTs+6:])/(downstreamNTs-3)
#         d.at[tr, 'TC3_ratio'] = d.loc[tr, 'TC3_down_dense']/d.loc[tr, 'TC3_up_dense']
		
		
		for i in range(len(ntcCounts)):
			averagegeneNTC[i] += ntcCounts[i]
		for i in range(len(tc3Counts)):
			averagegeneTC3[i] += tc3Counts[i]
		totalCountedTranscripts +=1
	
	avgGeneEqualNTC = [rpf/totalCountedTranscripts for rpf in averagegeneNTC] ### test average gene plot
	avgGeneEqualTC3 = [rpf/totalCountedTranscripts for rpf in averagegeneTC3] ### test average gene plot
	
#     print "avgGeneNTC!", avgGeneEqualNTC
#     print "avgGeneTC3!", avgGeneEqualTC3
	
	
			
		
			
			

	dout = d.drop(dropList, axis=0)
	
	print len(d)
	print len(dropList)
	print len(dout)
		
#         print exonsplicedcounts
		

	return avgGeneEqualNTC, avgGeneEqualTC3, dout
	
# retrieve_stop_density(d.head(5))
# retrieve_stop_density(pd.DataFrame(d.loc['ENST00000366695.3']).T)
# v1, v2 = retrieve_stop_density(d)
	
def calc_rt_ratios(ntcAvGene, tc3AvGene):

#     print len(dfout)
#     dftemp = dfout.loc[dfout['NTC_up_dense'] > 0]
#     print len(dftemp)
#     dfp = dftemp.loc[dftemp['TC3_up_dense']> 0]
#     print len(dfp)

#     dfratio = dfp.copy()

#     for tr in dfratio.index:

#         ntcRatio = dfp.loc[tr, 'NTC_down_dense']/dfp.loc[tr, 'NTC_up_dense']
#         tc3Ratio = dfp.loc[tr, 'TC3_down_dense']/dfp.loc[tr, 'TC3_up_dense']
#         dfratio.at[tr, 'NTC_ratio'] = ntcRatio
#         dfratio.at[tr, 'TC3_ratio'] = tc3Ratio

	# ddd['NTC_ratio'] = dfp['NTC_down_dense']/dfp['NTC_up_dense']
	
	v1 = ntcAvGene
	v2 = tc3AvGene
	
	ntcUP = sum(v1[:(upstreamNTs-3)])/(upstreamNTs-3)
	ntcDN = sum(v1[(downstreamNTs+6):])/(downstreamNTs-3)

	ntcRT = ntcDN/ntcUP


	tc3UP = sum(v2[:(upstreamNTs-3)])/(upstreamNTs-3)
	tc3DN = sum(v2[(downstreamNTs+6):])/(downstreamNTs-3)
	tc3RT = tc3DN/tc3UP
	
	return ntcRT, tc3RT


def plot_rt_ratio(ntcRT, tc3RT):
	rtvals = ntcRT + tc3RT
	print rtvals

	dfdict = {'rt': rtvals,
			 'sample':samps+samps,
			 'val':np.repeat(np.array(['ntc', 'tc3']),[len(samps),len(samps)]).tolist(),
			 }

	d = pd.DataFrame.from_dict(dfdict)
	d

	sns.catplot(data=d, x='val', y='rt', kind='swarm',  height=3, aspect=0.7, edgecolor="black", linewidth=1, s = 6)
	sns.boxplot(data=d, x='val', y='rt', showbox=False , width = 0.5, showcaps=False,
			   whiskerprops= dict(color="white"))

	ax = plt.gca()

	ax.set_ylim(0,0.5)


	# figout = '/home/jamie/misc/rt/utr3_rt_predict_vs_actual_update_191002.pdf'
	figout = '%s/figures/Fig5D.pdf' % rootDir
	# plt.savefig(figout, format='pdf', bbox_inches="tight")

def plot_avGenes(ntcAvGene, tc3AvGene):
	
	v1 = ntcAvGene
	v2 = tc3AvGene
	
	ntcUP = sum(v1[:(upstreamNTs-3)])/(upstreamNTs-3)
	ntcDN = sum(v1[(downstreamNTs+6):])/(downstreamNTs-3)

	ntcRT = ntcDN/ntcUP


	tc3UP = sum(v2[:(upstreamNTs-3)])/(upstreamNTs-3)
	tc3DN = sum(v2[(downstreamNTs+6):])/(downstreamNTs-3)
	tc3RT = tc3DN/tc3UP
	
	rtvals = ntcRT + tc3RT
	print rtvals

	### NTC
	fig, ax = plt.subplots(figsize=(9,6))
	pd.DataFrame(v1).plot.line(xticks=range(len(v1))[1::3], ax=ax)
	ax.set_xticklabels(range(-upstreamNTs, downstreamNTs+3)[::3])
	
	plt.plot((0,upstreamNTs-3-1), (ntcUP, ntcUP), color="black", linestyle="--")
	plt.plot((36,33+downstreamNTs-1), (ntcDN, ntcDN), color="black", linestyle="--")
#     plt.plot((30,32), (ntcUP,ntcUP), 'k-')
#     plt.plot((30,32), (ntcDN,ntcDN), 'k-')
#     plt.plot((31,31), (ntcUP, ntcDN), 'k-')

	figout = '%s/figures/Fig5S2A.pdf' % rootDir
	plt.savefig(figout, format='pdf', bbox_inches="tight")
	
	### TC3
	fig, ax = plt.subplots(figsize=(9,6))
	pd.DataFrame(v2).plot.line(xticks=range(len(v1))[1::3], ax=ax)
	ax.set_xticklabels(range(-upstreamNTs, downstreamNTs+3)[::3])
	
	plt.plot((0,upstreamNTs-3-1), (tc3UP, tc3UP), color="black", linestyle="--")
	plt.plot((36,33+downstreamNTs-1), (tc3DN, tc3DN), color="black", linestyle="--")
#     plt.plot((30,32), (tc3UP,tc3UP), 'k-')
#     plt.plot((30,32), (tc3DN,tc3DN), 'k-')
#     plt.plot((31,31), (tc3UP, tc3DN), 'k-')

	figout = '%s/figures/Fig5S2B.pdf' % rootDir
	plt.savefig(figout, format='pdf', bbox_inches = "tight")



def main():
	NTC_outlist = []
	TC3_outlist = []

	for s in samps:
		
		sampInd = s
		samp = samplelist[sampInd]
		print samp

		UTRdict, utr3adj, genomes = load_genomes(UTRfilestring, firstStopsCSV, twobitfile)
		df, dflist, namelist = load_samples() ### samples are filtered on 128 RAW reads
		stops = utr3adj.merge(stopPositions, left_index=True, right_index=True )
		filtStops = filter_transcript_features(stops)

		dfm = dflist[sampInd]

		dfm = dfm.set_index('#transcript', drop=True)

		d = filtStops.merge(dfm, left_index=True, right_index=True)

		print len(d)
		d.head()

		ntcAvGene, tc3AvGene, dfout = retrieve_stop_density(d, samp)
		
		ntcRT, tc3RT = calc_rt_ratios(ntcAvGene, tc3AvGene)
		print ntcRT, tc3RT
		
		NTC_outlist.append(ntcRT)
		TC3_outlist.append(tc3RT)
		
	# plot_rt_ratio(NTC_outlist, TC3_outlist)
	plot_avGenes(ntcAvGene, tc3AvGene)


if __name__ == '__main__':
	main()



