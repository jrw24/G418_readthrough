"""
Create a metagene around the first inframe stop codon
"""

### import dependencies
import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines


import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import struct
import pandas as pd
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
import rphelper as rph


rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)



### inputs


readingFrames = ['frameMinusOne', 'frameZero', 'framePlusOne']
figPositoin = ['left', 'center', 'right']

samplelist = samplelist[23:25]
samplenames = '--'.join(samplelist)
threshold = 10

pop = "fl"


assignment = "5"
ribosome_shift = "A"
norm_type = 'rpm'
densitystring = 'Density_rpm'

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


firstStopsCSV = stopcodons

upstreamNTs = 12
downstreamNTs = 12







def load_genomes(UTRfilestring, firstStopsCSV, twobitfile):
	"""
	make this a separate function so that these only need to be loaded a single time
	"""
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	utr3adj = pd.read_csv(firstStopsCSV, index_col=0)
	genome= twobitreader.TwoBitFile(twobitfile) # do we actually need to load this in here?
	return UTRdict, utr3adj, genome



def build_avggene_firstStop(UTRdict, utr3adj, threshold, file, frame_selection):
	"""
	frame_selection == ['frameMinusOne', 'frameZero', 'framePlusOne']
	Not equally weighted here
	"""
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, file)
	trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
		fp_assign_path, densitystring, assignment, ribosome_shift, 
		minlen, maxlen, file, minlen, maxlen, file, minlen, maxlen)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, file)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	
	### This is where all the count files are loaded into a dictionary
	trspdict= rph.readcountsf(trspfilestring) ### this takes a minute

	## create a list of 0's that is the length of the region of interest
		## add 3 to account for the stop codon
		## this will be added to for every transcript
	averagegene= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)] # add 3 for the inframe stop codon

	## add counters and set to zero
	noUTRentry = 0 ### discard transcripts not in UTRdict
	zeroUtrlen = 0 ### discard transcripts with zero 3'UTR length
	zeroCdsdense = 0 ### discard transcripts with zero reads in CDS
	lowCdsdense = 0 ## Optional CDS density thresholding, in RPKM
	noInframeStops = 0 ### discard transcripts without any inframe stop codons
	tooCloseToCDS = 0 ### not enough space between first inframe stop and Coding sequence
	utr3tooShort = 0 ### not enough 3'UTR in region of interest past first stop codon
	totalCountedTranscripts = 0 ## number included in final output

	### Iterated through transcripts one at a time, retrieving counts in region of interest:
	for trsp in trspdict:
		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue
		
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])
		assert mrnalen == cdslen + utr5len + utr3len ## check that this is true
		
		if utr3len == 0: # only check the 3'UTR here 
			zeroUtrlen += 1
			continue

		### set stop selection:
		if frame_selection == 'frameZero':

			utr3LenAdj = int(utr3adj.loc[trsp].frameZeroUtr3LenAdj)+3 # add plus 3 to include length of stop codon
			### utr3LenAdj is now the 0-based position one nucleotide past the first inframe stop codon relative to the 3'UTR
			### position of stop codon is utr3CountList[utr3LenAdj-3:utr3LenAdj]
			### possibly a stupid way to do this, but ill think about it more
			inframeStopCount = int(utr3adj.loc[trsp].frameZeroStopCount)

		elif frame_selection == 'framePlusOne':
			utr3LenAdj = int(utr3adj.loc[trsp].framePlusOneUtr3LenAdj)+3
			inframeStopCount = int(utr3adj.loc[trsp].framePlusOneStopCount)

		elif frame_selection == 'frameMinusOne':
			utr3LenAdj = int(utr3adj.loc[trsp].frameMinusOneUtr3LenAdj)+3
			inframeStopCount = int(utr3adj.loc[trsp].frameMinusOneStopCount)

		else:
			print "frame not set!"
			sys.exit()

		### filter out trsp's with no inframe stops, or stops too close to CDS or end of transcript
		if inframeStopCount == 0:
			noInframeStops += 1
			continue

		if utr3LenAdj-upstreamNTs-3 < 0: # taking everything up to the stop codon, includes first nt of utr3
			tooCloseToCDS += 1
			continue
	
		if utr3LenAdj+downstreamNTs > utr3len:
			utr3tooShort += 1
			continue

		### Load in count file for the transcript here
		exonsplicedcounts = trspdict[trsp]
		
		### define Coding sequence here
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len # cdsend == first position of utr3
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript "+ trsp
			sys.exit()

		cdscounts= exonsplicedcounts[cdsstart:cdsend]
		utr3CountList = exonsplicedcounts[cdsend:]

		cdsdensity= sum(cdscounts)/len(cdscounts)

		if cdsdensity == 0:
			zeroCdsdense += 1
			continue

		if cdsdensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
			lowCdsdense += 1
			continue

		### create the counts we want for our region
		avgGeneCounts = utr3CountList[(utr3LenAdj-upstreamNTs-3):(utr3LenAdj+downstreamNTs)]
		totalCountedTranscripts += 1

		for i in range(len(avgGeneCounts)):
			averagegene[i] += avgGeneCounts[i]

	print "Avggene run compolete for sample %s" % file
	print "Frame selection is: %s" % frame_selection
	print "Genes inculded %s" % totalCountedTranscripts
	print "Number of transcripts absent in UTRfile: %s" % noUTRentry
	print "Number of transcripts with zero UTR lengths: %s" % zeroUtrlen
	print "Number of transcripts with zero CDS density: %s" % zeroCdsdense
	print "Number of transcripts below CDS density threshold: %s" % lowCdsdense 
	print "Number of transcripts with no In-frame stops: %s" % noInframeStops
	print "Number of transcripts too close to normal stop codon: %s" % tooCloseToCDS
	print "Number of transcripts with insufficient 3'UTR: %s" % utr3tooShort
	print "- - - - - - - -"
	print "RPM: %s" % averagegene
	print "- - - - - - - -"

	return averagegene


def plot_avGenes (df_in, frameChoice):
	plotout = '%s/figures/Fig3E_%s.pdf' % (rootDir, frameChoice)

	mycolors = matplotlib.cm.get_cmap('tab20').colors + matplotlib.cm.get_cmap('tab20b').colors
	mycolors = ['black', 'orange']
	# mycolors = ['orange']
	dfplt = df_in

	xAxisRange = range(upstreamNTs+3+downstreamNTs)
	ticRange = xAxisRange[1::3]


	# ticRangeText = ticRange
	ticMid = (upstreamNTs+3+downstreamNTs)/2

	fig, ax = plt.subplots(figsize=(7,7))

	lineRange = ticRange
	for xval in lineRange:
		ax.axvline(xval, color="#D3D3D3", linestyle="--")

	for i in range(len(dfplt.columns)):
		dfplt.iloc[:,i].plot.line(use_index = True, ax = ax, legend = True, style = "-", color=mycolors[i])

	
	ax.legend(loc=2, prop={'size': 6})
	ax.set_xticks(ticRange)
	ax.set_xticklabels([-12,-9,-6,-3,0,3,6,9,12])
	ax.set_ylim(0, 600)

	plt.savefig(plotout, format="pdf", bbox_inches = "tight")


def main():


	for fr in readingFrames:
		frameChoice = fr
		csvOutfile = '%s/figures/Fig3E_%s.csv' % (rootDir, frameChoice)

		if not os.path.exists(csvOutfile):
			UTRdict, utr3adj, genomes = load_genomes(UTRfilestring, firstStopsCSV, twobitfile)
			avGeneList = []
			for file in samplelist:
				avGene = build_avggene_firstStop(UTRdict, utr3adj, threshold, file, frameChoice)
				avGeneList.append(avGene)

			df_avgene = pd.DataFrame(avGeneList).T
			df_avgene.columns = samplelist
		else:
			df_avgene = pd.read_csv(csvOutfile, index_col=0)


		# df_avgene.to_csv(csvOutfile)
		plot_avGenes(df_avgene, frameChoice)


if __name__ == '__main__':
	main()