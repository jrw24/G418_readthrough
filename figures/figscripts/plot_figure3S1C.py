


"""
Script for plotting read size distributions for each region of an mRNA
	Can then look at 5'UTR vs CDS vs 3'UTR vs mRNA
	Relies on inputs from densebuilder: Unnormalized, 5'Mapped

"""


import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
import pandas as pd
import pysam
import numpy as np
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


samplelist = [samplelist[18]]+[samplelist[20]]

### set inputs here:

defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
# insets= defaultInsets
insets= zeroInsets


UTRdict = rph.readindict(open(UTRfilestring, "rU"))


# samplelist = samplelist[5:9]

# for readsize in 

def region_size_dist_ftsize(readsize, sample):
	"""
	build the counts for a given region 
	"""
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, sample)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, sample)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	print "total reads for sample %s = %s" % (sample, totreads)

	# for readsize in ftsize:
	readsize = str(readsize) # convert to string
	trspdictfilestring = '%s/DensityUnnormalized/density5p_0shift_%s/%s_%sf/%s_%sf_' %(
		fp_assign_path, readsize, sample, readsize, sample, readsize)

	bamfilepath_readsize = '%s/%s_star_default/%s_%s_match.sorted.bam' % (
		fp_assign_path, sample, sample, readsize)
	bamfile = pysam.AlignmentFile(bamfilepath_readsize, 'rb')
	read_count_bam = bamfile.count()
	print "total reads in bamfile for sample: %s, read length: %s, equals == %s" % (sample, readsize, read_count_bam)

	## check total number of reads in this bamfile:

	## build the trspdict now for a given readlength:
	trspdict = rph.readcountsf(trspdictfilestring)

	## add counters 
	noUTRentry = 0
	zeroUtrlen = 0
	zeroUtrlenInsets = 0
	zeroCdsdense = 0
	lowCdsdense = 0
	lowCdsCounts = 0 # adding cds raw read counter

	### counters for output
	totUtr5Counts = 0
	totCdsCounts = 0
	totUtr3Counts = 0
	totMrnaCounts = 0

	## iterate through every transcript in the gtf file
	for trsp in trspdict:
		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue
	# from csv: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt
	# position in list: key, 0,		1,		2,		3,		4,		5,		6,		7,		8,			9

	# define base region sizes from UTRdict
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])

		### Not sure if I want to keep this here... see how many have lengths of zero first
		if utr5len == 0:
			zeroUtrlen +=1
			# print("transcript has zero utr5 len %s") % trsp
			# sys.exit() 
			continue
		if utr3len == 0:
			zeroUtrlen +=1
			continue

# get counts from density file
		exonsplicedcounts = trspdict[trsp]

		# set starts and ends 
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript %s" % trsp
			sys.exit()


		# modify region lengths using insets:
		
		utr5len = utr5len-insets['utr5Inset3']
		cdslen = cdslen-insets['cdsInset5']-insets['cdsInset3']
		utr3len = utr3len-insets['utr3Inset5']
		mrnalen = utr5len+cdslen+utr3len

		if utr5len == 0:
			zeroUtrlenInsets +=1
			# print "transcript has zero utr5 len %s" % trsp
			# sys.exit() 
			continue
		if utr3len == 0:
			zeroUtrlenInsets +=1
			continue

		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		mrnaCounts = utr5Counts+cdsCounts+utr3Counts

		totUtr5Counts += utr5Counts
		totCdsCounts += cdsCounts
		totUtr3Counts += utr3Counts
		totMrnaCounts += mrnaCounts
	print "UTR5total = %s, CDStotal = %s, UTR3total = %s, mRNAtotal = %s " % (totUtr5Counts, totCdsCounts, totUtr3Counts, totMrnaCounts)
	return [int(readsize), totUtr5Counts, totCdsCounts, totUtr3Counts, totMrnaCounts]



def main():

	counter = 0
	treatlist = ['Untr', 'G418']
	for sample in samplelist:

		treat = treatlist[counter]
		# ftsize = [31]
		# region_size_dist_ftsize(sample, 15)
		sl = [sample] * len(ftsize)
		p = Pool(nodes=40) 
		results = p.map(region_size_dist_ftsize, ftsize, sl)
		# print results

		dfcols = ['readLength', 'utr5', 'cds', 'utr3', 'mrna']
		df = pd.DataFrame.from_records(results, columns = dfcols)

		df = df.set_index('readLength')
		# print df

		utr5reads = df['utr5'].sum()
		cdsreads = df['cds'].sum()
		utr3reads = df['utr3'].sum()
		mrnareads = df['mrna'].sum()

		df['utr5'] = df['utr5']/utr5reads
		df['cds'] = df['cds']/cdsreads
		df['utr3'] = df['utr3']/utr3reads
		df['mrna'] = df['mrna']/mrnareads

		df.drop('utr5', axis=1, inplace=True)
		df.drop('mrna', axis=1, inplace=True)

		# print df
		outfig = '%s/figures/Fig3S1C_%s.pdf' % (rootDir, treat)

		fig, ax = plt.subplots()

		df.plot.line()
		plt.xticks([15,20,25,30,35,40])
		plt.title(sample)

		plt.savefig(outfig, format='pdf', bbox_inches="tight")
		plt.close()
		counter +=1
		# utr5df = pd.DataFrame()

		# for rl in results:
			


# 	###




if __name__ == '__main__':
	main()
