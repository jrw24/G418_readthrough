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
import rphelper as rph
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


rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)


### inputs
print samplelist

threshold = 10
customSize = 30
pop = "fl"
region_choice = "start"

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


upstreamNTs = 100 ### keep this value as a positive int
downstreamNTs = 150


mRNAdf = pd.read_csv('%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_mRNAseqs.csv' % rootDir,
						sep=',', index_col=0)


def load_genomes(UTRfilestring, twobitfile):
	"""
	make this a separate function so that these only need to be loaded a single time
	"""
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	genome= twobitreader.TwoBitFile(twobitfile) # do we actually need to load this in here?
	return UTRdict, genome


def build_avggene_firstStop(UTRdict, threshold, samp, region_choice, inset_choice):
	"""
	frame_selection == ['frameMinusOne', 'frameZero', 'framePlusOne']
	"""
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, samp)
	trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
		fp_assign_path, densitystring, assignment, ribosome_shift, 
		minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, samp)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	
	### This is where all the count files are loaded into a dictionary
	trspdict= rph.readcountsf(trspfilestring) ### this takes a minute

	## create a list of 0's that is the length of the region of interest
		## add 3 to account for the stop codon
		## this will be added to for every transcript
	averagegene= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)] # add 3 for start or stop codon

	# print averagegene

	## add counters and set to zero
	noUTRentry = 0 ### discard transcripts not in UTRdict
	zeroCdsdense = 0 ### discard transcripts with zero reads in CDS
	lowCdsdense = 0 ## Optional CDS density thresholding, in RPKM
	regionTooShort = 0 ### not enough 3'UTR in region of interest past first stop codon
	totalCountedTranscripts = 0 ## number included in final output


	### calculate mRNA-region densities 

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

	### Iterated through transcripts one at a time, retrieving counts in region of interest:
	for trsp in trspdict:
		
		# if trsp != 'ENST00000319248.13': ## testing with PRDX1
		# 	continue

		### Load in count file for the transcript here
		exonsplicedcounts = trspdict[trsp]

		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue
		
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])
		assert mrnalen == cdslen + utr5len + utr3len ## check that this is true
		

		### define Coding sequence here
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len # cdsend == first position of utr3
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript "+ trsp
			sys.exit()

		# cdscounts= exonsplicedcounts[cdsstart:cdsend] ### counts are already in rpm
		# utr3CountList = exonsplicedcounts[cdsend:]

		# cdsdensity= sum(cdscounts)/len(cdscounts) ### unmodified CDS density

		# if cdsdensity == 0:
		# 	zeroCdsdense += 1
		# 	continue

		# if cdsdensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
		# 	lowCdsdense += 1
		# 	continue


		### Calculate Region Densities ###
		utr5lenMod = utr5len-insets['utr5Inset3']
		cdslenMod = cdslen-insets['cdsInset5']-insets['cdsInset3']
		utr3lenMod = utr3len-insets['utr3Inset5']
		mrnalenMod = utr5lenMod+cdslenMod+utr3lenMod

		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		mrnaCounts = utr5Counts+cdsCounts+utr3Counts


		### RAW counts
		RAWutr5Counts = int(utr5Counts*(totreads/1E6))
		RAWutr3Counts = int(utr3Counts*(totreads/1E6))
		RAWcdsCounts = int(cdsCounts*(totreads/1E6))
		RAWmrnaCounts = int(mrnaCounts*(totreads/1E6))

		### denisites
		# mrnaDensity = (mrnaCounts/mrnalenMod) 
		cdsDensity = (cdsCounts/cdslenMod) 
		# utr5Density = (utr5Counts/utr5lenMod)
		# utr3Density = (utr3Counts/utr3lenMod)

		#### RPKM densities
		# mrnaDensity_rpkm = (mrnaCounts/mrnalenMod) * 1000
		cdsDensity_rpkm = (cdsCounts/cdslenMod) * 1000
		# utr5Density_rpkm = (utr5Counts/utr5lenMod) * 1000
		# utr3Density_rpkm = (utr3Counts/utr3lenMod) * 1000

		### throw out zero's
		if cdsDensity == 0:
			zeroCdsdense += 1
			continue

		if cdsDensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
			lowCdsdense += 1
			continue



		### define vector in valid CDS region, normalize by cdsDensity
		cdsSplicedCounts = exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']]
		cdsNormCounts = [rpf/cdsDensity for rpf in cdsSplicedCounts] ## just region of cds within insets for counts, sum/len == 1

		# print cdslenMod, "modified cds length"
		# print sum(cdsNormCounts), "total normalized counts, should equal length of cds"
		# print sum(cdsNormCounts)/len(cdsNormCounts) ## this should == 1

		### define vector for whole transcript, normlaized by cdsDensity
		exonNormCounts = [rpf/cdsDensity for rpf in exonsplicedcounts] ## counts normalized by cds density, using only region within insets
		# print exonNormCounts
		# print sum(exonNormCounts)
		# print sum(exonNormCounts)/len(exonNormCounts) ## should typically be less than 1, unless greater density from utrs and start/stop codons



		### for start codon metagenes - 
		if region_choice == 'start':
			### check boundaries:
			# print mRNAdf.loc[trsp]['mRNAseqs'][cdsstart:cdsstart+3] ### start codon
			
			if len(exonNormCounts[cdsstart-upstreamNTs:cdsstart+3+downstreamNTs]) < upstreamNTs+3+downstreamNTs:
				regionTooShort +=1
				continue
			else:
				totalCountedTranscripts +=1
				avggene_counts = exonNormCounts[cdsstart-upstreamNTs:cdsstart+3+downstreamNTs]
				for i in range(len(avggene_counts)): ### add these counts to the running total
					averagegene[i] += avggene_counts[i]

		if region_choice == 'stop':
			### check boundaries
			# print mRNAdf.loc[trsp]['mRNAseqs'][cdsend-3:cdsend] ### stop codon

			if len(exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs]) < upstreamNTs+3+downstreamNTs:
				regionTooShort +=1
				continue
			else:
				totalCountedTranscripts +=1
				# print upstreamNTs+3+downstreamNTs
				# print len(exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs])
				avggene_counts = exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs]
				for i in range(len(avggene_counts)): ### add these counts to the running total
					averagegene[i] += avggene_counts[i]


	averagegene_equal = [rpf/totalCountedTranscripts for rpf in averagegene] ### divide by total number of valid transcripts
	# print averagegene
	# print averagegene_equal

	positions = range(-upstreamNTs-1, downstreamNTs+2) # start or stop codon is [-1,0,1]

	df = pd.DataFrame({'position':positions, 'avg':averagegene_equal})
	df = df[['position', 'avg']]

	# print df['avg'].sum()/len(df)

	# print df 

	print "Avggene run compolete for sample %s" % samp
	print "Transcripts inculded %s" % totalCountedTranscripts
	print "Number of transcripts absent in UTRfile: %s" % noUTRentry
	print "Number of transcripts with zero CDS density: %s" % zeroCdsdense
	print "Number of transcripts below CDS density threshold: %s" % lowCdsdense 
	print "Number of transcripts too short avggene region: %s" % regionTooShort
	print "- - - - - - - -"
	print df 
	print "- - - - - - - -"



	if region_choice == 'start':
		alignpos = "1"
	elif region_choice == 'stop':
		alignpos = "2"
	else:
		print "alignpos not set!!!"

	### write csv file
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (
		rootpath, genome_name, experiment, samp)
	avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (
		fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) 
	if not os.path.exists(avggene_csv_path):	os.makedirs(avggene_csv_path)

	csv_outfile = "%s/%s_%s_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv" %( # adding cdsNorm to indicate normalized densities
		avggene_csv_path, samp, pop, threshold, minlen, maxlen, alignpos)

	df.to_csv(csv_outfile, index=False)




def main():

	print UTRfilestring
	print twobitfile 

	UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)

	for samp in samplelist:
		build_avggene_firstStop(UTRdict, threshold, samp, region_choice, inset_choice)

if __name__ == '__main__':
	main()





