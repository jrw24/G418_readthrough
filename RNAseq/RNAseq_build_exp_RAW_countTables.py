"""
190122
Processing ribosome profiling data for input to DEseq2

Goal is to generate a count table for all transcripts in the reference file 
Counts must be in RAW reads, not normalized to RPM values

function is build from riboseq_buildDenseTables_rpkm_utr3adj
"""


import sys
import os
import csv
from datetime import datetime
import pandas as pd
import argparse
import importlib

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/RNAseq" % args.rootDir)
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/RNAseq/libsettings" % args.rootDir)


rootDir = args.rootDir
libsetName = args.libSetFile

libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)
import rphelper as rph



### function inputs

inset_choice = 'zero'
# samplelist = samplelist


### load UTRdict to be used for all samples
UTRdict= rph.readindict(open(UTRfilestring, "rU"))
# countsIDlist, countsOutdict = build_count_tables(trspdict, UTRdict, utr3adj, args.inset_choice, 
# 												bool(args.cdsDenThresh), args.norm_type, float(args.raw_dense_thresh), 
# 												float(args.rpm_dense_thresh), args.outfilestring, int(args.totreads),
# 												minUtr3len)

### sample inputs:




# def build_count_table_single(sample, UTRdict = UTRdict, inset_choice = 'default'):
def build_count_table_single(sample, UTRdict = UTRdict,):

	"""
	This is the function that returns total counts within the following regions:
		mRNA, cds, utr5, utr3

	*** Do not filter out any transcripts in output file ***

	Length of each region is modified by insets to account for biases around start and stop codons

	incorporating trspdict into this function

	### using totreads hack to get back to unnormalized, integer values without going thru the trouble of also mapping RAW A site reads

	"""


	densitystring = "Density_rpm" 
	ribosome_shift = "A"
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, sample)
	
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, sample)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()

	trspdictfilestring = '%s/%s/density5p_rnaseq/%s/%sf_' %(
		fp_assign_path, densitystring, sample, sample)

	outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store count tables
	if not os.path.exists(outfolder):   os.makedirs(outfolder)

	
	### load the counts for all transcripts into a dictionary
	trspdict= rph.readcountsf(trspdictfilestring)

	### define the insets, and select which insets to use when assigning footprint counts
	### purpose of insets is to avoid counting ribosomes around the start and stop codons
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
	
	countsIDlist = []
	countsOutdict = {}

	noUTRentry = 0
	zeroUtrlen = 0
	zeroCdsdense = 0
	totalCountedTranscripts = 0

	for trsp in trspdict:
		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue
		
		### Schema for subsetting UTRdict is below:
		# from csv: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt
		# position in list: key, 0,		1,		2,		3,		4,		5,		6,		7,		8,			9

		# define base region sizes from UTRdict
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])

		## count transcripts with zero lengths for utr's
		if utr5len == 0 or utr3len == 0:
			zeroUtrlen +=1
			# print("transcript has zero utr5 len %s") % trsp
			# sys.exit() 
		# if utr3len == 0:
		# 	zeroUtrlen +=1

		### get counts from density file
		exonsplicedcounts = trspdict[trsp]

		### set starts and ends of the coding sequence
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript %s" % trsp
			sys.exit()


		### modify region lengths using insets:
		## not sure if these modified lengths are necessarily informative here since I am not calcing densities
		utr5len = utr5len-insets['utr5Inset3']
		cdslen = cdslen-insets['cdsInset5']-insets['cdsInset3']
		utr3len = utr3len-insets['utr3Inset5']
		mrnalen = utr5len+cdslen+utr3len

		
		### This is the place to add insets for counts
		## assigning RAW counts to a given region
		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		mrnaCounts = int(utr5Counts+cdsCounts+utr3Counts)


		### 181208 - adding raw counts back to data tables, use totreads to adjust these back to the appropriate value
		RAWutr5Counts = int(utr5Counts*(totreads/1E6))
		RAWutr3Counts = int(utr3Counts*(totreads/1E6))
		RAWcdsCounts = int(cdsCounts*(totreads/1E6))
		RAWmrnaCounts = RAWutr5Counts+RAWcdsCounts+RAWutr3Counts ## nevermind, already had a calculation for RAWcounts, this is not necessary anymore...
			## although, this could be a problem that I am filtering out certain transcripts still
			## I'm not sure If I can just correct this to be unfiltered...  

		### define variables to output for final list
		transcriptOutlist = [RAWmrnaCounts, RAWcdsCounts, RAWutr5Counts, RAWutr3Counts] #adding row for utr3 adjust:
		countsIDlist.append(transcriptOutlist[0]) # contains list of keys for the dictionary
		countsOutdict[trsp] = transcriptOutlist # add this transcript to the dictionary
		totalCountedTranscripts += 1
	
	print "***********************************************************"
	print "Count table for complete for %s" % sample
	print "Number of transcripts absent in UTRfile: %s" % noUTRentry
	print "Number of transcripts with zero UTR lengths: %s" % zeroUtrlen
	print "Number of transcripts with zero CDS density: %s" % zeroCdsdense
	print "Total transcripts included in count table: %s" % totalCountedTranscripts
	print "***********************************************************"	

	#### need to think about writing up these things here, might not be strictly necessary
	# this is all for writing to the output text file
	# fc = open(outfilestring+"_output.txt", "w") # open in appending mode
	# fc.write("---Summary of DenseTables Run--- \n")
	# fc.write(str(datetime.now())+"\n")
	# fc.write("Total reads in library used for normalization: %s \n" % totreads)
	# fc.write("Number of transcripts absent in UTRfile: %s \n" % noUTRentry)
	# fc.write("Number of transcripts with zero UTR lengths: %s \n" % zeroUtrlen)
	# fc.write("Number of transcripts with zero CDS density: %s \n" % zeroCdsdense)
	# fc.write("Number of transcripts below min 3'UTR length to first inframe stop: %s \n" % tooShortUtr3AdjLenCalc)
	# fc.write("Total transcripts included in count table: %s \n" % totalCountedTranscripts)
	# fc.write("%s insets were used with the following values: %s \n" % (inset_choice, insets))
	# fc.write("Normalization type: %s \n" % norm_type)
	# if cdsDenThresh == True and norm_type == 'raw':
	# 	fc.write("CDS density threshold on raw reads: %s \n" % str(raw_dense_thresh))
	# elif cdsDenThresh == True and norm_type == 'rpm':
	# 	fc.write("CDS density threshold on rpm: %s \n" % str(rpm_dense_thresh))
	# elif cdsDenThresh == False:
	# 	fc.write("No threshold was set on cds density \n")
	# else:
	# 	print "cdsDenThresh not set!"
	# 	fc.close()
	# 	sys.exit()
	# fc.close()

	# print countsOutdict['uc001qop.2'] # gapdh

	dfout = pd.DataFrame.from_dict(countsOutdict, orient='index') #upgrade to pandas v 0.23 to use the columns keyword
	cols= [sample+'_mRNACounts', sample+'_cdsCounts', sample+'_utr5Counts', sample+'_utr3Counts']
	dfout.columns = cols
	# dfout.rename(columns=cols)
	# print dfout.head()
	df_mRNA = dfout.iloc[:,0]
	df_cds = dfout.iloc[:,1]
	df_utr5 = dfout.iloc[:,2]
	df_utr3 = dfout.iloc[:,3]

	return df_mRNA, df_cds, df_utr5, df_utr3

	# return countsIDlist, countsOutdict
	#### this still needs editing

# def 

def build_output_tables():

	### define out put directory
	outdir = "%s/FPassignment/%s/%s/analysis/countTablesRAW" % (rootpath, genome_name, experiment)
	if not os.path.exists(outdir):   os.makedirs(outdir)

	df_mRNA_list = []
	df_cds_list = []
	df_utr5_list = []
	df_utr3_list = []

	for sample in samplelist:

		df_mRNA_temp, df_cds_temp, df_utr5_temp, df_utr3_temp = build_count_table_single(sample)

		df_mRNA_list.append(df_mRNA_temp)
		df_cds_list.append(df_cds_temp)
		df_utr5_list.append(df_utr5_temp)
		df_utr3_list.append(df_utr3_temp)

	# print df_mRNA_list

	df_mRNA_out = pd.concat(df_mRNA_list, axis=1)
	df_cds_out = pd.concat(df_cds_list, axis=1)
	df_utr5_out = pd.concat(df_utr5_list, axis=1)
	df_utr3_out = pd.concat(df_utr3_list, axis=1)

	# print df_mRNA_out.head()

	df_mRNA_out.to_csv("%s/%s_mRNA_RAWcounts.csv" % (outdir, experiment), index_label = "gene_id")
	df_cds_out.to_csv("%s/%s_cds_RAWcounts.csv" % (outdir, experiment), index_label = "gene_id")
	df_utr5_out.to_csv("%s/%s_utr5_RAWcounts.csv" % (outdir, experiment), index_label = "gene_id")
	df_utr3_out.to_csv("%s/%s_utr3_RAWcounts.csv" % (outdir, experiment), index_label = "gene_id")

	print "Count tables are saved in: %s" % outdir

def main():
	build_output_tables()



if __name__ == '__main__':
	# execute only if run as a script
	main()



# def write_countTable_to_csv(countsIDlist, countsOutdict, outfilestring):
# 	"""
# 	write the output of the count tables to a csv file
# 	"""
# 	t=[]
# 	headers= 	['#transcript','chrom','featnum','strand','mrnalen','cdslen','utr5len','utr3len',
# 				'gene_name','stopcodon','stop4nt', 'mrnaCounts', 'cdsCounts', 'utr5Counts', 'utr3Counts', 
# 				'RAWmrnaCounts', 'RAWcdsCounts', 'RAWutr5Counts', 'RAWutr3Counts',
# 				'mrnaDensity_rpkm','cdsDensity_rpkm','utr5Density_rpkm','utr3Density_rpkm','utr3_occupancy',
# 				'utr3LenAdj', 'inframeStopCount', 'utr3AdjCounts', 'utr3LenAdjCalc', 'utr3AdjDensity_rpkm', 'RRTS']
# 	t.append(headers)

# 	for i in countsIDlist: # position starts at 0
# 		newline= countsOutdict[i]
# 		t.append(newline)

# 	fa = open(outfilestring+".csv", "w")
# 	writer = csv.writer(fa)
# 	writer.writerows(t)
# 	fa.close()


# def main():
# 	## using argparser to load arguments from workflow
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('--trspdictfilestring', help= 'input transcript density files')
# 	parser.add_argument('--UTRfilestring', help= 'UTRs file')
# 	parser.add_argument('--cdsDenThresh', help= 'boolean value- should density filter be used?')
# 	parser.add_argument('--norm_type', help= 'type of normalizaion, should be raw reads or rpm')
# 	parser.add_argument('--raw_dense_thresh', help= 'threshold for CDS density for raw normalization')
# 	parser.add_argument('--rpm_dense_thresh', help= 'threshold for CDS density for rpm normalization')
# 	parser.add_argument('--inset_choice', help= 'inset values to be used to avoid start and stop codon peaks')
# 	parser.add_argument('--outfilestring', help= 'output file name')
# 	parser.add_argument('--totreads', help= 'total number of reads used after raw densebuilder run')
# 	# parser.add_argument
# 	args = parser.parse_args()

# 	utr3adj = pd.read_csv(stopframe_csv, index_col=0)
# 	trspdict= rph.readcountsf(args.trspdictfilestring)
# 	UTRdict= rph.readindict(open(args.UTRfilestring, "rU"))
# 	countsIDlist, countsOutdict = build_count_tables(trspdict, UTRdict, utr3adj, args.inset_choice, 
# 													bool(args.cdsDenThresh), args.norm_type, float(args.raw_dense_thresh), 
# 													float(args.rpm_dense_thresh), args.outfilestring, int(args.totreads),
# 													minUtr3len)
# 	write_countTable_to_csv(countsIDlist, countsOutdict, args.outfilestring)

# if __name__ == '__main__':
# 	# execute only if run as a script
# 	main()