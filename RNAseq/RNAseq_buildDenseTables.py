"""
This script takes input density files (binary) 
	and builds a count table for each transcript in an annotation 
	Annotation is defined by UTRfilestring:	
		a csv with mrna, cds, utr5, and utr3 sizes

Output is a csv with the total number of raw reads OR rpm values 
	mapped to each region in the csv annotation file

Design is based off of Colin's makeavggene_main.py
"""
### import dependencies
import sys
import csv
import rphelper as rph
from datetime import datetime
import argparse


def build_count_tables(trspdict, UTRdict, inset_choice, cdsDenThresh, norm_type, raw_dense_thresh, rpm_dense_thresh, outfilestring):
	"""
	This is the function that returns total counts within the following regions:
		mRNA, cds, utr5, utr3

	Length of each region is modified by insets to account for biases around start and stop codons
	3'UTR occupance is calculated to look at readthrough of stop codons, to be plotted in a downstream function
	"""
	
	print "Cds Density threshold???", cdsDenThresh


	defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 } # keep these the same as the profiling insets
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
	zeroMrnaDense = 0
	lowCdsdense = 0
	totalCountedTranscripts = 0

	for trsp in trspdict:
		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue
		# from csv: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt
		# position in list: key, 0,		1,		2,		3,		4,		5,			6,		7,		8,			9

		# define base region sizes from UTRdict
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])

		if utr5len == 0:
			zeroUtrlen +=1
			# print("transcript has zero utr5 len %s") % trsp
			# sys.exit() 
			# continue
		if utr3len == 0:
			zeroUtrlen +=1
			# continue

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
		# mrnalen = utr5len+cdslen+utr3len
		mrnalen = len(exonsplicedcounts) # independent of insets

		if utr5len == 0:
			zeroUtrlen +=1
			# print "transcript has zero utr5 len %s" % trsp
			# sys.exit() 
			# continue
		if utr3len == 0:
			zeroUtrlen +=1
			# continue

		### This is the place to add insets for counts

		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		# mrnaCounts = utr5Counts+cdsCounts+utr3Counts
		mrnaCounts = sum(exonsplicedcounts) ## want this to be independent of insets

		# if ignoring insets, better to just set insets to zeros
		# mrnaCounts = sum(exonsplicedcounts)
		# cdsCounts = sum(exonsplicedcounts[cdsstart:cdsend])
		# utr5Counts = sum(exonsplicedcounts[:cdsstart])
		# utr3Counts = sum(exonsplicedcounts[cdsend:])

		mrnaDensity = mrnaCounts/mrnalen
		cdsDensity = cdsCounts/cdslen
		try: 
			utr5Density = utr5Counts/utr5len
		except ZeroDivisionError:
			utr5Density = 0
		try: 
			utr3Density = utr3Counts/utr3len
		except ZeroDivisionError:
			utr3Density = 0
		

		# if cdsDensity == 0: # discard transcripts with no reads in the CDS
		# 	zeroCdsdense += 1
		# 	continue

		if mrnaDensity == 0: # discard transcripts with no reads in the CDS
			zeroMrnaDense += 1
			continue

		if cdsDenThresh == True: ## generally want this to be false
			if norm_type == 'raw':
				if cdsDensity < raw_dense_thresh: # must have at least one ribosome per 100 bp if raw reads are used
					lowCdsdense += 1
					continue # discard these reads
			elif norm_type == 'rpm':
				if cdsDensity < rpm_dense_thresh:
					lowCdsdense += 1
					continue # discard these reads
			else:
				print "Error, normalization is undefined"
				sys.exit()
		# we would expect this to equal ~1 for mRNA-seq
		try:
			utr3_occupancy = utr3Density/cdsDensity
		except ZeroDivisionError:
			utr3_occupancy = 0

		### use the region sizes defined by insets here instead of those in the UTRfile
		#	this is more representative of how the data was actually calculated
		# #transcript,chrom,featnum,strand,mrnalen,cdslen,utr5len,utr3len,gene_name, mrnaDensity, cdsDensity, utr5Density, utr3Density, utr3_occupancy
		transcriptOutlist = [trsp, UTRdict[trsp][0], UTRdict[trsp][1], UTRdict[trsp][2], 
							mrnalen, cdslen, utr5len, utr3len, UTRdict[trsp][7], UTRdict[trsp][8], UTRdict[trsp][9],
							mrnaCounts, cdsCounts, utr5Counts, utr3Counts,
							mrnaDensity, cdsDensity, utr5Density, utr3Density, utr3_occupancy]
		countsIDlist.append(transcriptOutlist[0]) # contains list of keys for the dictionary
		countsOutdict[trsp] = transcriptOutlist
		totalCountedTranscripts += 1
	
	print "Number of transcripts absent in UTRfile: %s" % noUTRentry
	print "Number of transcripts with zero UTR lengths: %s" % zeroUtrlen
	print "Number of transcripts with zero mRNA density: %s" % zeroMrnaDense
	print "Number of transcripts below CDS density threshold: %s" % lowCdsdense
	print "Total transcripts included in count table: %s" % totalCountedTranscripts
	
	# this is all for writing to the output text file
	fc = open(outfilestring+"_output.txt", "w") # open in appending mode
	fc.write("---Summary of DenseTables Run--- \n")
	fc.write(str(datetime.now())+"\n")
	fc.write("Number of transcripts absent in UTRfile: %s \n" % noUTRentry)
	fc.write("Number of transcripts with zero UTR lengths: %s \n" % zeroUtrlen)
	fc.write("Number of transcripts with zero CDS density: %s \n" % zeroCdsdense)
	fc.write("Number of transcripts below CDS density threshold: %s \n" % lowCdsdense)
	fc.write("Total transcripts included in count table: %s \n" % totalCountedTranscripts)
	fc.write("%s insets were used with the following values: %s \n" % (inset_choice, insets))
	fc.write("Normalization type: %s \n" % norm_type)
	if cdsDenThresh == True and norm_type == 'raw':
		fc.write("CDS density threshold on raw reads: %s \n" % str(raw_dense_thresh))
	elif cdsDenThresh == True and norm_type == 'rpm':
		fc.write("CDS density threshold on rpm: %s \n" % str(rpm_dense_thresh))
	elif cdsDenThresh == False:
		fc.write("No threshold was set on cds density \n")
	else:
		print "cdsDenThresh not set!"
		fc.close()
		sys.exit()
	fc.close()

	return countsIDlist, countsOutdict
		

def write_countTable_to_csv(countsIDlist, countsOutdict, outfilestring):
	"""
	write the output of the count tables to a csv file
	"""
	t=[]
	headers= 	['#transcript','chrom','featnum','strand','mrnalen','cdslen','utr5len','utr3len',
				'gene_name','stopcodon','stop4nt', 'mrnaCounts', 'cdsCounts', 'utr5Counts', 'utr3Counts',
				'mrnaDensity','cdsDensity','utr5Density','utr3Density','utr3_occupancy']
	t.append(headers)

	for i in countsIDlist: # position starts at 0
		newline= countsOutdict[i]
		t.append(newline)

	fa = open(outfilestring+".csv", "w")
	writer = csv.writer(fa)
	writer.writerows(t)
	fa.close()


def main():
	## using argparser to load arguments from workflow
	parser = argparse.ArgumentParser()
	parser.add_argument('--trspdictfilestring', help= 'input transcript density files')
	parser.add_argument('--UTRfilestring', help= 'UTRs file')
	parser.add_argument('--cdsDenThresh', help= 'boolean value- should density filter be used?')
	parser.add_argument('--norm_type', help= 'type of normalizaion, should be raw reads or rpm')
	parser.add_argument('--raw_dense_thresh', help= 'threshold for CDS density for raw normalization')
	parser.add_argument('--rpm_dense_thresh', help= 'threshold for CDS density for rpm normalization')
	parser.add_argument('--inset_choice', help= 'inset values to be used to avoid start and stop codon peaks')
	parser.add_argument('--outfilestring', help= 'output file name')
	# parser.add_argument
	args = parser.parse_args()

	args.cdsDenThresh = args.cdsDenThresh == 'True'

	print "ARGS: ", args.cdsDenThresh
	# print bool(args.cdsDenThresh)

	trspdict= rph.readcountsf(args.trspdictfilestring)
	UTRdict= rph.readindict(open(args.UTRfilestring, "rU"))
	countsIDlist, countsOutdict = build_count_tables(trspdict, UTRdict, args.inset_choice, args.cdsDenThresh, args.norm_type, float(args.raw_dense_thresh), float(args.rpm_dense_thresh), args.outfilestring)
	write_countTable_to_csv(countsIDlist, countsOutdict, args.outfilestring)

if __name__ == '__main__':
	# execute only if run as a script
	main()