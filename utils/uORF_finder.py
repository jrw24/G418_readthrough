"""
190402 - Jamie Wangen
This script goes through a gff3 file and identifies all of the uORFs
2 data frames are output:
	the 1st has all identified uORFs in the genome (those with zero uORFs are NOT included)
	the 2nd is a summary with the number of uORFs for each transcript

Based of riboseq_stopcodon_finder

"""

### Import dependencies


import sys
import GFF
import pandas as pd
import twobitreader
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import importlib


validChrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
			'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
			'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 
			'chrX', 'chrY', 'chrM', 'chrSinV', 'chrLUC']

parser = argparse.ArgumentParser()
parser.add_argument('--gtfInFilePrefix', help= 'the input annotation downloaded from gencode without ".gtf"')
parser.add_argument('--rootDir', help = 'root directory')
parser.add_argument('--twoBitGenome', help = 'genome file in 2bit format')

args = parser.parse_args()

### specify annotation files here
GenAnPath = '%s/genomes' % args.rootDir
GTFfile = "%s.gtf" % args.gtfInFilePrefix # input annotation file in gff3 format
# genome file in 2bit format
twobitfile = args.twoBitGenome
genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file

# output file name 
uORFtableOutfile = '%s_uORFtable.csv' % (args.gtfInFilePrefix)
uORFsummaryOutfile = '%s_uORFsummary.csv' % (args.gtfInFilePrefix)



#### functions ####

### functions imported from Colin's densbuilder
def makeGFFlist(GTFinput):
	"""
	Create a dictionary with a key for each chromosome in the GFF file
	"""
	GTFlist={}
	for chr in GTFinput:
		GTFlist[chr.id]=chr
	return GTFlist

def chrpostomrnapos(chrpos,chrom,featnum,GFFlist):
	"""
	This funciton takes a genomic query position (chrpos) and a chromosome number (ex 'chr6')
		along with the feature number (the entry in the gff file for that transcript) defined by build_utr_table()
		and the dictionary of transcript from the GFFlist

	The output is mrnapos which is the transcript relative position (position along the mRNA)
		of the original genomic query postion (chrpos)
	"""
	trsp_id= GFFlist[chrom].features[featnum].id
	trsp_strand= GFFlist[chrom].features[featnum].strand
	trsp_chromstart= int(GFFlist[chrom].features[featnum].location.start.position)  # 0-based
	trsp_chromend= int(GFFlist[chrom].features[featnum].location.end.position)
	sublist=[]

	for subfeature in GFFlist[chrom].features[featnum].sub_features:     # Make list of features
		if subfeature.type== 'exon':
			start= subfeature.location.start.position
			end= subfeature.location.end.position
			sublist.append([start,end])

	if trsp_strand== -1:    
		sublist.reverse()
	assert len(sublist)!= 0, ("transcript %s has a sublist length of zero!" % trsp_id)

	prevexonlen= 0 
	for item in sublist:
		exonstart= item[0]
		exonend= item[1]
		exonlen= exonend- exonstart

		if trsp_strand== 1:
			if chrpos>= exonstart and chrpos< exonend:      
				mrnapos= prevexonlen+ chrpos- exonstart
				return mrnapos
			else:   prevexonlen+= exonlen
		else:
			if chrpos< exonend and chrpos>= exonstart:
				mrnapos= prevexonlen+ (exonend-1)- chrpos       # Need -1 because end is not actual end, it is 1 beyond end.
				return mrnapos
			else:   prevexonlen+= exonlen 

def find_uORFs(GFFlist):
	"""
	using the same basic structure as denesbuilder_main, this function identifies all uORFs and write csv files
	"""
	### define start codon
		## could possibly change this to look at non canonical start codons
	startCodon = Seq('ATG')

	### build empty data frames, rows will be appended as function iterates over transcripts
	dfCols = ['trxname', 'symbol', 'strand', 'uORFCounter', 
			  'startPosition', 'cdsExtension',
			 'utr5len', 'cdslen', 'utr3len',
			 'uORFlen', 'uORFseq', 'uORFaa']
	uORFdf = pd.DataFrame(columns=dfCols)

	summaryCols = ['trxname', 'symbol', 'chr', 'tr_number', 
				   'strand', 'uORFCounter', 'cdsExtension']
	summarydf = pd.DataFrame(columns=summaryCols)

	####

	total_transcripts = 0 
	nonvalidchorms = 0
	nonATGstart = 0
	wrongstopcodon = 0
	validchroms = 0
	excluded_chroms = []
	included_chroms = []
	for chrom in GFFlist:
		if not chrom in validChrs:
			excluded_chroms.append(chrom)
			nonvalidchorms += 1
	#         print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		# print chrom
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
	#         print transcript
			tr_attribute_list = []
			transcriptnum+=1
			trsp_id= transcript.id # it is a number 
			trsp_strand= transcript.strand
			trsp_genename= transcript.qualifiers['Name'][0]
			trsp_chromstart= int(transcript.location.start.position)  # 0-based
			trsp_chromend= int(transcript.location.end.position) 
			transcriptlist= [0.0 for x in range(abs(trsp_chromend- trsp_chromstart))] # a list for transcript (pre-mRNA), not CDS

			exonsplicedseq= SeqIO.Seq('')
			transcriptseq= SeqIO.Seq(genome[chrom][trsp_chromstart: trsp_chromend])

			### handling transcripts with no start or stop codon:
			# startcodonmrnapos = 'absent'
			# stopcodonmrnapos = 'absent'

			startCodonMrnaList = []
			stopCodonMrnaList = []
			
			
			for item in GFFlist[chrom].features[transcriptnum].sub_features:
				if trsp_strand== 1:
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position) # not 0-based
						exonstart_feat= exonstart- trsp_chromstart
						exonend_feat= exonend- trsp_chromstart # Not 0-based, it is fine for length....next line. 
						exonsplicedseq+= transcriptseq[exonstart_feat:exonend_feat] # takes from exonstart to exonend-1
					if item.type== 'start_codon':
						startcodonpos= item.location.start.position # 0-based position
						# startcodonmrnapos=  chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)  # spliced mRNA position
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.end.position- 1 # 0-based position
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))

				if trsp_strand== -1:
					# reverse_complement() # this comes from seqIO
					transcriptseq_rev= transcriptseq.reverse_complement() 
					if item.type== 'exon': # or item.type== 'CDS':  # For yeast, use 'CDS'
						exonstart= int(item.location.start.position)  # 0-based position
						exonend= int(item.location.end.position)    # not 0-based 
						exonstart_feat= (trsp_chromend-1)- (exonend- 1)         # 0-based
						exonend_feat= (trsp_chromend-1)- exonstart      # 0-based
						exonseq= transcriptseq_rev[exonstart_feat:exonend_feat+ 1] 
						exonsplicedseq= exonseq+ exonsplicedseq
					if item.type== 'start_codon':
						startcodonpos= item.location.end.position- 1 # Need to -1 to be 0-based.
						# startcodonmrnapos= chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist)
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))
					if item.type== 'stop_codon':
						stopcodonpos= item.location.start.position  # start.position is 0-based already. 
						# stopcodonmrnapos= chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist)
						stopCodonMrnaList.append(chrpostomrnapos(stopcodonpos,chrom,transcriptnum,GFFlist))
			

			if len(startCodonMrnaList) > 0:
				# print "MORE THAN 1 START", startCodonMrnaList
				startcodonmrnapos = min(startCodonMrnaList)
			else:
				print "!!! no start codon for %s" % (trsp_id)
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)


	# 		if startcodonmrnapos == 'absent' or stopcodonmrnapos == 'absent':
	# #             print "no start of stop for trsp  %s" % trsp_id
	# 			continue
			
			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]

			if str(cdsseq[:3].upper())!= "ATG": 
				nonATGstart += 1
				continue    # ignore non-AUG start codons
			stopcodon= str(cdsseq[-3:].upper())
			# if len(utr3seq) > 0:
			# 	stop4nt = stopcodon +str(utr3seq[0].upper())
			# elif len(utr3seq) == 0: 
			# 	stop4nt = '0'
			# else:
			# 	print "there is a 3'UTR with negative length..."
			# 	sys.exit()
			if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
				wrongstopcodon += 1
				continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			#### Counting of uORFs ####
			uORFcounter = 0
			cdsExtension = 0

			for i in range(len(utr5seq)):
				### iterate over every nucleotide in the 5'UTR

				codon = utr5seq[i:i+3] # define the codon at each position
				if str(codon) == str(startCodon): # check if it is a start codon
					uORFcounter += 1

					startPosition = i
					seqIndex = i
					uORFaa = []
					uORFseq = []
					# print "found start codon at pos %s" % startPosition
					aminoAcid = codon.translate()
					uORFseq.append(str(codon))
					uORFaa.append(str(aminoAcid))
					
					while str(aminoAcid) != "*": # continue this loop until a stop codon is encoutered
						seqIndex +=3 # advance by 3 nt's each time (1 codon)
						nextCodon = utr5seq[seqIndex:seqIndex+3]
						aminoAcid = nextCodon.translate()
						if len(nextCodon) == 3: # ensure that a full codon is still present, do not want 1 or 2 nts
							uORFseq.append(str(nextCodon))
							uORFaa.append(str(aminoAcid))

						if seqIndex > len(utr5seq)-2: # if uORF continues into cds, retreive sequences from here
							# -2 is because this will not yeild a full codon (only 2 nt's)
							# print "end of UTR"
							cdsExtension = 1

							utrCdsSeq = utr5seq+cdsseq
							nextCodon = utrCdsSeq[seqIndex:seqIndex+3]
							aminoAcid = nextCodon.translate()
							uORFseq.append(str(nextCodon))
							uORFaa.append(str(aminoAcid))
							# print nextCodon, aminoAcid

							if seqIndex > len(utrCdsSeq): ## if uORF exceeds coding region, stop counting this,
								### could eventually extend to the 3'UTR if any transcript exists here
								print 'end of CDS for trsp %s' % trsp_id
								break

					uORFseqCat = "".join(uORFseq) # remove seperate list entries and concat to a string
					uORFaaCat = "".join(uORFaa)

					### save all uORF features to a list, and build into a dataframe
					uORF_features = [trsp_id, trsp_genename, trsp_strand, uORFcounter, startPosition, 
									 cdsExtension, len(utr5seq), len(cdsseq), len(utr3seq),
									 len(uORFseqCat), uORFseqCat, uORFaaCat]

					dftemp = pd.DataFrame([uORF_features], columns=dfCols) ## 
					# print dftemp

					uORFdf = pd.concat([uORFdf, dftemp], ignore_index=True)

				if i == (len(utr5seq)-1): # at the end of the 5'UTR, do this ...
					# print i
					uORFsummary = [trsp_id, trsp_genename, chrom, transcriptnum, trsp_strand, uORFcounter, cdsExtension]
					# print uORFsummary
					dfSummaryTemp = pd.DataFrame([uORFsummary], columns=summaryCols)
					# print dfSummaryTemp
					summarydf = pd.concat([summarydf, dfSummaryTemp], ignore_index=True)
		
	uORFdf.to_csv(uORFtableOutfile)
	summarydf.to_csv(uORFsummaryOutfile)
	print summarydf.head()



def main():
	GTFgen= GFF.parse(GTFfile)
	GFFlist = makeGFFlist(GTFgen)
	find_uORFs(GFFlist)

if __name__ == '__main__':
	main()




