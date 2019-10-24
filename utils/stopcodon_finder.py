"""
181106 - Jamie Wangen
This script goes through all entries in a gff3 file and finds all stop codons after the cannonical termination codon
Identify number of stop codons present in each reading
Calculate distance from cannonical stop to 1st in frame stop

Basing this off of riboseq_makeUTRtable.py

"""


### import dependencies:

import sys
import GFF
import twobitreader
import argparse
from Bio import SeqIO
import csv 

### Inputs for functions here:

# Add list of acceptable chromosomes that will be output to the table
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

# output file name (without the .csv)
outfilestring = '%s_stopcodons' % (args.gtfInFilePrefix)


### default is False
inculde_noncanon_start = False
include_noncanon_stop = False


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

def build_stopcodon_table(GFFlist, inculde_noncanon_start, include_noncanon_stop):
	"""
	This is a function to get the cds and utr sizes for an mRNA from a GFF file
	returns a list with: #transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name
	Includes most of the functions from densebuilder_main but does not return counts
	"""
	# GFFlist = GFFinput

	transcriptdict={}
	ucscIDlist = []
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
			# print chrom
			continue    # check that only valid choromosomes are used
		validchroms+=1
		included_chroms.append(chrom)
		transcriptnum= -1 # set to negative one so first transcript is == to 0
		for transcript in GFFlist[chrom].features:  # this is where the SeqFeatures are actually stored
			### testing
			# if transcript == 'uc002zku.3':	print transcript
			###
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

			### use lists to handle transcripts with multiple start and stop codons
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
						startCodonMrnaList.append(chrpostomrnapos(startcodonpos,chrom,transcriptnum,GFFlist))	# spliced mRNA position
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
				startcodonmrnapos = 0 ### adding for transcripts without start codon
			# if len(stopCodonMrnaList)
			if len(stopCodonMrnaList) > 0:
				stopcodonmrnapos = max(stopCodonMrnaList)
			else:
				print "!!! no stop codon for %s" % (trsp_id)
				stopcodonmrnapos = len(exonsplicedseq) - 3 ### leave 3nt's in "3'UTR"

			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:] ## start utr3 after the stop codon 
			
			if inculde_noncanon_start == False:
				if str(cdsseq[:3].upper())!= "ATG": 
					nonATGstart += 1
					continue    # ignore non-AUG start codons
			
			stopcodon= str(cdsseq[-3:].upper())
			if len(utr3seq) > 0:
				stop4nt = stopcodon +str(utr3seq[0].upper())
			elif len(utr3seq) == 0: 
				stop4nt = '0'
			else:
				print "there is a 3'UTR with negative length..."
				sys.exit()
			
			if include_noncanon_stop == False:
				if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
					wrongstopcodon += 1
					continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length


			###### Finding inframe stop codons ######

			### Frame zero for loop, 
				### count each codon into 3'UTR using 0-based counting
				### With zero-based counting, next stopcodon * 3 == adjusted 3'UTR length
			frameZeroTrans = utr3seq.translate()
			frameZeroPos = -1
			frameZeroStopCounter = 0
			frameZeroUtr3LenAdj = 0
			for codon in frameZeroTrans:
				frameZeroPos +=1 
				if codon == '*':
					frameZeroStopCounter +=1
				if codon == '*' and frameZeroStopCounter == 1:
					frameZeroUtr3LenAdj = frameZeroPos*3
			if frameZeroUtr3LenAdj == 0 and frameZeroStopCounter == 0:
				frameZeroUtr3LenAdj = len(utr3seq)

			### Frame +1 for loop, 
			framePlusOneTrans = utr3seq[1:].translate() # start one nucleotide into 3'UTR for +1 frameshift
			framePlusOnePos = -1
			framePlusOneStopCounter = 0
			framePlusOneUtr3LenAdj = 0
			for codon in framePlusOneTrans:
				framePlusOnePos +=1 
				if codon == '*':
					framePlusOneStopCounter +=1
				if codon == '*' and framePlusOneStopCounter == 1:
					framePlusOneUtr3LenAdj = (framePlusOnePos*3)+1
			if framePlusOneUtr3LenAdj == 0 and framePlusOneStopCounter == 0:
				# framePlusOneUtr3LenAdj = len(utr3seq[1:])
				framePlusOneUtr3LenAdj = len(utr3seq)-1
				if framePlusOneUtr3LenAdj == -1:
					framePlusOneUtr3LenAdj = 0 ## avoid negative length utr's
			### checking for bugs with 3'utr length adjustments in the +1 frame
			if utr3len !=0 and framePlusOneStopCounter == 0 and framePlusOneUtr3LenAdj == 0:
				if utr3len !=1:
					print "oh shit, lenadj is wrong in +1 frame for %s" % trsp_id
					print "utr3len: %s" % utr3len
					print "framePlusOneStopCounter: %s" % framePlusOneStopCounter
					print "framePlusOneUtr3LenAdj: %s" % framePlusOneUtr3LenAdj
					sys.exit()

			### Frame -1 for loop, 
			frameMinusOneTrans = (cdsseq[-1]+utr3seq).translate() # include last nucleotide of cds for -1 frameshift
			frameMinusOnePos = -1
			frameMinusOneStopCounter = 0
			frameMinusOneUtr3LenAdj = 0
			for codon in frameMinusOneTrans:
				frameMinusOnePos +=1 
				if codon == '*':
					frameMinusOneStopCounter +=1
				if codon == '*' and frameMinusOneStopCounter == 1:
					frameMinusOneUtr3LenAdj = (frameMinusOnePos*3)-1
			if frameMinusOneUtr3LenAdj == 0 and frameMinusOneStopCounter == 0:
				# frameMinusOneUtr3LenAdj = len(cdsseq[-1]+utr3seq)
				frameMinusOneUtr3LenAdj = len(utr3seq)+1

			trsp_attr_list = [trsp_id, chrom, transcriptnum, trsp_strand, 
							  mRNAlen, cdslen, utr5len, utr3len, trsp_genename, stopcodon, stop4nt,
							 frameZeroStopCounter, frameZeroUtr3LenAdj,
							 framePlusOneStopCounter, framePlusOneUtr3LenAdj,
							 frameMinusOneStopCounter, frameMinusOneUtr3LenAdj]

			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1
			#transcript,chrom,featnum,strand,mrna_len,cds_len,5utr_len,3utr_len,gene_name,stopcodon,stop4nt 
	print "total number of transcripts in data table: %s" % total_transcripts
	print "Number of included chromosomes chr: %s" % validchroms
	print "Number of excluded chromosomes chr: %s" % nonvalidchorms
	print "included chroms: ", included_chroms
	print "excluded chroms: ", excluded_chroms
	print "transcripts discarded due to non-AUG start codon %s" % nonATGstart
	print "transcripts discarded due to noncanonical stop codon %s" % wrongstopcodon
	return ucscIDlist, transcriptdict
		
def write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict):
	"""
	For writing dictionary to csv file
		first headers are written
		Then one line at a time is added to t, 
			This is a list that will hold all of the trsp_attr_lists 
		Finally these are written to each line of a csv file
	"""
	t=[]
	headers= ['#transcript','chrom','featnum','strand',
			  'mrna_len','cds_len','5utr_len','3utr_len','gene_name','stopcodon','stop4nt',
			  'frameZeroStopCount', 'frameZeroUtr3LenAdj',
			  'framePlusOneStopCount', 'framePlusOneUtr3LenAdj',
			  'frameMinusOneStopCount', 'frameMinusOneUtr3LenAdj']
	t.append(headers)

	for i in ucscIDlist: # position starts at 0
		newline= transcriptdict[i]
		t.append(newline)

	fa = open(outfilestring+".csv", "w")
	writer = csv.writer(fa)
	writer.writerows(t)
	fa.close()


def main():
	GTFgen= GFF.parse(GTFfile)
	GFFlist = makeGFFlist(GTFgen)
	ucscIDlist, transcriptdict = build_stopcodon_table(GFFlist, inculde_noncanon_start, include_noncanon_stop)
	write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict)

	
if __name__ == '__main__':
	# execute only if run as a script
	main()
