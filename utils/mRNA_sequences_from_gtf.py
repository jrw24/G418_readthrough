"""
Script for writing sequences of spliced mRNAs from a gff3 file
	modified from stop_context_seq.py
	output a csv file with #transcript, genename, sequence

"""

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
outfilestring = '%s_mRNAseqs' % (args.gtfInFilePrefix)



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

def get_mRNA_sequence(GFFlist):

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
			
			### choose start and stop codons			
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

			mRNAseq = exonsplicedseq
			cdsseq= exonsplicedseq[startcodonmrnapos: stopcodonmrnapos+ 1] # take from startcodonmrnapos to stopcodonmrnapos
			utr5seq = exonsplicedseq[:startcodonmrnapos]
			utr3seq = exonsplicedseq[stopcodonmrnapos+1:]

			outseq = utr5seq.lower()+cdsseq.upper()+utr3seq.lower()

			
			if str(cdsseq[:3].upper())!= "ATG": 
				nonATGstart += 1
				# continue    # ignore non-AUG start codons
			
			### stopcodon is included in cdsseq, represnted by the last 3nt's
			stopcodon= str(cdsseq[-3:].upper())

			if stopcodon!= "TGA" and stopcodon!= "TAG" and stopcodon!= "TAA":   
				wrongstopcodon += 1
				# continue    # ignore weird stop codons

			# build itmes in transcript attribute list
			mRNAlen = len(exonsplicedseq)
			cdslen = len(cdsseq)
			utr5len = len(utr5seq)
			utr3len = len(utr3seq)
			assert mRNAlen == utr3len+cdslen+utr5len # check that sum of features equals mRNA length

			trsp_attr_list = [trsp_id, trsp_genename, outseq]
			ucscIDlist.append(trsp_attr_list[0])
			transcriptdict[trsp_id] = trsp_attr_list
			total_transcripts += 1

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
	headers= ['#transcript', 'genename', 'mRNAseqs']
	t.append(headers)

	for i in ucscIDlist: # position starts at 0
		newline= transcriptdict[i]
		t.append(newline)

	fa = open(outfilestring+".csv", "w")
	writer = csv.writer(fa)
	writer.writerows(t)
	fa.close()

def main():
	GTFgen = GFF.parse(GTFfile)
	GFFlist = makeGFFlist(GTFgen)
	ucscIDlist, transcriptdict = get_mRNA_sequence(GFFlist)
	write_utr_stopcodon_csvfile(ucscIDlist, transcriptdict)

if __name__ == '__main__':
	# execute only if run as a script
	main()



