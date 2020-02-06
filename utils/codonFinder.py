"""
Script to identify all of the codons in an annotation files

Requires:
	1) UTR table with mRNA region information for all transcripts
	2) mRNA sequence csv file

"""



import sys, os
import pandas as pd
import argparse
import twobitreader
from collections import OrderedDict
from pathos.multiprocessing import ProcessingPool as Pool

parser = argparse.ArgumentParser()
parser.add_argument('--gtfInFilePrefix', help= 'the input annotation downloaded from gencode without ".gtf"')
parser.add_argument('--rootDir', help = 'root directory')
parser.add_argument('--twoBitGenome', help = 'genome file in 2bit format')
parser.add_argument('--threadNumb', help = 'number of threads')

args = parser.parse_args()

### specify annotation files here
GenAnPath = '%s/genomes' % args.rootDir
GTFfile = "%s.gtf" % args.gtfInFilePrefix # input annotation file in gff3 format
# genome file in 2bit format
twobitfile = args.twoBitGenome
genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file


### Inputs

UTRfilepath = "%s_UTRs.csv" % args.gtfInFilePrefix
mRNAfilepath = "%s_mRNAseqs.csv" % args.gtfInFilePrefix
gtfPre = args.gtfInFilePrefix
genAnName = gtfPre.split('/')[-1]
# genAnName = "%s" % args.gtfInFilePrefix

csvOutDir = "%s_codons" % args.gtfInFilePrefix
if not os.path.exists(csvOutDir):  
	os.makedirs(csvOutDir)
	
utrdf = pd.read_csv(UTRfilepath)
mRNAdf = pd.read_csv(mRNAfilepath)

codonList = [
	'AAA', 'AAC', 'AAG', 'AAT',
	'ACA', 'ACC', 'ACG', 'ACT',
	'AGA', 'AGC', 'AGG', 'AGT',
	'ATA', 'ATC', 'ATG', 'ATT',
	'CAA', 'CAC', 'CAG', 'CAT',
	'CCA', 'CCC', 'CCG', 'CCT',
	'CGA', 'CGC', 'CGG', 'CGT',
	'CTA', 'CTC', 'CTG', 'CTT',
	'GAA', 'GAC', 'GAG', 'GAT',
	'GCA', 'GCC', 'GCG', 'GCT',
	'GGA', 'GGC', 'GGG', 'GGT',
	'GTA', 'GTC', 'GTG', 'GTT',
	'TAA', 'TAC', 'TAG', 'TAT',
	'TCA', 'TCC', 'TCG', 'TCT',
	'TGA', 'TGC', 'TGG', 'TGT',
	'TTA', 'TTC', 'TTG', 'TTT'  
]

codonList1 = [
	'AAA', 'AAC', 'AAG', 'AAT',
	'ACA', 'ACC', 'ACG', 'ACT',
	'AGA', 'AGC', 'AGG', 'AGT',
	'ATA', 'ATC', 'ATG', 'ATT',
	'CAA', 'CAC', 'CAG', 'CAT',
	'CCA', 'CCC', 'CCG', 'CCT',
	'CGA', 'CGC', 'CGG', 'CGT',
	'CTA', 'CTC', 'CTG', 'CTT',
]

codonList2 = [
	'GAA', 'GAC', 'GAG', 'GAT',
	'GCA', 'GCC', 'GCG', 'GCT',
	'GGA', 'GGC', 'GGG', 'GGT',
	'GTA', 'GTC', 'GTG', 'GTT',
	'TAA', 'TAC', 'TAG', 'TAT',
	'TCA', 'TCC', 'TCG', 'TCT',
	'TGA', 'TGC', 'TGG', 'TGT',
	'TTA', 'TTC', 'TTG', 'TTT'  
]



def find_cds_seq(dfin):

	df = dfin.copy()
	
	for tr in df.index:
		mrnaseq = df.loc[tr]['mRNAseqs']

		utr5len = df.loc[tr]['5utr_len']
		cdslen = df.loc[tr]['cds_len']

		cdsStart = utr5len
		cdsEnd = utr5len + cdslen

		cdsseq = mrnaseq[cdsStart:cdsEnd]

		df.loc[tr, 'cdsSeq'] = cdsseq

	return df
	

def find_codon_positions(df, codonList):
	### iterate through codons one at a time

	for cod in codonList:
		
		print cod
		csv_outpath = "%s/%s_%s_1.csv" % (csvOutDir, genAnName, cod)
		print csv_outpath
		dfcolumns = ['headers', 'gene', 'chrom', 'trsp_num', 'cds_pos']
		tempdf = pd.DataFrame(columns = dfcolumns)
		
		for tr in df.index:
			seq = df.loc[tr, 'cdsSeq']
			
			counter = 0
			while counter < len(seq):
				codon = seq[counter:counter+3]

				if codon == cod:
					
					hdr = df.loc[tr, '#transcript']
					gene = df.loc[tr, 'genename']
					chrom = df.loc[tr, 'chrom']
					trspNum = df.loc[tr, 'featnum']
					pos = counter
					
					entryList = [[hdr, gene, chrom, trspNum, pos]]
					entryDf = pd.DataFrame(entryList, columns = dfcolumns)
					
					tempdf = tempdf.append(entryDf, ignore_index=True)


				counter +=3
				
		print tempdf
		tempdf.to_csv(csv_outpath, index=False)
		
def find_codon_positions_multi(cod, df):
	### iterate through codons one at a time

	
	print cod
	csv_outpath = "%s/%s_%s_1.csv" % (csvOutDir, genAnName, cod)
	print csv_outpath
	dfcolumns = ['headers', 'gene', 'chrom', 'trsp_num', 'cds_pos']
	tempdf = pd.DataFrame(columns = dfcolumns)
	
	for tr in df.index:
		seq = df.loc[tr, 'cdsSeq']
		
		counter = 0
		while counter < len(seq):
			codon = seq[counter:counter+3]

			if codon == cod:
				
				hdr = df.loc[tr, '#transcript']
				gene = df.loc[tr, 'genename']
				chrom = df.loc[tr, 'chrom']
				trspNum = df.loc[tr, 'featnum']
				pos = counter
				
				entryList = [[hdr, gene, chrom, trspNum, pos]]
				entryDf = pd.DataFrame(entryList, columns = dfcolumns)
				
				tempdf = tempdf.append(entryDf, ignore_index=True)


			counter +=3
			
	# print tempdf
	tempdf.to_csv(csv_outpath, index=False)


def main():
	dfref = utrdf.merge(right = mRNAdf, on='#transcript')
	df = find_cds_seq(dfref)
	dfl = [df]*len(codonList) ### replicate the dateframe as a list for the length of the input codon list
	# find_codon_positions(df, codonList)
	### split the codon list into two sets to ease the computational burder
	p = Pool(nodes=int(args.threadNumb))
	p.map(find_codon_positions_multi, codonList, dfl)
	# p.map(find_codon_positions_multi, codonList1, dfl)
	# p.map(find_codon_positions_multi, codonList2, dfl)


if __name__ == '__main__':
	main()

#     
		
	