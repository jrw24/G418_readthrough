__author__ = "Jamie Wangen"

"""
Building noncoding RNA database for ncRNA subtraction to use with ribosome proifling/ RNAseq experiments

	Current (190604) annotation is Gencode.v30 with hg38 genome build

	Annotation files used:
		* gencode.v30.annotation.gtf.gz - contains most rRNAs
		* gencode.v30.tRNAs.gtf.gz - this has the tRNAs, similar to tRNAdb but for hg38,
			unfortunately these tRNAs don't really have gene names in this annotation file, 
			412/420 tRNAs are present in this database of 649 tRNAs
		* rRNA is from NCBI refseq, manually downloaded fasta files
			28S: RNA28SN1 - RNA28SN5
				[NR_145822.1, NR_146148.1, NR_146154.1, NR_146118.1, NR_003287.4]
			18S: RNA18SN1 - RNA18SN5
				[NR_145820.1, NR_146146.1, NR_146152.1, NR_146119.1, NR_003286.4]
			5.8S: RNA5-8SN1 - RNA5-8SN5
				[NR_145821.1, NR_146147.1, NR_146153.1, NR_146120.1, NR_003285.3]
			5S: RNA5S1 - RNA5S17
				[NR_023363.1, NR_023364.1, NR_023365.1, NR_023366.1, NR_023367.1,
				 NR_023368.1, NR_023369.1, NR_023370.1, NR_023371.1, NR_023372.1,
				 NR_023373.1, NR_023374.1, NR_023375.1, NR_023376.1, NR_023377.1,
				 NR_023378.1, NR_023379.1]
				*** some 5S may be redundant with gencode annotation

	Classes of RNA included in Gencode annotation:
		['Mt_rRNA', 'Mt_tRNA', 'rRNA', 'miRNA', 'scRNA', 'scaRNA', 'snoRNA', 'snRNA', 'sRNA', 'vaultRNA']

"""

import sys, os
import pandas as pd 
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import twobitreader
import argparse

##### inputs #####

parser = argparse.ArgumentParser()
parser.add_argument('--gtfInFile', help= 'the input annotation downloaded from gencode without ".gtf"')
parser.add_argument('--gtfInFileTrna', help= 'gencode annotation for tRNAs')
parser.add_argument('--rRNAncbi', help = 'ncbi sequences for human rRNAs')
parser.add_argument('--rootDir', help = 'root directory')
parser.add_argument('--twoBitGenome', help = 'genome file in 2bit format')

args = parser.parse_args()


### specify annotation files here
GenAnPath = '%s/genomes' % args.rootDir

# genome file in 2bit format
twobitfile = args.twoBitGenome
genome = twobitreader.TwoBitFile(twobitfile) # create handler to open the 2bit file

### ncRNA gtf file
gtf_file_path = args.gtfInFile # input annotation file in gff3 format
gtf_outdir = GenAnPath
gtf_outfile = "%s/gencodeV30_noncodingRNAs.gtf.gz" % (gtf_outdir)
fasta_outfile = "%s/gencodeV30_noncodingRNAs.fa" % (gtf_outdir)

### tRNA gtf file
tRNA_gtf_file_path = args.gtfInFileTrna
tRNA_fasta_outfile = "%s/gencodeV30_tRNAs.fa" % (gtf_outdir)

### rRNA fasta file
rRNA_fasta_file = args.rRNAncbi

### ncRNA output fasta
ncRNA_output_fasta = "%s/gencodeV30_ncRNA_all.fa" % (gtf_outdir)


### functions

def read_in_gtf(infile_path, gtf_rows_to_skip):
	""" 
	read in a gtf file
	""" 
	df = pd.read_csv(infile_path, compression='gzip', sep='\t', dtype=str, header=None, skiprows=range(gtf_rows_to_skip))

	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	df.columns = cols

	return df


def parse_entry(tr):
	"""
	given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
	this works for default gencode gtf file formats
	"""
	tr = tr.replace('"', '')
	trl = tr.split("; ")
	trdict = OrderedDict()

	for j in trl:
		k = j.split(" ")

		if k[0] in trdict:
#             print "%s already in dict" % (k[0])
			trdict[k[0]].append(k[1])
		else:       
			trdict[k[0]]=[k[1]]
	return trdict

def parse_mod_entry(tr):
	"""
	given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
		*** Use this to parse gtf files that have modified column 8
		*** this is for files compatible with downstream GFF package
	"""
	trl = tr.split(";")
	trdict = OrderedDict()

	for j in trl:
		k = j.split("=")

		if k[0] in trdict:
#             print "%s already in dict" % (k[0])
			trdict[k[0]].append(k[1])
		else:       
			trdict[k[0]]=[k[1]]
	return trdict


def build_gene_indexes(df):
	"""
	Take an input dataframe and output an ordered dict with the index slices for every gene
	"""
	geneDict = OrderedDict()

	geneCount = 0
	previousGeneIndex = 0

	current_id=""
	current_gene=""

	for i in range(len(df)):

		if df.loc[i,'feature'] == 'gene':
			trdict = parse_entry(df.loc[i,'transcript_id'])

			curGeneID = trdict['gene_id'][0]
		
			if geneCount != 0:
				newGeneIndex = i
				geneDict[current_id] = [previousGeneIndex,newGeneIndex]
				previousGeneIndex = i
				current_id = trdict['gene_id'][0]
				geneCount += 1

			else:
				newgeneIndex = 0
				geneCount +=1
				current_id = trdict['gene_id'][0]
		if i == (len(df)-1):
			newGeneIndex = i+1
			current_id = trdict['gene_id'][0]
			geneDict[current_id] = [previousGeneIndex,newGeneIndex]
	return geneDict


def find_noncoding_genes(geneDict, df):
	"""
	3) filter geneDict on noncoding genes
	"""

	### tRNA's not part of the annotation, separate GTF file for these... 
	noncoding_types = ['Mt_rRNA', 
					   'Mt_tRNA',
					   'rRNA',
					   'miRNA',
					   'scRNA',
					   'scaRNA',
					   'snoRNA',
					   'snRNA',
					   'sRNA',
					   'vaultRNA'
					  ]

	geneDictNonCoding = OrderedDict() ### store noncoding transcripts only

	for gene in geneDict:
		tr = df.loc[geneDict[gene][0], 'transcript_id']

		trdict = parse_entry(tr)

		if trdict['gene_type'][0] in noncoding_types:
			geneDictNonCoding[gene] = geneDict[gene]

	return geneDictNonCoding


def build_transcript_indexes(geneDictCoding, df):
	"""
	take an input geneDict and find the indexes for all transcripts associated with each gene
	"""

	TR_index_dict = OrderedDict()

	for gene in geneDictCoding:

		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
	
		trPrev = -1
		trNamePrev = ""
		
		### iterate through a slice of the data frame for each gene
		### search for transcripts over that slice
		### find transcript slices
		for i in range(geneDictCoding[gene][0], geneDictCoding[gene][1]):
			if trDF.loc[i,'feature'] == 'transcript':
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				trCur = i
				trNameCur = trdict['transcript_id'][0]
				
				if trPrev != -1: # do not make an entry for the first transcript
					TR_index_dict[trNamePrev] = [trPrev, trCur]

				trPrev = trCur
				trNamePrev = trNameCur
			
			### for the final transcript
			if i == geneDictCoding[gene][1]-1:
				trdict = parse_entry(trDF.loc[i,'transcript_id'])
				TR_index_dict[trdict['transcript_id'][0]] = [trCur, i+1]
	return TR_index_dict


def convert_trsp_index(geneDictNonCoding, df, TR_index_dict):
	"""
	take input geneDict and output a single transcript for each gene
		use transcript slices defined in TR_index_dict
		for noncoding genes, take the longest transcripts

		There are ~7 snoRNAs that have multiple noncoding transcripts, but the longest of these
			always encompasses the entire shorter transcripts
	"""


	geneDictCanon = OrderedDict()
	
	for gene in geneDictNonCoding:
		trDF = df.iloc[geneDictNonCoding[gene][0]:geneDictNonCoding[gene][1]]
		trDFz = trDF.reset_index(drop=True)
		
		trCount = 0
		trDictLoc = OrderedDict()
		
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				trName = trdict['transcript_id'][0]
				trDictLoc[trName] = [trDFz.loc[i, 'chromStart'], trDFz.loc[i, 'chromEnd']]
				trCount += 1
		
		if trCount > 1:
#             print gene, "more than 1 trsp !!! \n"
			
			rangeDict = OrderedDict() ## store the ranges, and take the longest
			for key in trDictLoc:
				trRange = len(range(int(trDictLoc[key][0]),int(trDictLoc[key][1])))
				rangeDict[key] = trRange
				
			v=list(rangeDict.values())
			k=list(rangeDict.keys())
			trOut = k[v.index(max(v))]
#             print trOut
			geneDictCanon[trOut] = [gene, TR_index_dict[trOut]]
			
			

		else: ## for genes with single transcripts
			trOut = trDictLoc.keys()[0]
			geneDictCanon[trOut] = [gene, TR_index_dict[trOut]]
	return geneDictCanon


def build_df_dict(df, geneDictCanon):
	"""
	for each transcript selected in geneDictCanon, 
		take a dataframe slice from df that has all entries for the transcript

	currently keying the output dictionary on 'geneIDs', 
		could be switched to transcript id's later if desired
	"""

	outDict = OrderedDict()

	for tr in geneDictCanon:
		outDict[geneDictCanon[tr][0]] = df.iloc[geneDictCanon[tr][1][0]:geneDictCanon[tr][1][1]]

	return outDict


def edit_col8(dfIn):
	"""
	Edit the format of column 8 to remove spaces and add:
		ID= ;Name= ; for transcript
		Parent= ;Name= ; for other entires

	Sort the entries based on the value of chromStart
		This will be compaitble with downstream GTF parser from GFF package

	"""

	for i in dfIn.index:
		tr = parse_entry(dfIn.loc[i,'transcript_id'])

		### if line is a transcript, add Id= Name=
		if dfIn.loc[i,"feature"] == 'transcript':

			### make entry of new identifies to work with GFF parser 
			newline = "ID=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			### remove spaces and add semicolon to items in col 8
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			### Merge these into a single string
			outline = newline+";"+line8 
			outline = outline[0:-1] ## strip off last ';'

			### set the value of the cell to the edited line
			dfIn.at[i,"transcript_id"] = outline


		### for all other lines, add Parent= Name= 
		else:
			newline = "Parent=%s;Name=%s" % ("".join(tr['transcript_id']), "".join(tr['gene_name']))
			line8 = ";".join("=".join((str(k),str(",".join(v)))) for k,v in tr.items())
			outline = newline+";"+line8
			outline = outline[0:-1]
			dfIn.at[i,"transcript_id"] = outline

	dfIn['sort_vals'] = dfIn['chromStart'].astype(int)
	dfIn = dfIn.sort_values(by=['sort_vals'], axis=0, ascending=True)
	dfIn = dfIn.drop(['sort_vals'], axis=1)
	dfOut = dfIn.copy()

	return dfOut


def mod_last_column(outDictExclu):
	"""
	Edit the final column of each dataframe
	"""
	
	outDictMod = OrderedDict()

	for key in outDictExclu:
		k = key
		dfIn = outDictExclu[key].copy()
		dfMod = edit_col8(dfIn)
		outDictMod[k] = dfMod
	
	return outDictMod


def output_df(outdict, out_file):
	"""
	compile all modified dataframes for each transcript into a master dataframe
	build the output dataframe from the modified dictionary and write this to a file:
	"""
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	colOut = ['#chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'transcript_id']
	gtfDF = pd.DataFrame(columns=cols)

	for trsp in outdict:
		gtfDF = gtfDF.append(outdict[trsp], ignore_index=True)
		
	gtfDF.columns = colOut
	# print gtfDF.head()
	gtfDF.to_csv(out_file, compression='gzip', sep='\t', index=False)


def tr_nc_dict(dfin):
	"""
	similar to build_gene_indexes()
	use this to define transcript indexes from GTF file written for ncRNA


	"""

	tr_nc_index_dict = OrderedDict()
	
	trCount = 0
	previousTrIndex = 0

	current_id=""
	current_tr=""

	for i in range(len(dfin)):
#         print dfin.loc[i]

		if dfin.loc[i,'feature'] == 'transcript':
			trdict = parse_mod_entry(dfin.loc[i,'transcript_id'])

			curGeneID = trdict['gene_id'][0]

			if trCount != 0:
				newTrIndex = i
				tr_nc_index_dict[current_id] = [previousTrIndex,newTrIndex]
				previousTrIndex = i
				current_id = trdict['gene_id'][0]
				trCount += 1

			else:
				newTrIndex = 0
				trCount +=1
				current_id = trdict['gene_id'][0]
				
		if i == (len(dfin)-1):
			newTrIndex = i+1
			current_id = trdict['gene_id'][0]
			tr_nc_index_dict[current_id] = [previousTrIndex,newTrIndex]
			
	return tr_nc_index_dict


def check_multi_exon(tr_nc_index_dict, ncdf):
	"""
	simple function that checks for multiple exons for each transcript in tr_nc_index_dict
	"""

	for gene in tr_nc_index_dict:
	
		tempdf = ncdf.iloc[tr_nc_index_dict[gene][0]:tr_nc_index_dict[gene][1]]
		exon_count = 0
		
		for i in tempdf.index:
			if tempdf.loc[i,'feature'] == 'exon':
				exon_count += 1
	#     print exon_count
		if exon_count >1 :
			print " more than one exon for %s" % gene
			sys.exit()	# prevent writing fasta if there is multi exon transcript


def make_fasta_dict(ncdf):
	
	fasta_outdict = OrderedDict() 

	for i in ncdf.index:

		if ncdf.loc[i,'feature'] == 'transcript':
			chrom = ncdf.loc[i,'#chrom']
			chrStart = int(ncdf.loc[i,'chromStart'])
			chrEnd = int(ncdf.loc[i,'chromEnd'])
			strand = ncdf.loc[i,'strand']
			
			if strand == "+":
				chrStart = chrStart-1 ## gtf files are 1 based, convert to 0-based for python
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trdict = parse_mod_entry(ncdf.loc[i,'transcript_id'])
			
			else: # for neg strand
				chrStart = chrStart-1
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trSeq = trSeq.reverse_complement() # negative strand
				trdict = parse_mod_entry(ncdf.loc[i,'transcript_id'])

			### add output annotation line features
			trID = trdict['ID'][0]
			desc = "| "+trdict['gene_type'][0]+" | "+trdict['gene_name'][0]+ " | %s; %s; %s:%s" % (chrom, strand, chrStart, chrEnd)

			trSeqRec = SeqRecord(trSeq, id=trID, name=trdict['gene_name'][0], description=desc)
			fasta_outdict[trID] = trSeqRec
	
	return fasta_outdict


def write_output_fasta_ncRNA(fasta_outdict, fasta_outfile):
	"""
	takes a dictionary of all fasta entries and write an output fasta file
	"""

	ncFasta_iterator = (record for record in fasta_outdict.values())
	SeqIO.write(ncFasta_iterator, fasta_outfile, "fasta")


### tRNA database specific functions:

def make_tRNA_fasta_dict(tRNAdf):
	"""
	similar to make_fasta_dict, but for the tRNA database
	"""


	tRNA_fasta_outdict = OrderedDict()

	for i in tRNAdf.index:

		if tRNAdf.loc[i,'feature'] == 'tRNA':
			chrom = tRNAdf.loc[i,'#chrom']
			chrStart = int(tRNAdf.loc[i,'chromStart'])
			chrEnd = int(tRNAdf.loc[i,'chromEnd'])
			strand = tRNAdf.loc[i,'strand']
			
			if strand == "+":
				chrStart = chrStart-1 ### gtf files are 1-based, convert to 0-based
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trdict = parse_entry(tRNAdf.loc[i,'transcript_id'])
			
			else: # for neg strand
				chrStart = chrStart-1
				trSeq = SeqIO.Seq(genome[chrom][chrStart:chrEnd])
				trSeq = trSeq.reverse_complement()
				trdict = parse_entry(tRNAdf.loc[i,'transcript_id'])

			trID = "tRNA_"+trdict['gene_id'][0]
			desc = "| tRNA | "+trdict['gene_type'][0] + " | %s; %s; %s:%s" % (chrom, strand, chrStart, chrEnd)

			trSeqRec = SeqRecord(trSeq, id=trID, name=trdict['gene_name'][0], description=desc)
			tRNA_fasta_outdict[trID] = trSeqRec
	
	return tRNA_fasta_outdict
		
def write_output_fasta_tRNA(tRNA_fasta_outdict, tRNA_fasta_outfile):

	tRNA_ncFasta_iterator = (record for record in tRNA_fasta_outdict.values())
	SeqIO.write(tRNA_ncFasta_iterator, tRNA_fasta_outfile, "fasta")

def merge_ncRNA_fastas():
	"""
	command to merge all fasta files into a single fasta
	"""
	merge_cmnd = "cat %s %s %s > %s" % (rRNA_fasta_file, fasta_outfile, tRNA_fasta_outfile, ncRNA_output_fasta)
	print merge_cmnd
	os.system(merge_cmnd)

	index_cmnd = "samtools faidx %s" % (ncRNA_output_fasta)
	print index_cmnd
	os.system(index_cmnd)


def main():
	
	### make the genecode ncRNA gtf file
	df = read_in_gtf(gtf_file_path, gtf_rows_to_skip=5)
	geneDict = build_gene_indexes(df)
	geneDictNonCoding = find_noncoding_genes(geneDict, df)
	TR_index_dict = build_transcript_indexes(geneDictNonCoding, df)
	geneDictCanon = convert_trsp_index(geneDictNonCoding, df, TR_index_dict)
	outDict = build_df_dict(df, geneDictCanon)
	outDictMod = mod_last_column(outDictExclu=outDict)
	output_df(outDictMod, gtf_outfile)

	### make fasta file
	ncdf = read_in_gtf(gtf_outfile, gtf_rows_to_skip=1)
	print ncdf.head()
	tr_nc_index_dict = tr_nc_dict(ncdf)
	check_multi_exon(tr_nc_index_dict, ncdf)
	fasta_outdict = make_fasta_dict(ncdf)
	write_output_fasta_ncRNA(fasta_outdict, fasta_outfile)

	### for tRNA fasta file
	tRNAdf = read_in_gtf(tRNA_gtf_file_path, gtf_rows_to_skip=5)
	tRNA_fasta_outdict = make_tRNA_fasta_dict(tRNAdf)
	write_output_fasta_tRNA(tRNA_fasta_outdict, tRNA_fasta_outfile)

	### write the output file:
	merge_ncRNA_fastas()



if __name__ == '__main__':
	main()

