"""
Custom GTF parser to retireive single appris transcripts for each gene
	Working with gencode.v30.annotation.gtf - this is the most uptodate hg38 annotation
	Similar workflow as wangen_gtf2_tr.py

### gtf files are 1 based, not 0 based with respect to genomic coordinates

"""


import sys, os
import pandas as pd 
import time
import argparse
from collections import OrderedDict

### Function Inputs: ###

parser = argparse.ArgumentParser()
parser.add_argument('--gtfInFile', help= 'the input annotation downloaded from gencode, in Gzipped format')
parser.add_argument('--rootDir', help = 'root directory')

args = parser.parse_args()



gtf_file_path = args.gtfInFile # not sure whether I can keep this gzipped...
gtf_rows_to_skip = 5

gtf_outdir = "%s/genomes/" % args.rootDir
gtf_outfile = "%s/gencodeV30_protCode_TermStopCodon_validUTRs.gtf" % (gtf_outdir)


"""
Overall workflow:

Step
	Input
	Output

1) Read in gtf file as a data.frame using 
	read_in_gtf()
	df
2) Build gene dictionary with positions of each gene in the data.frame
	build_gene_indexes()
	geneDict -> {'geneID':[geneStartIndex:geneEndIndex]}
3) Filter geneDict on coding genes only
	find_coding_genes()
	geneDictCoding -> subset of geneDict
4) Build transcript index dictionary with position of every transcript
	build_transcript_indexes()
	TR_index_dict -> {'trspID': [trStartIndex:trEndIndex]}
5) Select single transcript per gene
	Currently using ccds filtering and then appris filtering, taking isoform with shortest UTRs from appris consensus CDS
	choose_appris_canonical()
	... under constructions
	geneDictCanon -> {'trspID':[geneID,[trStartIndex:trStopIndex]]}
6) Build dictionary with a dataframe corresponding to chosen transcript for each gene
	outDict -> {'geneID': DFslice}
7) check each transcript for overlaps with neighboring transcripts on the same strand
	add all genes with overlaps to the exclusion_overlaps dictionary
	exlusion_overlaps -> {'geneID':[start,end,strand,geneID,chrom]}
8) Remove all genes with overlaps from outDict
	outDictExclu -> subset of outDict
9) Modify column8 so that it can be parsed by downstream GTF parser
	outDictMod, changed format of column 8 from outDictExclu
		Also *** sorted all entries for each transcript by chromStart - required by downstream functions (densebuilder ect...)
10) append all dataframes back into master data.frame and write this as an output gtf file


"""


def read_in_gtf():
	"""
	1) read in the gtf file here
	"""
	df = pd.read_csv(gtf_file_path, compression = 'gzip', sep='\t', dtype=str, header=None, skiprows=range(gtf_rows_to_skip))
	cols = ['#chrom', 'source', 'feature', 'chromStart', 'chromEnd', 'score', 'strand', 'frame', 'transcript_id']
	df.columns = cols
	return df


def parse_entry(tr):
	"""
	given a single row of a gtf data.frame, parse column 8 into an ordered dictionary
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


def build_gene_indexes(df):
	"""
	2) take an input dataframe and output an ordered dict with the index slices for every gene
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

### -------- ###
			
def find_coding_genes(geneDict, df):
	"""
	3) filter geneDict on coding genes only
	"""

	geneDictCoding = OrderedDict() ### store coding transcripts only

	for gene in geneDict:
		tr = df.loc[geneDict[gene][0], 'transcript_id']
		
		trdict = parse_entry(tr)
		
		if trdict['gene_type'][0] == 'protein_coding':
			geneDictCoding[gene] = geneDict[gene]

	return geneDictCoding

def find_processed_pseudogenes(geneDict, df):
	"""
	locate all of the processes pseudogenes
	"""

	pseudoGeneDict = OrderedDict()

	for gene in geneDict:

		tr = df.loc[geneDict[gene][0], 'transcript_id']
		
		trdict = parse_entry(tr)
		
		if trdict['gene_type'][0] == 'processed_pseudogene':
			pseudoGeneDict[gene] = geneDict[gene]
	return pseudoGeneDict



def find_overlapping_loci(geneDictCoding, df):
	"""
	a huge fraction of genes have 'overlaps'
	this is not super useful
	"""

	overlapLociDict = OrderedDict()

	for gene in geneDictCoding:
		tr = df.loc[geneDictCoding[gene][0], 'transcript_id']

		trdict = parse_entry(tr)

		if 'tag' in trdict:
			if 'overlapping_locus' in trdict['tag']:
				overlapLociDict[gene] = trdict
				# overlapLociDict[trdict['transcript_id'][0]]= trdict['tag']

	print overlapLociDict
	print len(geneDictCoding)
	print len(overlapLociDict)



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
		### search for transcripts ofver that slice
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



def choose_appris_canonical(geneDictCoding, df, TR_index_dict):
	"""
	take gene and transcript dictionaries, and choose the best transcript accoding these criteria:
	
	re-working this to include all protein coding, and all ccds genes

	"""

	geneDictCanon = OrderedDict()
	geneDictChrLoc = OrderedDict() # [trxStart, trxEnd, tr_for_start, tr_for_stop]

	ap1 = 0
	ap2 = 0
	ap3 = 0
	ap4 = 0
	ap5 = 0

	ap1_alt = 0
	ap2_alt = 0

	single_iso_pri = 0
	single_iso_alt = 0
	multi_iso = 0
	total_iso = 0
	total_iso_pri = 0
	total_iso_alt = 0

	identical_trsps = 0
	noAppris = 0
	
	noCCDS = 0
	singleCCDS = 0
	multiCCDS = 0

	noProtCode = 0
	singleProtCode = 0
	multiProtCode = 0

	protCodeCount = 0
	cds_NF_count = 0
	noValidUTRsRemaining = 0
	invalidUTRtotal = 0

	### now we have all genes that are coding, with start and stop indexes:
	for gene in geneDictCoding:


		trDF = df.iloc[geneDictCoding[gene][0]:geneDictCoding[gene][1]]
		trDFz = trDF.reset_index(drop=True) # z is for zero based conversion here   

		### Now we are finding the CCDS transcript only
		### and retrieving the 'best' transcript according to appris_principal
		
		protCodeTrDict = OrderedDict()
		ccdsTrDict = OrderedDict() ### store all transcripts that have CCDS annotation
		trCount = 0
		# protCodeCount = 0
		ccdsCount = 0
		
		### find all protein coding transcripts
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				if trdict['transcript_type'][0] == 'protein_coding':

					### discard transcripts without completed coding regions:
					# if 'cds_start_NF' in trdict['tag'] or 'cds_end_NF' in trdict['tag']:
					if 'cds_start_NF' in trdict['tag'] or 'cds_end_NF' in trdict['tag'] or 'mRNA_start_NF' in trdict['tag'] or 'mRNA_end_NF' in trdict['tag']:
						cds_NF_count +=1
						# print "cds_NF for %s" % (trdict['transcript_id'])
						# print trdict['tag']
					else:
						protCodeTrDict[trdict['transcript_id'][0]] = trdict['tag'] 
						protCodeCount += 1

		# print protCodeCount
		# print protCodeTrDict.keys()


		### find all CCDS transcripts
		for i in range(len(trDFz)):
			if trDFz.loc[i, 'feature'] == 'transcript':
				trCount +=1
				### check for ccds
				tr = trDFz.loc[i, 'transcript_id']
				trdict = parse_entry(tr)
				if 'tag' in trdict:
					if 'CCDS' in trdict['tag']:
						ccdsTrDict[trdict['transcript_id'][0]]= trdict['tag']
						ccdsCount += 1

		# print ccdsCount
		# print ccdsTrDict.keys()
		# print ccdsTrDict

		### find all appris transcripts for the gene
		principal_iso_count = 0 
		alt_iso_count = 0
		
		ap1dict = OrderedDict()
		ap2dict = OrderedDict()
		ap3dict = OrderedDict()
		ap4dict = OrderedDict()
		ap5dict = OrderedDict()

		ap1dict_alt = OrderedDict()
		ap2dict_alt = OrderedDict()
		
		apDictAll_pri = OrderedDict()
		apDictAll_alt = OrderedDict()
		
		for tran in ccdsTrDict:

			if 'appris_principal_1' in ccdsTrDict[tran]:
				ap1dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 1
				ap1 += 1
				principal_iso_count += 1

			elif 'appris_principal_2' in ccdsTrDict[tran]:
				ap2dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 2
				ap2 += 1
				principal_iso_count += 1

			elif 'appris_principal_3' in ccdsTrDict[tran]:
				ap3dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 3
				ap3 += 1
				principal_iso_count += 1

			elif 'appris_principal_4' in ccdsTrDict[tran]:
				ap4dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 4
				ap4 += 1
				principal_iso_count += 1

			elif 'appris_principal_5' in ccdsTrDict[tran]:
				ap5dict[tran] = TR_index_dict[tran]
				apDictAll_pri[tran] = 5
				ap5 += 1
				principal_iso_count += 1

			elif 'appris_alternative_1' in ccdsTrDict[tran]:
				ap1dict_alt[tran] = TR_index_dict[tran]
				apDictAll_alt[tran] = 6
				ap1_alt += 1
				alt_iso_count += 1

			elif 'appris_alternative_2' in ccdsTrDict[tran]:
				ap2dict_alt[tran] = TR_index_dict[tran]
				apDictAll_alt[tran] = 7
				ap2_alt += 1
				alt_iso_count += 1

		apDictAll = OrderedDict(apDictAll_pri.items()+apDictAll_alt.items())

		# ### add a dictionary to output all appri transcripts with start and stop indexes from df
		# for trnscrpt in apDictAll:
		# 	geneDictAllAppris[trnscrpt] = [gene, TR_index_dict[trnscrpt]] 

		total_iso_pri += principal_iso_count
		total_iso_alt += alt_iso_count ### ALT work in progress
		total_iso += principal_iso_count+alt_iso_count
		
		"""	
		schema for choosing a single transcript from appris annotations
			check for single isoforms first
			then take transcripts with highest appris values first, looks like only one val allowed per gene
			if multiple transcripts with same appris values:
				check for overlaps, if none include both transcripts
					all appris principals have overlaps
				if overlaps but on different strands, include both
					all transcripts are on same strand, as expected
				if overlap and on same strand, take terminal stop codon
					if these share a stop codon, take longest 3'UTR
		"""

		###############
		## Prot Code ##
		###############


		if len(protCodeTrDict) == 0: 
			# print "no protein coding genes here!"
			noProtCode +=1 ## these should have already been filtered out
			continue



		###
		### Build dfcomp first
		### first build dfcomp for comparison features of each transcript for comparison
		dfcomp = pd.DataFrame(columns=['geneID', 'trsp', 'strand', 'trxStart', 'trxEnd',
									  'scStart', 'scEnd', 'startCodonStart', 'startCodonEnd',
									   'cdsLen', 'utr3Len', 'utr5Len', 'apri_val', 'apri_class', 'CCDS'])
		
		### define whether transcript is a primary or alternative 
		pri_vals = [1,2,3,4,5]
		alt_vals = [6,7]


		for tr in protCodeTrDict.keys():
			# print tr

			### check for appris entry
			if tr in apDictAll.keys():
				# print apDictAll[tr]
				apri_val = apDictAll[tr]

				if apri_val in pri_vals:
					apri_class = 'pri'
				elif apri_val in alt_vals:
					apri_class = 'alt'
				else:
					apri_class = 'unknown'
			else:
				apri_val = 0
				apri_class = "none"

			### check for CCDS entry
			if tr in ccdsTrDict.keys():
				ccds_present = 'yes'
			else:
				ccds_present = 'no'


			### make a df slice of each transcript to retrieve data
			tempdf = df.iloc[TR_index_dict[tr][0]:TR_index_dict[tr][1]]
			tempdf.reset_index(drop=True, inplace=True) ## set to 0-based index
			
			strand = tempdf.loc[0]['strand']
			cdsLen = 3 ## account for inclusion of stop codon, this is not included in gtf file
			utr3Len = 0
			utr5Len = 0
			
			utrVals = []
		
			for row in range(len(tempdf)):
				## define transcript start and end positions
				if tempdf.loc[row]['feature'] == 'transcript':
					trxStart = tempdf.loc[row]['chromStart']
					trxEnd = tempdf.loc[row]['chromEnd']

				## define stop_codon start and end positions
				if tempdf.loc[row]['feature'] == 'stop_codon':
					scStart = tempdf.loc[row]['chromStart']
					scEnd = tempdf.loc[row]['chromEnd']

				## define start_codon start and end positions
				if tempdf.loc[row]['feature'] == 'start_codon':
					startCodonStart = tempdf.loc[row]['chromStart']
					startCodonEnd = tempdf.loc[row]['chromEnd']
			
				## define length of coding sequence
				if tempdf.loc[row]['feature'] == 'CDS':
					cdsLen += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
				
				## build list of all UTR values
				if tempdf.loc[row]['feature'] == 'UTR':
					### 3'UTR
					if strand == '+' and tempdf.loc[row]['chromStart'] >= scStart:
						utr3Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
					if strand == '-' and tempdf.loc[row]['chromEnd'] <= scEnd:
						utr3Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
					### 5'UTR
					if strand == '+' and tempdf.loc[row]['chromEnd'] < startCodonStart:
						utr5Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
					if strand == '-' and tempdf.loc[row]['chromStart'] > startCodonEnd:
						utr5Len += int(tempdf.loc[row]['chromEnd']) - int(tempdf.loc[row]['chromStart']) + 1
			### assemble all features into a list, in order of columns
			df_tr_entry = [[gene, tr, strand, trxStart, trxEnd, scStart, scEnd, 
							startCodonStart, startCodonEnd, cdsLen, utr3Len, utr5Len, 
							apri_val, apri_class, ccds_present]]

			 ## build a temparary df to hold each transcript
			df_to_add = pd.DataFrame(df_tr_entry, columns = dfcomp.columns)
			
			## add transcript df to gene df
			dfcomp = dfcomp.append(df_to_add, ignore_index=True)

			# print dfcomp

		### output the minimum and maximum chromosome start and end positions for all coding transcripts
		# print dfcomp
		trxStartmin = dfcomp['trxStart'].min()
		trxEndmax = dfcomp['trxEnd'].max()
		tr_min = dfcomp.loc[pd.to_numeric(dfcomp['trxStart']).idxmin()]['trsp']
		tr_max = dfcomp.loc[pd.to_numeric(dfcomp['trxEnd']).idxmax()]['trsp']

		geneDictChrLoc[gene] = [trxStartmin, trxEndmax, tr_min, tr_max]

		### Check for valid UTR lengths
			## throw out transcripts with UTR lengths of zero
		pc_tr_count = len(dfcomp)

		dfcomp = dfcomp.loc[dfcomp['utr3Len']>3] # 3 is minmum UTR length here
		dfcomp = dfcomp.loc[dfcomp['utr5Len']>0]

		valUTR_tr_count = len(dfcomp)

		invalid_UTR_count = pc_tr_count - valUTR_tr_count
		invalidUTRtotal += invalid_UTR_count
		# print invalid_UTR_count, "INVAL UTR"

		### discard this gene if no valid transcripts remain
		if len(dfcomp) == 0:
			noValidUTRsRemaining +=1
			continue

		### selecting isoform: 
		if len(dfcomp) == 1:
			singleProtCode +=1
			# tr = protCodeTrDict.keys()[0] ### I think this is wrong here
			dfcomp = dfcomp.reset_index(drop=True)
			tr = dfcomp.loc[0, 'trsp']
			# print tr, "single ccds"

			### !!!! output transcript here
			geneDictCanon[tr] = [gene, TR_index_dict[tr]] 

			### build min and max genome locations for this gene
			tempdf = df.iloc[TR_index_dict[tr][0]:TR_index_dict[tr][1]]
			tempdf.reset_index(drop=True, inplace=True) ## set to 0-based index

			tempdf = tempdf.loc[tempdf['feature'] == 'transcript'] ## should only be one entry

			if len(tempdf) != 1:
				print "not single transcript for %s" % gene
				sys.exit()

			trxStart = tempdf.loc[0]['chromStart']
			trxEnd = tempdf.loc[0]['chromEnd']

			# for row in range(len(tempdf)): ## safest way to make sure that it is the whole transcript entry
			# 		## define transcript start and end positions
			# 		if tempdf.loc[row]['feature'] == 'transcript':
			# 			trxStart = tempdf.loc[row]['chromStart']
			# 			trxEnd = tempdf.loc[row]['chromEnd']

			geneDictChrLoc[gene] = [trxStart, trxEnd, tr, tr]

		elif len(protCodeTrDict) > 1:
			
			multiProtCode += 1

			### QC1 - verify that transcripts map to one strand, these have overlaping transcripts:
			strnds = dfcomp['strand'].unique()
			if len(strnds) > 1:
				print "MULTI STRAND for %s" % gene ## There are none in genecode annotation file
				sys.exit()


			# if gene == 'ENSG00000189195.13':
			# 	print dfcomp

			# ### QC2 - Check for overlaps between transcripts in dfcomp ###
			# overlaps = []
			# for j in range(len(dfcomp)):
			# 	trxspan = set(range(int(dfcomp.loc[j]['trxStart']), int(dfcomp.loc[j]['trxEnd'])))
			# 	overlaps.append(trxspan)

			# interIndex = []
			# maxcount = len(overlaps)
			# ## iteratively check each transcript for overlap with all other transcripts
			# for x in range(len(overlaps)):

			# 	y = x +1
			# 	while y < maxcount:
			# 		intersect = overlaps[x].intersection(overlaps[y])
			# 		if len(intersect) > 0:
			# 			interIndex.append([x,y])
			# 		y = y+ 1
			# if len(interIndex) == 0:
			# 	print "no overlaps for %s" % gene ### all genes in annotation have overlaps, must be a requirement for appris principals...
			# 	### take only transcripts from principal isoform
			# 	dfcomp = dfcomp.loc[dfcomp['apri_class'] == 'pri'] ## chane this to only the principal values and continue with procesing
			# 	print dfcomp
			#### -- end overlaps -- ###



			### QC3 - Take 3' most stop codon
			unique_stops = dfcomp['scEnd'].unique()

			### check for multiple stop codons, there is one case of this for ENSG00000108395.14 (TRIM37), dif exons that give same cds
			if len(unique_stops) > 1:
				# print "MULTIPLE STOP CODONS for %s" % gene 
				# print dfcomp
				# print unique_stops
				# print strnds 

				### take subset of transcripts with terminal stop codon
				if strnds == "+":
					dfcomp = dfcomp.loc[dfcomp['scEnd'] == unique_stops.max()]
				if strnds == "-":
					dfcomp = dfcomp.loc[dfcomp['scEnd'] == unique_stops.min()]
				# print dfcomp

			# if dfcomp.loc[0]['geneID'] == 'ENSG00000108395.14':
			# 	print dfcomp
			
			remaining_classes = dfcomp['apri_class'].unique()
			# print remaining_classes

			### QC4 - filter by appris classification

			if 'pri' in remaining_classes:
				dfcomp = dfcomp.loc[dfcomp['apri_class'] == 'pri']
			elif 'alt' in remaining_classes:
				dfcomp = dfcomp.loc[dfcomp['apri_class'] ==  'alt']
			elif 'none' in remaining_classes:
				dfcomp = dfcomp.loc[dfcomp['apri_class'] == 'none']
			else:
				print 'no valid classes for %s' % gene
				sys.exit()

			### filter by CCDS classification

			ccds_remaining = dfcomp['CCDS'].unique()

			if 'yes' in ccds_remaining:
				dfcomp =dfcomp.loc[dfcomp['CCDS'] == 'yes']
			elif 'no' in ccds_remaining:
				dfcomp = dfcomp.loc[dfcomp['CCDS'] == 'no']
			else:
				print 'ambiguous CCDS classification for %s' % gene
				sys.exit()


			### QC5 - take longest coding sequences

			if len(dfcomp) > 1:
				cdsLenUnique = dfcomp['cdsLen'].unique()
				# print cdsLenUniqueAlt.max()
				dfcomp = dfcomp.loc[dfcomp['cdsLen'] == cdsLenUnique.max()]

			# print dfcomp

			### >>
			###  Select transcripts with MIN 3'UTR length ###
			### >>

			### take minimum utr3 if longer than 3 nt's, unless only 1 transcript and it is that length

			min_utr3_thresh = 3

			utr3lengths = dfcomp['utr3Len'].unique() ## list of 

			if utr3lengths.min() > min_utr3_thresh: # if longer than threshold, simply take the minimum
				utr3min = utr3lengths.min()
			elif utr3lengths.min() <= min_utr3_thresh and len(utr3lengths) == 1: # if only one utr3 value and below the min, just use this
				utr3min = utr3lengths.min()
			else:
				print utr3lengths
				utr3len_thresh = [x for x in utr3lengths if x > min_utr3_thresh]
				print utr3len_thresh
				if len(utr3len_thresh) > 0:
					utr3min = min(utr3len_thresh)
				else:
					utr3min = utr3lengths.min()
				print utr3min

			# utr3min = dfcomp['utr3Len'].min() ### set minimum value for 3'UTR
			# print utr3min
			utr3valCounts = dfcomp['utr3Len'].value_counts() ### get number of occarnaces of each 3'utr length

			utr3minOccurances = utr3valCounts[utr3min]
			
			### single transcript with this 3'UTR
			if utr3minOccurances == 1:
				minInd = pd.to_numeric(dfcomp['utr3Len']).idxmin()
				outTr = dfcomp.loc[minInd,'trsp']
				
				### !!!!! output this transcript here
				geneDictCanon[outTr] = [gene, TR_index_dict[outTr]]
			
			if utr3minOccurances > 1:
				# print "multiple 3'UTRs"


				dfcomp = dfcomp.loc[dfcomp['utr3Len'] == utr3min] ### reset dfcomp to only include min utr3's
				# print dfcomp
				### next take transcript with shortest 5'UTR
				utr5min = dfcomp['utr5Len'].min()
				# if utr5min == 3:
				# 	print "no utr for %s " % gene 
				utr5valCounts = dfcomp['utr5Len'].value_counts()
				utr5minOccurances = utr5valCounts[utr5min]
				
				if utr5minOccurances == 1:
					minUtr5Ind = pd.to_numeric(dfcomp['utr5Len']).idxmin()
					outTr = dfcomp.loc[minUtr5Ind, 'trsp'] ### selected transcript
					
					### !!!!! output a transcript here
					geneDictCanon[outTr] = [gene, TR_index_dict[outTr]]
					
				if utr5minOccurances > 1:
					### some transcripts are completely identical in genecode annotation, must be historical reason
					### only difference is transcript name, and some of the annotations in last column
					# print "MULTIPLE VALID trsps for %s" % gene
					dfcomp = dfcomp.reset_index(drop=True)
					print dfcomp
					identical_trsps += 1
					outTr = dfcomp.loc[0, 'trsp'] # take the first one, since these are identical in terms of genomic features
					### !!!!! output transcript here
					geneDictCanon[outTr] = [gene, TR_index_dict[outTr]] 


	print "" # print out some summaries 
	print "Number of identical transcripts: ", identical_trsps #20 of these in the genome that fall out at the end
	print "Total number of appris isoforms: ", total_iso
	print "Number of genes with multiple appris isoforms: ", multi_iso
	# print "Number of genes with single appris isofrom: ", single_iso
	# print geneDictCanon   
	print "Count for appris principals 1/2/3/4/5/6/7 isoforms: ", ap1, ap2, ap3, ap4, ap5, ap1_alt, ap2_alt
	# print geneDictCanon
	print "CCDS genes without appris isofomrs: ", noAppris
	print "genes included == %s" % (len(geneDictCanon))
	print "genes in protCodeTrDict == %s" % (protCodeCount)
	print "cds_NF count %s" % (cds_NF_count)
	print "protein coding transcripts == %s" % (protCodeCount)
	print "genes discarded due to zero length UTRs == %s" % (noValidUTRsRemaining)
	print "total transcripts with invalid UTRs == %s" % (invalidUTRtotal)
	# print "length of geneDictAllAppris == %s" % (len(geneDictAllAppris))

	return geneDictCanon, geneDictChrLoc


def build_df_dict(df, geneDictCanon):
	"""
	for each transcript selected in geneDictCanon, take a dataframe slice from df that has all entries for the transcript

	currently keying the output dictionary on 'geneIDs', could be switched to transcript id's later if desired
	"""

	outDict = OrderedDict()

	for tr in geneDictCanon:
		outDict[geneDictCanon[tr][0]] = df.iloc[geneDictCanon[tr][1][0]:geneDictCanon[tr][1][1]]

	return outDict

# def build_df_dict_byTr(df, genDict):
# 	"""
# 	for each transcript selected in geneDictCanon, take a dataframe slice from df that has all entries for the transcript

# 	currently keying the output dictionary on 'geneIDs', could be switched to transcript id's later if desired
# 	"""

# 	outDict = OrderedDict()

# 	for tr in genDict:
# 		outDict[tr] = [genDict[tr][0], df.iloc[genDict[tr][1][0]:genDict[tr][1][1]]] # list with [geneName, dataframe]

# 	return outDict


def overlap_features(genedf):
	"""
	for a given data.frame of a transcript, extract info to check for overlaps
	Used in find_transcript_overlaps()
	output: [chrStart, chrEnd, strand, geneName, chrom]
	"""
	chrStart = genedf.loc[0,'chromStart']
	chrEnd = genedf.loc[0, 'chromEnd']
	strand = genedf.loc[0, 'strand']
	trdict = parse_entry(genedf.loc[0, 'transcript_id'])
	geneName = trdict['gene_id'][0]
	trName = trdict['transcript_id'][0]
	chrom = genedf.loc[0, '#chrom']
	outlist = [chrStart, chrEnd, strand, geneName, chrom, trName]
	return outlist

def find_transcript_overlaps(outDict, geneDictChrLoc):
	"""
	For each transcript, look at nearest neighbor on same strand and check if transcripts have overlaps
	add all transcripts with overlaps to exclusion_overlaps dictionary
	"""
	
	exclusion_overlaps = OrderedDict()
	excluded_trsp_count = 0
	gene_key_list = outDict.keys() ### ordered list of gene names in the dictionary


	trsp_entry = -1
	for gene in outDict:

		trsp_entry += 1
		genedf = outDict[gene]
		genedf = genedf.reset_index(drop=True)
		
		# chrStart = genedf.loc[0,'chromStart']
		# chrEnd = genedf.loc[0, 'chromEnd']
		# strand = genedf.loc[0, 'strand']
		# trdict = parse_entry(genedf.loc[0, 'transcript_id'])
		# geneName = trdict['gene_id'][0]
		
		overlap_feats = overlap_features(genedf)
		
		cur_strand = overlap_feats[2]
		cur_chrom = overlap_feats[4]
		
		##### DOWNSTREAM Overlaps 
		next_strand = "0"
		search_index = 1
		break_signal = 0
		
		while next_strand != cur_strand:
		
			next_tr = trsp_entry+search_index ### only for the end of the list
			if next_tr == len(gene_key_list):
				print "reached end of downstream list for %s " % gene
				break_signal = 1
				break
			
			next_gene = gene_key_list[next_tr] ### adding this
			nextdf = outDict[gene_key_list[next_tr]]
			nextdf = nextdf.reset_index(drop=True)
			next_over = overlap_features(nextdf)

			### adding genome position
			
			next_chrom = next_over[4]
			
			### check to make sure that this is still the same chr
			if next_chrom != cur_chrom:
				print "search downstream to next chr for %s" % gene
				break_signal = 1
				break

			next_strand = next_over[2]
			search_index +=1
			
		next_gene_start = geneDictChrLoc[next_gene][0]

		# if int(overlap_feats[1]) > int(next_over[0]) and break_signal !=1:
		if int(overlap_feats[1]) > int(next_gene_start) and break_signal !=1:
			# print "******"
			# print "OVERLAPING DOWNSTREAM TRSP %s, %s" % (gene, overlap_feats[5])
			# print "******"
			# print overlap_feats
			# # print next_over
			# print next_gene

			exclusion_overlaps[gene] = overlap_feats
			excluded_trsp_count += 1
			
			
		##### UPSTREAM Overlaps
		next_strand = "0"
		search_index = 1
		break_signal = 0
		
		while next_strand != cur_strand:
			next_tr = trsp_entry-search_index ### this time subtract search index
			if next_tr == -1:
				print "end of upstream list for %s" % gene
				break_signal = 1
				break
				
			next_gene = gene_key_list[next_tr] ### adding this
			nextdf = outDict[gene_key_list[next_tr]]
			nextdf = nextdf.reset_index(drop=True)
			next_over = overlap_features(nextdf)
			
			next_chrom = next_over[4]

			### check to make sure that this is still the same chr
			if next_chrom != cur_chrom:
				print "search upstream to prev chr for %s" % gene
				break_signal = 1
				break

			next_strand = next_over[2]
			search_index +=1
			
		next_gene_end = geneDictChrLoc[next_gene][1]

		# if int(overlap_feats[0]) < int(next_over[1]) and break_signal !=1:
		if int(overlap_feats[0]) < int(next_gene_end) and break_signal !=1:
			# print "******"
			# print "OVERLAPING UPSTREAM TRSP %s, %s" % (gene, overlap_feats[5])
			# print "******"
			# print overlap_feats
			# # print next_over
			# print next_gene
			
			exclusion_overlaps[gene] = overlap_feats
			excluded_trsp_count += 1

	print "Total Overlapping transcripts revised: %s" % excluded_trsp_count
	return exclusion_overlaps


def remove_overlapping_transcripts(outDict, exclusion_overlaps):
	"""
	simply remove all of the transcripts in exclusion_overlaps from outDict
	"""

	outDictExclu = outDict

	for exclu in exclusion_overlaps:
	#     print exclu
		outDictExclu.pop(exclu)
		
	print len(outDictExclu)
	return outDictExclu




def edit_col8(dfIn):
	"""
	Edit the format of column 8 to remove spaces and add:
		ID= Name= for transcript
		Parent= Name= for other entires

	Sort the entries based on the value of chromStart

		This will be compaitble with downstream GTF parser
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


def define_pseudogene_positions(pseudoGeneDict, df):

	pseudoChromDict = OrderedDict()

	for gene in pseudoGeneDict:

		# print gene
		tempdf = df.iloc[pseudoGeneDict[gene][0]:pseudoGeneDict[gene][1]]
		tempdf = tempdf.loc[tempdf['feature']=='transcript']
		tempdf = tempdf.reset_index(drop=True)

		# print tempdf
		overlap_feats = overlap_features(tempdf)
		# print overlap_feats

		over_range = overlap_feats[0:2]
		over_range = [int(x) for x in over_range]
		# over_range = map(int, over_range)
		# cur_chrom = overlap_feats[4]
		cur_chrom = overlap_feats[4]+"_"+overlap_feats[2] ## account for strandedness
		# print cur_chrom

		if cur_chrom in pseudoChromDict:
			pseudoChromDict[cur_chrom].append(over_range)
		else:
			pseudoChromDict[cur_chrom] = [over_range]

			# pseudoChromDict[cur_chrom] = pseudoChromDict[cur_chrom].append(over_range)

	# print pseudoChromDict
	return pseudoChromDict

	# print pseudoChromDict


def check_pseudo_overlap(x1, x2, y1, y2):
	"""
	check for overlaps [x1, x2] and [y1, y2]
	should only be true if these overlap
	"""
	return max(x1, y1) <= min(x2, y2) 


def find_pseudo_overlaps(outDictExclu, pseudoChromDict):

	pseudo_exclude = OrderedDict()
	pseudo_exclude_count = 0

	for gene in outDictExclu:
		# print gene
		# print outDictExclu[gene]

		genedf = outDictExclu[gene]
		genedf = genedf.reset_index(drop=True)

		for feat in genedf.index: ### only check UTR entries, not perfect, but can be modified, maybe to exons?
			if genedf.loc[feat, 'feature'] == 'UTR':
				chrStart = genedf.loc[feat,'chromStart']
				chrEnd = genedf.loc[feat, 'chromEnd']
				strand = genedf.loc[feat, 'strand']
				trdict = parse_entry(genedf.loc[feat, 'transcript_id'])
				geneName = trdict['gene_id'][0]
				trName = trdict['transcript_id'][0]
				chrom = genedf.loc[feat, '#chrom']
				outlist = [chrStart, chrEnd, strand, geneName, chrom, trName]

				overlap_feats = outlist
				# overlap_feats = overlap_features(genedf)
				
				gene_over_range = overlap_feats[0:2]
				gene_over_range = [int(x) for x in gene_over_range]
				# gene_over_range = map(int, gene_over_range)
				cur_chrom = overlap_feats[4]+"_"+overlap_feats[2]

				# print gene_over_range
				# print gene_over_range[0]
				# print gene_over_range[1]

				try: # 'chrM_+' is not in pseudoChromDict ... 
					chr_ref_list = pseudoChromDict[cur_chrom]
				except KeyError as e:
					# pseudoChromDict[cur_chrom]
					print(e.message)
					print "no key in pseudoChromDict"
					continue

				for pseudo in chr_ref_list:
					# print pseudo
					overTest = check_pseudo_overlap(gene_over_range[0], gene_over_range[1], pseudo[0], pseudo[1])

					if overTest == True:
						print "overlap found here for %s with pseudo %s" % (gene, pseudo)
						pseudo_exclude[gene] = [cur_chrom, pseudo]
						pseudo_exclude_count += 1

	print pseudo_exclude_count, "excluded pseudogene count"

	return pseudo_exclude



# 		if cur_chrom in pseudo_exclude:
# #             print "%s already in dict" % (k[0])
# 			pseudo_exclude[cur_chrom].append()
# 		else:       
# 			trdict[k[0]]=[k[1]]

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
	print gtfDF.head(5)
	gtfDF.to_csv(out_file, sep='\t', index=False)


def main():

	time_start = time.time()
	df = read_in_gtf()
	### chr22 only ###
	# df = df.iloc[2613031:2674008]
	# df = df.iloc[2213031:2674008] ## expanded to chr17
	# df = df.reset_index(drop=True)
	### ---------- ###
	print df.head()
	geneDict = build_gene_indexes(df)
	geneDictCoding = find_coding_genes(geneDict, df)
	pseudoGeneDict = find_processed_pseudogenes(geneDict, df)
	pseudoChromDict = define_pseudogene_positions(pseudoGeneDict, df)


	# print pseudoGeneDict.items()[0:3]

	### find_overlapping_loci(geneDictCoding, df)

	TR_index_dict = build_transcript_indexes(geneDictCoding, df)
	geneDictCanon, geneDictChrLoc = choose_appris_canonical(geneDictCoding, df, TR_index_dict)


	# # print geneDictChrLoc.items()[0:3]
	outDict = build_df_dict(df, geneDictCanon)
	exclusion_overlaps = find_transcript_overlaps(outDict, geneDictChrLoc)

	# # for trsp in exclusion_overlaps:
	# # 	if trsp not in exclusion_overlaps_classic:
	# # 		print trsp


	outDictExclu = remove_overlapping_transcripts(outDict, exclusion_overlaps)
	pseudo_exclude = find_pseudo_overlaps(outDictExclu, pseudoChromDict)

	outDictExcluPseudo = remove_overlapping_transcripts(outDict=outDictExclu, exclusion_overlaps=pseudo_exclude)

	outDictMod = mod_last_column(outDictExcluPseudo)
	output_df(outDictMod, gtf_outfile)
	time_end = time.time()
	print "total time == ", (time_end - time_start)


if __name__ == '__main__':
	main()