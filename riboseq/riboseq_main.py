__author__ = 'Jamie Wangen'

"""
Workflow for ribosome profling

Adapted from Boris Zinshteyn's and Colin Wu's pipeline 

Updated: 191009
"""

### Import dependencies
import matplotlib
matplotlib.use('Agg') ### set backend
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import argparse
import subprocess
import collections
from collections import OrderedDict
import pysam
import glob
import pandas as pd 
import seaborn as sns 
import numpy as np 
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
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

print libsetName

libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)

threadNumb = str(args.threadNumb)

### reset rootpath
# rootpath = "%s/%s" % (rootDir, rootpath)




# ### enter libsettings file for chosen library
# from riboseq_libsettings_Granzyme import *

##### Functions:

class RPexperiment(object):
	
	### instantiate all of the things required for the experiment
	def __init__(self, samp, experiment, rootpath, fastqPath, fivePrimeNTstoTrim, 
				threadNumb, genome_name, ncRNAstarGenome, starGenome, 
				linker, GTFfile, twobitfile, threshold, UTRfilestring, ftsize,
				flmin, flmax, eAmin, eAmax, eEmin, eEmax, cdsDenThresh, raw_dense_thresh, 
				rpm_dense_thresh, norm_type, inset_choice, deduplicate_reads,
				read_proc_dir, align_dir, stopcodons):
		self.samp = samp 
		self.experiment = experiment
		self.rootpath = rootpath
		self.fastqPath = fastqPath
		self.fivePrimeNTstoTrim = fivePrimeNTstoTrim
		self.threadNumb = threadNumb
		self.genome_name = genome_name
		self.ncRNAstarGenome = ncRNAstarGenome
		self.starGenome = starGenome
		self.linker = linker
		self.GTFfile = GTFfile
		self.twobitfile = twobitfile
		self.threshold = threshold
		self.UTRfilestring = UTRfilestring
		self.ftsize = ftsize
		self.flmin = flmin
		self.flmax = flmax
		self.eAmin = eAmin
		self.eAmax = eAmax
		self.eEmin = eEmin
		self.eEmax = eEmax
		self.cdsDenThresh = cdsDenThresh
		self.raw_dense_thresh = raw_dense_thresh
		self.rpm_dense_thresh = rpm_dense_thresh
		self.norm_type = norm_type
		self.inset_choice = inset_choice
		self.deduplicate_reads = deduplicate_reads
		self.read_proc_dir = read_proc_dir
		self.align_dir = align_dir
		self.stopcodons = stopcodons

	##### read_processing and alignment:
	def deduplicate_input_reads(self):
		"""
		use tally to deduplicate identical reads
		This must be done before trimming 5p ends and adapters,
			both should have randomized sequences to identify PCR duplicates
		"""
		if self.deduplicate_reads == True:

			fastqPath = self.fastqPath
			fastq_name = '%s/%s*fastq.gz' % (fastqPath, self.samp)

			dedupe_outpath = '%s/%s' % (self.rootpath, self.read_proc_dir)
			dedupe_outfile = '%s/%s_unique.fastq.gz' % (dedupe_outpath, self.samp)
			dedupe_logfile = '%s/%s_dedupe.log' % (dedupe_outpath, self.samp)
			if not os.path.exists(dedupe_outpath): os.makedirs(dedupe_outpath)

			if os.path.isfile(dedupe_outfile): 
				print "deduplication already performed"
			else:
				### dedupe command
				dedupe = 'tally -i %s -o %s --with-quality 2>>%s' % (
					fastq_name, dedupe_outfile, dedupe_logfile)
				print dedupe
				subprocess.Popen(dedupe, shell=True).wait()
				print "***** deduplication complete *****"
				print ""
		elif self.deduplicate_reads == False:

			print "No deduplication for sample: %s" % self.samp

		else:
			print "deduplication is not set!"
			sys.exit()


	def cut_5p_end(self):
		"""
		using seqtk trimfq to trim nucleotides from 5' end of reads
		these come from addition in RT primer: 'RNNN'
		check to see whether reads have been deduplicated before cutting the 5 prime end
		seqtk does not produce log files
		"""

		if self.deduplicate_reads == True:
			print "using deduplicated reads"
			trimming_out_path = '%s/%s' % (self.rootpath, self.read_proc_dir)
			trimming_out_file = '%s/%s_5ptrimmed.fastq' % (trimming_out_path, self.samp)
			trimming_out_filegz = '%s.gz' % (trimming_out_file)

			if not os.path.exists(trimming_out_path): os.makedirs(trimming_out_path) ### this should already exist, but check anyway

			fastq_name = '%s/%s_unique.fastq.gz' % (trimming_out_path, self.samp)

			if os.path.isfile(trimming_out_filegz):
				print "fiveprime end already trimmed"
			else:
				### trimming command
				fiveP_trimming = 'seqtk trimfq -b %s %s > %s' % (
					self.fivePrimeNTstoTrim, fastq_name, trimming_out_file)
				print fiveP_trimming
				subprocess.Popen(fiveP_trimming, shell=True).wait()
				### compression of trimmed file
				compression_command = 'pigz -p %s %s' % (self.threadNumb, trimming_out_file)
				print compression_command
				subprocess.Popen(compression_command, shell=True).wait()
				print "***** 5' end trimming complete *****"
				print ""

		else:
			print "reads were not deduplicated previously"
			fastqPath = self.fastqPath
			trimming_out_path = '%s/%s' % (self.rootpath, self.read_proc_dir)
			trimming_out_file = '%s/%s_5ptrimmed.fastq' % (trimming_out_path, self.samp)
			trimming_out_filegz = '%s.gz' % (trimming_out_file)

			# verify that directory exists for writing files
			if not os.path.exists(trimming_out_path): os.makedirs(trimming_out_path)
			# sample entry must be identical to start of fastq file name
				# must also end in fastq.gz
			fastq_name = '%s/%s*fastq.gz' % (fastqPath, self.samp)

			if os.path.isfile(trimming_out_filegz):
				print "fiveprime end already trimmed"
			else:
				### trimming command
				fiveP_trimming = 'seqtk trimfq -b %s %s > %s' % (self.fivePrimeNTstoTrim, fastq_name, trimming_out_file)
				print fiveP_trimming
				subprocess.Popen(fiveP_trimming, shell=True).wait()
				### compression of trimmed file
				compression_command = 'pigz -p %s %s' % (self.threadNumb, trimming_out_file)
				print compression_command
				subprocess.Popen(compression_command, shell=True).wait()


	def remove_adapter(self):
		"""
		-x specifies the 3' adapter on on the read
		-l specifies the minimum post-trimming read length
		-L specifies the maximum post-trimming read length
		-Q specifies the minimum quality score 
		"""
		trimming_out_path = '%s/%s' % (self.rootpath, self.read_proc_dir)
		trimming_out_file = '%s/%s' % (trimming_out_path, self.samp)
		trimming_out_filegz = '%s-trimmed.fastq.gz' % (trimming_out_file)
		if not os.path.exists(trimming_out_path): os.makedirs(trimming_out_path)

		fiveP_trimmed_infile = '%s_5ptrimmed.fastq.gz' % (trimming_out_file)
		if os.path.isfile(trimming_out_filegz):
			print "adapter previously removed"
		else:
			read_trimming = 'skewer -x %s -Q 10 -l 15 -L 51 --quiet -o %s --threads %s %s' % (self.linker, trimming_out_file, self.threadNumb, fiveP_trimmed_infile)
			print read_trimming
			subprocess.Popen(read_trimming, shell=True).wait()
			
			compression_command = 'pigz -p %s %s-trimmed.fastq' % (self.threadNumb, trimming_out_file)
			print compression_command
			subprocess.Popen(compression_command, shell=True).wait()


	def ncRNA_subtract_STAR(self):
		"""
		using --alignIntronMax 1 to prevent splicing
		set to 1 mismatch using --outFilterMismatchNmax 1
		"""
		ncRNA_subSTAR_outpath = '%s/%s' % (self.rootpath, self.read_proc_dir)
		ncRNA_subSTAR_outfile = '%s/%s' % (ncRNA_subSTAR_outpath, self.samp)
		if not os.path.exists(ncRNA_subSTAR_outpath): os.makedirs(ncRNA_subSTAR_outpath)
		
		ncRNA_sub_star = 'STAR \
						  --runThreadN %s \
						  --genomeDir %s \
						  --readFilesIn %s-trimmed.fastq.gz \
						  --readFilesCommand gunzip -c \
						  --outReadsUnmapped Fastx \
						  --limitBAMsortRAM 20000000000 \
						  --outSAMtype BAM SortedByCoordinate \
						  --outFilterMismatchNmax 1 \
						  --outFilterMultimapNmax 100 \
						  --alignIntronMax 1 \
						  --outWigType wiggle read1_5p \
						  --outFileNamePrefix %s' % (
			self.threadNumb, self.ncRNAstarGenome, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)

		print ncRNA_sub_star
		subprocess.Popen(ncRNA_sub_star, shell=True).wait()
		
		index_command = 'samtools index %sAligned.sortedByCoord.out.bam' % (ncRNA_subSTAR_outfile)
		subprocess.Popen(index_command, shell=True).wait()

		compression_command = 'pigz -p %s -c %sUnmapped.out.mate1 > %s_no_ncRNA.fastq.gz' % (self.threadNumb, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
		subprocess.Popen(compression_command, shell=True).wait()
		subprocess.Popen('rm %sUnmapped.out.mate1' % (ncRNA_subSTAR_outfile), shell=True).wait() 


	### STAR alignments
	def build_STAR_index(self):
		"""
		Build STAR index for genomic alignments if index does not currently exist for chosen genome
		Still need to revise this to make it more adaptable
		Not part of typical workflow
		"""
		CMMD= 'STAR \
				--runThreadN %s \
				--runMode genomeGenerate \
				--genomeSAsparseD 1 \
				--genomeDir /home/jamie/genomes/human/hg19/star \
				--genomeFastaFiles /home/jamie/genomes/human/hg19/hg19.fa \
				--sjdbGTFfile /home/jamie/genomes/human/hg19/gene_annotation/knownGene.gtf \
				--sjdbOverhang 100' % (self.threadNumb)
		print CMMD
		os.system(CMMD)


	def STAR_align(self):
		"""
		Align reads to the genome using STAR
		Make sure to allow for soft-clipping to permit alignment of reads with A or T added by superscript III
		45 seems to be the max number of supported threads for alignment
		"""
		star_infile = '%s/%s/%s_no_ncRNA.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
		star_out_path = '%s/%s/%s/STAR/%s' % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		star_out_file = '%s/%s' % (star_out_path, self.samp)
		if not os.path.exists(star_out_path): os.makedirs(star_out_path)
		
		star_command = 'STAR \
			--runThreadN %s \
			--genomeDir %s \
			--readFilesIn %s \
			--readFilesCommand gunzip -c \
			--outSAMtype BAM SortedByCoordinate \
			--alignSJDBoverhangMin 1 \
			--alignSJoverhangMin 8 \
			--outFilterType BySJout \
			--outFilterMultimapNmax 1 \
			--outWigType wiggle read1_5p \
			--outFileNamePrefix %s \
			--quantMode TranscriptomeSAM \
			--outReadsUnmapped Fastx' % (
			self.threadNumb, self.starGenome, star_infile, star_out_file)
		print star_command
		subprocess.Popen(star_command, shell=True).wait()
		
		## sort and index bam outputs
		## genome
		index_genome = 'samtools index %sAligned.sortedByCoord.out.bam' % (star_out_file)
		print index_genome
		subprocess.Popen(index_genome, shell=True).wait()
		## transcriptome
		sort_transcriptome = 'samtools sort -@ %s %sAligned.toTranscriptome.out.bam %sAligned.toTranscriptome.out.sorted' % (self.threadNumb, star_out_file, star_out_file)
		index_transcriptome = 'samtools index %sAligned.toTranscriptome.out.sorted.bam' % (star_out_file)
		print sort_transcriptome
		subprocess.Popen(sort_transcriptome, shell=True).wait()
		print index_transcriptome
		subprocess.Popen(index_transcriptome, shell=True).wait()
		subprocess.Popen('rm %sAligned.toTranscriptome.out.bam' % (star_out_file), shell = True).wait() 


	def STAR_align_Single_Random_Mulitmapper(self):
		"""
		Align reads to the genome using STAR
		Make sure to allow for soft-clipping to permit alignment of reads with A or T added by superscript III
		45 seems to be the max number of supported threads for alignment

		Here we are allowing multimapping to the genome, and then randomly selecting one of those single reads with the highest mapping score
		From Alex Dobin on google groups:
			--outFilterMultimapNmax N controls the max number of loci N a read can map to in the genome 
			--outSAMmultNmax max number of multiple alignments for a read that will be output to the SAM/BAM files
			--outMultimapperOrder Random -> random order of alignments for each multi-mapper
			--runRNGseed : if single thread is used, and a seed is set, then order will remain fixed for random outputs

		"""
		star_infile = '%s/%s/%s_no_ncRNA.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
		star_out_path = '%s/%s/%s/STAR/%s' % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		star_out_file = '%s/%s' % (star_out_path, self.samp)
		if not os.path.exists(star_out_path): os.makedirs(star_out_path)
		
		star_command = 'STAR \
						--runThreadN %s \
						--genomeDir %s \
						--readFilesIn %s \
						--readFilesCommand gunzip -c \
						--outSAMtype BAM SortedByCoordinate \
						--alignSJDBoverhangMin 1 \
						--alignSJoverhangMin 8 \
						--outFilterType BySJout \
						--outFilterMultimapNmax 200 \
						--outSAMmultNmax 1 \
						--outMultimapperOrder Random \
						--outWigType wiggle read1_5p \
						--outFileNamePrefix %s \
						--quantMode TranscriptomeSAM \
						--outReadsUnmapped Fastx' % (
			self.threadNumb, self.starGenome, star_infile, star_out_file)
		print star_command
		subprocess.Popen(star_command, shell=True).wait()
		
		## sort and index bam outputs
		## genome
		index_genome = 'samtools index %sAligned.sortedByCoord.out.bam' % (star_out_file)
		print index_genome
		subprocess.Popen(index_genome, shell=True).wait()
		## transcriptome
		sort_transcriptome = 'samtools sort -@ %s %sAligned.toTranscriptome.out.bam %sAligned.toTranscriptome.out.sorted' % (self.threadNumb, star_out_file, star_out_file)
		index_transcriptome = 'samtools index %sAligned.toTranscriptome.out.sorted.bam' % (star_out_file)
		print sort_transcriptome
		subprocess.Popen(sort_transcriptome, shell=True).wait()
		print index_transcriptome
		subprocess.Popen(index_transcriptome, shell=True).wait()
		subprocess.Popen('rm %sAligned.toTranscriptome.out.bam' % (star_out_file), shell = True).wait()  

	# def read_count_tracking(self):
	# 	"""

	# 	NEEDS TO BE FIXED

	# 	perform qc following read processing to track where every read went
	# 	Steps: input-fastq > deduplication > 5p-trimming > linker-trimming > ncRNA-subtraction > alginmment-to-genome

	# 	### 190405 - looks like this has a bug, first 3 steps are above 100% somehow
	# 	"""
	# 	fq_input_reads = 0
	# 	fq_dedup_reads = 0
	# 	fiveP_trim_reads = 0
	# 	linker_trim_reads = 0
	# 	ncRNA_sub_fq_reads = 0
	# 	genome_aligned_reads = 0
	# 	density_assigned_reads = 0 

	# 	fastqPath = self.fastqPath
	# 	fastq_input = glob.glob('%s/%s*fastq.gz' % (fastqPath, self.samp))
	# 	fastq_dedup = '%s/%s/%s_unique.fastq.gz' % (rootpath, self.read_proc_dir, self.samp)
	# 	fastq_5ptrimmed = '%s/%s/%s_5ptrimmed.fastq.gz' % (rootpath, self.read_proc_dir, self.samp)
	# 	fastq_link_trimmed = '%s/%s/%s-trimmed.fastq.gz' % (rootpath, self.read_proc_dir, self.samp)
	# 	fastq_ncRNA_dep = '%s/%s/%s_no_ncRNA.fastq.gz' % (rootpath, self.read_proc_dir, self.samp)
	# 	bam_genome_aligned = '%s/%s/%s/STAR/%s/%sAligned.sortedByCoord.out.bam' % (rootpath, self.align_dir, genome_name, self.samp, self.samp)
	# 	fp_assigned = '%s/FPassignment/%s/%s/%s/%s_FPassigned_counts.txt' % (rootpath, genome_name, experiment, self.samp, self.samp)
		
	# 	with pysam.FastxFile(fastq_input[0]) as fh1:
	# 		for entry in fh1:
	# 			fq_input_reads +=1
	# 	with pysam.FastxFile(fastq_dedup) as fh2:
	# 		for entry in fh2:
	# 			fq_dedup_reads +=1
	# 	with pysam.FastxFile(fastq_5ptrimmed) as fh3:
	# 		for entry in fh3:
	# 			fiveP_trim_reads +=1
	# 	with pysam.FastxFile(fastq_link_trimmed) as fh4:
	# 		for entry in fh4:
	# 			linker_trim_reads +=1
	# 	with pysam.FastxFile(fastq_ncRNA_dep) as fh5:
	# 		for entry in fh5:
	# 			ncRNA_sub_fq_reads +=1
	# 	with pysam.AlignmentFile(bam_genome_aligned, 'rb') as bh1:
	# 		genome_aligned_reads = bh1.count()
	# 	with open(fp_assigned, 'r') as th1:
	# 		density_assigned_reads = int(th1.read())

	# 	read_count_list = [fq_input_reads, fq_dedup_reads, fiveP_trim_reads, linker_trim_reads, ncRNA_sub_fq_reads, genome_aligned_reads, density_assigned_reads]
	# 	read_count_list_percent = [float(x) for x in read_count_list]
	# 	read_count_list_percent = [(x/fq_input_reads)*100 for x in read_count_list_percent]

	# 	### write the read counts and percents to a csv file
	# 	csv_out = "%s/%s/%s_read_accounting.csv" % (rootpath, self.read_proc_dir, self.samp)
	# 	plot_labels = ['input', 'dedupe', '5ptrim', 'linker', 'ncRNAsub', 'STAR', 'FPassign']

	# 	read_dict = OrderedDict(zip(plot_labels, read_count_list))
	# 	read_df = pd.DataFrame.from_dict(read_dict, orient='index')
	# 	read_df['percent_remaining'] = read_count_list_percent
	# 	read_df.columns = ['reads_remaining', 'percent_remaining']

	# 	read_df.to_csv(csv_out)


	# 	### set x and y values for bar plot
	# 	ys = read_count_list_percent
	# 	xs = np.arange(len(read_count_list_percent))

	# 	plt.bar(xs, ys)
	# 	plt.xticks(xs, plot_labels)
	# 	for i, v in enumerate(ys):
	# 		plt.text(i-0.3, v, "{0:.2f}".format(v), color='black', fontweight='bold')

	# 	figout = '%s/%s/%s_readLoss.pdf' % (rootpath, self.read_proc_dir, self.samp)
	# 	plt.savefig(figout, format='pdf')
	# 	plt.close()

	# 	# print read_count_list
	# 	# print read_count_list_percent


	# def ncRNA_distribution(self):
	# 	"""
	# 	plot the distribution of reads mapping to the noncoding RNA
	# 	currently counts multimapping reads more than once, will need to find a work around later if this becomes important
	# 	Requires generation of csv file with a defined class for each ncRNA
	# 	"""

	# 	### import file with annotated classses of each noncoding RNA
	# 	ncRNA_csv_file = "%s_classes.csv" % (ncRNAgenome)
	# 	ncRNA_classes_df = pd.read_csv(ncRNA_csv_file)

	# 	### use samtools to count total reads that aligned to each ncRNA, write this to a table
	# 	ncRNAalignedBAM = "%s/%s/%sAligned.sortedByCoord.out.bam" % (rootpath, self.read_proc_dir, self.samp)
	# 	idxstats_table_out = "%s/%s/%s_ncRNA_idxstats_countTable.txt" % (rootpath, self.read_proc_dir, self.samp)
	# 	samtools_cmnd = "samtools idxstats %s > %s" % (ncRNAalignedBAM, idxstats_table_out)
	# 	subprocess.Popen(samtools_cmnd, shell=True).wait()

	# 	### read in table and sum counts in each region
	# 	colnames = ['ID', 'seq_length', 'mapped_reads', 'unmapped_reads']
	# 	ncRNA_countTab = pd.read_csv(idxstats_table_out, sep='\t', header=None, names=colnames)
	# 	df_ncRNA = ncRNA_classes_df.merge(ncRNA_countTab, how='left', left_on='key', right_on='ID')
	# 	# print df_ncRNA.head()
	# 	df_group = df_ncRNA.groupby('class')['mapped_reads'].sum()

	# 	### read in read accounting files:
	# 	incsv = "%s/%s/%s_read_accounting.csv" % (rootpath, self.read_proc_dir, self.samp)
	# 	indf = pd.read_csv(incsv, index_col = 0)

	# 	reads_df = pd.DataFrame(df_group)
	# 	reads_df.loc['unassigned'] = indf.loc['ncRNAsub']['reads_remaining']-indf.loc['FPassign']['reads_remaining']
	# 	reads_df.loc['assigned'] = indf.loc['FPassign']['reads_remaining']
	# 	reads_df.loc['total'] = indf.loc['linker']['reads_remaining']
	# 	reads_df['percent'] = [(x/reads_df.loc['total']['mapped_reads']*100) for x in reads_df['mapped_reads']]

	# 	df_plt = reads_df
	# 	df_plt.drop(df_plt.tail(1).index,inplace=True)

	# 	df_plt['percent'].plot.pie(startangle=0, autopct='%1.1f%%', figsize = (5,5))
	# 	plotout = '%s/%s/%s_read_accounting_piePlot.pdf' % (rootpath, self.read_proc_dir, self.samp)
	# 	plt.savefig(plotout, format='pdf')
	# 	plt.close()


	### Read size distribution

	def read_size_dist_single(self):
		"""
		updated readsize distribution for a single library
		relies on size distribution from densebuilder -
			this is more accurate as it only includes transcripts that are actually included in profiling data
			Currently (181208) not completely sure how this is handling read size calculations after soft clipping
		"""
		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		bamfile_prefix = '%s/%s_star_default/%s' % (fp_assign_path, self.samp, self.samp)
		df = pd.DataFrame(columns = ['count'])
		# print df

		for read_size in ftsize:
			inbamfile = '%s_%s_match.sorted.bam' % (bamfile_prefix, read_size)
			bamfile = pysam.AlignmentFile(inbamfile, 'rb')
			read_count = bamfile.count()
			df.loc[read_size] = [read_count]

		total_reads = df['count'].sum()
		colname = "percent"
		df[colname] = (df['count'] / total_reads)*100

		plot_outfile = '%s/%s_readsizedist_single.pdf' % (fp_assign_path, self.samp)
		# fig, ax = plt.subplots()
		df.plot.line(y='percent', style='.-', legend=False) #ax=ax,
		plt.xlim(xmin = min(ftsize), xmax = max(ftsize))
		plt.xticks(ftsize)
		plt.savefig(plot_outfile, format='pdf')
		plt.close()


	##### Densebuilder #####
		## build separate density files for each individual read length in ftsize with densebuilder_multi()
		## build density for group of read sizes with shift to A or P site using densebuilder_riboshift()
			# creates density files for full length (fl), empty A site (eA), empty E site (eE), and all lengths (al == fl+eA+eE)
			# avggene_multi(fiveprime, start) and get_riboshift() must run before this function
		# using ribosome_shift to define position of ribosome where density is assinged, 0shift for no adjustment of reads

	def densebuilder_multi_raw(self, fiveorthreeprime, ribosome_site="0"):
		"""
		Function for building footprint density according to 5' end of reads
		Setting this up for multiprocessing using ftsize inputs for iterator
		wrapper function for densebuilder_5p_onesize()
		"""
		### ensure that output directory exists
		assignment = fiveorthreeprime # "5" or "3" for fiveprime or threeprime alignment
		ribosome_shift = ribosome_site
		bam_inpath = "%s/%s/%s/STAR/%s" % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		bamfileinput = "%s/%sAligned.sortedByCoord.out.bam" % (bam_inpath, self.samp)
		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
		if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

		read_length_list = self.ftsize

		densitystring = "DensityUnnormalized"
		totreads = 1E6
		
		def densebuilder_onesize(readin):
			"""
			build dense for one size fragment
			"""
			readsize = readin 
			
			## create a dictionary for input read size with no offset
			riboshiftdict = "{"+str(readsize)+':[0]'+"}"
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.samp, readsize)
			
			## choose where to write density files:
			densityfilepath = '%s/%s/density%sp_%sshift_%s/%s_%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, readsize, self.samp, readsize)
			densityfileout = '%s/%s_%sf_' % (densityfilepath, self.samp, readsize)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
				rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)
			
			### sorting and indexing output bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.samp, readsize)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.samp, readsize)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)
		
		p = Pool(nodes=int(self.threadNumb))
		p.map(densebuilder_onesize, read_length_list)

		#######
		"""
		now we want to calculate the total number of reads from the bamfiles in ftsize
		use pysam to calculate the number of reads present for all read lengths
		"""
		readcounter = 0 
		for readlen in read_length_list:
			bamfileinput = '%s/%s_%s_match.sorted.bam' %(bamfileoutpath, self.samp, readlen)
			bamfile = pysam.AlignmentFile(bamfileinput, 'rb')
			readnumber = bamfile.count()
			readcounter += readnumber
		
		readcount_outfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.samp)
		readcountout = open(readcount_outfile, "w")
		readcountout.write(str(readcounter))
		readcountout.close()


	def densebuilder_multi_norm(self, fiveorthreeprime, ribosome_site="0"):
		"""
		Function for building footprint density according to 5' end of reads, with normalization for total valid reads
		Setting this up for multiprocessing using ftsize inputs for iterator
		wrapper function for densebuilder_5p_onesize()
		"""
		### set input variables, use 5' alignment by default
		assignment = fiveorthreeprime # "5" or "3" for fiveprime or threeprime alignment
		ribosome_shift = ribosome_site
		bam_inpath = "%s/%s/%s/STAR/%s" % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		bamfileinput = "%s/%sAligned.sortedByCoord.out.bam" % (bam_inpath, self.samp)
		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
		if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

		densitystring = "Density_rpm"
		totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.samp)
		totreadcountf = open(totreads_countfile, "r")
		totreads = int(totreadcountf.read())
		totreadcountf.close()

		read_length_list = self.ftsize
		
		def densebuilder_onesize(readin):
			"""
			build dense for one size fragment
			"""
			readsize = readin 
			
			## create a dictionary for input read size with no offset
			riboshiftdict = "{"+str(readsize)+':[0]'+"}"
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.samp, readsize)
			
			## choose where to write density files:
			densityfilepath = '%s/%s/density%sp_%sshift_%s/%s_%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, readsize, self.samp, readsize)
			densityfileout = '%s/%s_%sf_' % (densityfilepath, self.samp, readsize)

			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
				rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)
			
			### sorting and indexing output bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.samp, readsize)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.samp, readsize)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)
		
		p = Pool(nodes=int(self.threadNumb))
		p.map(densebuilder_onesize, read_length_list)


	def densebuilder_riboshift(self, ribosome_site):
		"""
		Function to map densities for ranges of read lengths to single sites of the ribosome
		Normally map to the A site of the ribosome
		"""

		assignment = "5" # should be five prime assignment to use shift dicts
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
		
		bam_inpath = "%s/%s/%s/STAR/%s" % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		bamfileinput = "%s/%sAligned.sortedByCoord.out.bam" % (bam_inpath, self.samp)
		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
		if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

		densitystring = "Density_rpm"
		totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.samp)
		totreadcountf = open(totreads_countfile, "r")
		totreads = int(totreadcountf.read())
		totreadcountf.close()

		## Set paths to riboshift dicts for each species of interest
		fl_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eA_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eE_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)


		fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
		fl_ribo_shift.columns = ['riboshift']
		eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
		eA_ribo_shift.columns = ['riboshift']
		eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
		eE_ribo_shift.columns = ['riboshift']


		fl_sizes = fl_ribo_shift.index.tolist()
		eA_sizes = eA_ribo_shift.index.tolist()
		eE_sizes = eE_ribo_shift.index.tolist()

		fl_shiftdict = collections.OrderedDict()
		eA_shiftdict = collections.OrderedDict()
		eE_shiftdict = collections.OrderedDict()

		for i in fl_sizes:
			fl_shiftdict.update( {str(i):'['+str(fl_ribo_shift.loc[i,'riboshift'])+']'})
		for i in eA_sizes:
			eA_shiftdict.update( {str(i):'['+str(eA_ribo_shift.loc[i,'riboshift'])+']'})
		for i in eE_sizes:
			eE_shiftdict.update( {str(i):'['+str(eE_ribo_shift.loc[i,'riboshift'])+']'})

		allLenghts_shiftdict = collections.OrderedDict()
		allLenghts_shiftdict.update(eE_shiftdict)
		allLenghts_shiftdict.update(eA_shiftdict)
		allLenghts_shiftdict.update(fl_shiftdict)

		def build_fl_dense(fldict):

			readsize_dict = fldict
			pop = "fl"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# fl_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'} ## manually set read lengths
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.samp, pop)

			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.samp, minlen, maxlen)
			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, self.samp, minlen, maxlen)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
				rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)
			
			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.samp, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.samp, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

		def build_eA_dense(eAdict):

			readsize_dict = eAdict
			pop = "eA"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# eA_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.samp, pop)

			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.samp, minlen, maxlen)
			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, self.samp, minlen, maxlen)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
				rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)

			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.samp, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.samp, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

		def build_eE_dense(eEdict):

			readsize_dict = eEdict
			pop = "eE"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# eE_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.samp, pop)
			
			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.samp, minlen, maxlen)
			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, self.samp, minlen, maxlen)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
				rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)

			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.samp, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.samp, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

		def build_allLengths_dense(alldict):
			
			readsize_dict = alldict
			pop = "aL"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# aL_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.samp, pop)

			densityfilepath = '%s/%s/density%sp_%sshift_%s/%s_%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, pop, self.samp, pop)
			densityfileout = '%s/%s_%sf_' % (densityfilepath, self.samp, pop)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
				rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)

			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.samp, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.samp, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)


		p1 = multiprocessing.Process(target=build_fl_dense, args=(fl_shiftdict,))
		p1.start()
		p2 = multiprocessing.Process(target=build_eA_dense, args=(eA_shiftdict,))
		p2.start()
		p3 = multiprocessing.Process(target=build_eE_dense, args=(eE_shiftdict,))
		p3.start()
		p4 = multiprocessing.Process(target=build_allLengths_dense, args=(allLenghts_shiftdict,))
		p4.start()
		p1.join()
		p2.join()
		p3.join()
		p4.join()


	### Avggene functions:
	def avggene_multi(self, alignposition, fiveorthreeprime, normalizaiton, ribosome_site="0"):
		"""
		Create avggenes for all readlengths defined by ftsize
		Run all avggene scripts in parallel using pathos

		"""
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
		read_length_list = self.ftsize

		# this is for normalizaing within a specific transcript I think... need to check on this...
		# From Colin: equalweight = '0'    # 1 gives equal weight to all genes.
		if norm == 'uneq':
			equalweight = '0'
		if norm == 'eq':
			equalweight = '1'
		
		filtermodule = str(0)
		exclusionmodule = str(0)
		thresh = str(0)  # rpkm threshold
		if alignpos == '1':
			regionlength5 = '50'
			regionlength3 = '150'
			cdslength = regionlength3
		if alignpos == '2':
			regionlength5 = '150'
			regionlength3 = '100'
			cdslength = regionlength5
		if self.norm_type == "raw":
			densitystring = "DensityUnnormalized"
		elif self.norm_type == "rpm":
			densitystring = "Density_rpm"
		else:
			print "Normalization is not set for densebuilder!"
			sys.exit()

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)

		outfolder =  "%s/avggene%s_ORF%s_%sshift_%s%s%s" % (fp_assign_path, alignpos, self.norm_type, ribosome_shift, assignment, norm, cdslength)
		if not os.path.exists(outfolder):   os.makedirs(outfolder)

		def avggene_onesize(readin):
			"""
			Calculate average gene for one read size
			"""
			readsize = str(readin) # set read size
			trspdictfilestring = '%s/%s/density%sp_%sshift_%s/%s_%sf/%s_%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, readsize, self.samp, readsize, self.samp, readsize)
			
			outfilebase = "%s/%s_%sf" % (outfolder, self.samp, readsize)
			print ""
			commandstring = 'python2 %s/riboseq/makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(
				rootDir,regionlength5,regionlength3,trspdictfilestring,self.UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
			print commandstring
			os.system(commandstring)

		p = Pool(nodes=int(self.threadNumb))
		p.map(avggene_onesize, read_length_list)


	### Avgerage gene function for riboshifted density files:
	def avggene_riboshift(self, alignposition, ribosome_site, fiveorthreeprime='5', normalizaiton='uneq', threshold='0'):
		"""
		Create avggenes for all readlengths defined by ftsize
		Run all avggene scripts in parallel using pathos
		"""
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # should be 5' mapped at this point
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
		
		if norm == 'uneq':
			equalweight = '0'
		if norm == 'eq':
			equalweight = '1' 	#equalweight = '0'    # 1 gives equal weight to all genes.

		filtermodule = str(0)
		exclusionmodule = str(0)
		thresh = threshold # rpkm threshold
		
		if alignpos == '1':
			regionlength5 = '50'
			regionlength3 = '150'
			cdslength = regionlength3
		if alignpos == '2':
			regionlength5 = '150'
			regionlength3 = '100'
			cdslength = regionlength5
		if self.norm_type == "raw":
			densitystring = "DensityUnnormalized"
		elif self.norm_type == "rpm":
			densitystring = "Density_rpm"
		else:
			print "Normalization is not set for densebuilder!"
			sys.exit()

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		

		outfolder =  "%s/avggene%s_ORF%s_%sshift_%s%s%s" % (fp_assign_path, alignpos, self.norm_type, ribosome_shift, assignment, norm, cdslength) # norm should be 'uneq' for now
		if not os.path.exists(outfolder):   os.makedirs(outfolder)

		### define read sizes for fl, eA, and eE - need minlen and maxlen
		fl_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eA_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eE_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)


		fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
		fl_ribo_shift.columns = ['riboshift']
		eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
		eA_ribo_shift.columns = ['riboshift']
		eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
		eE_ribo_shift.columns = ['riboshift']

		# print eE_ribo_shift.index.tolist()
		fl_sizes = fl_ribo_shift.index.tolist()
		eA_sizes = eA_ribo_shift.index.tolist()
		eE_sizes = eE_ribo_shift.index.tolist()

		minlen_fl = fl_sizes[0]
		maxlen_fl = fl_sizes[-1]
		minlen_eA = eA_sizes[0]
		maxlen_eA = eA_sizes[-1]
		minlen_eE = eE_sizes[0]
		maxlen_eE = eE_sizes[-1]
		#####

		def avggene_fl():
			"""
			Calculate average gene for one read size
			"""
			# readsize = str(readin)
			# change this to read in fl reads
			pop = "fl"
			minlen = minlen_fl
			maxlen = maxlen_fl
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.samp, minlen, maxlen, self.samp, minlen, maxlen)
			outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, self.samp, pop, thresh, minlen, maxlen)
			print ""
			commandstring = 'python2 %s/riboseq/makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(
				rootDir,regionlength5,regionlength3,trspdictfilestring,self.UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)

			print commandstring
			os.system(commandstring)
		
		def avggene_eA():
			"""
			Calculate average gene for one read size
			"""
			# readsize = str(readin)
			# change this to read in fl reads
			pop = "eA"
			minlen = minlen_eA
			maxlen = maxlen_eA
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.samp, minlen, maxlen, self.samp, minlen, maxlen)
			outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, self.samp, pop, thresh, minlen, maxlen)
			print ""
			commandstring = 'python2 %s/riboseq/makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(
				rootDir,regionlength5,regionlength3,trspdictfilestring,self.UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)

			print commandstring
			os.system(commandstring)

		def avggene_eE():
			"""
			Calculate average gene for one read size
			"""
			# readsize = str(readin)
			# change this to read in fl reads
			pop = "eE"
			minlen = minlen_eE
			maxlen = maxlen_eE
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.samp, minlen, maxlen, self.samp, minlen, maxlen)
			outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, self.samp, pop, thresh, minlen, maxlen)
			print ""
			commandstring = 'python2 %s/riboseq/makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(
				rootDir,regionlength5,regionlength3,trspdictfilestring,self.UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)

			print commandstring
			os.system(commandstring)

		def avggene_aL():
			pop = "aL"
			# minlen = minlen_eE
			# maxlen = maxlen_eE
			trspdictfilestring = '%s/%s/density%sp_%sshift_%s/%s_%sf/%s_%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, pop, self.samp, pop, self.samp, pop)
			outfilebase = "%s/%s_rpkmThresh%s_%s_%sf" % (outfolder, self.samp, thresh, pop, pop)
			print ""
			commandstring = 'python2 %s/riboseq/makeavggene_main.py --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(
				rootDir,regionlength5,regionlength3,trspdictfilestring,self.UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)

			print commandstring
			os.system(commandstring)

		## writing multiprocessing fix: 
		p1 = multiprocessing.Process(target=avggene_fl)
		p1.start()
		p2 = multiprocessing.Process(target=avggene_eA)
		p2.start()
		p3 = multiprocessing.Process(target=avggene_eE)
		p3.start()
		p4 = multiprocessing.Process(target=avggene_aL)
		p4.start()
		p1.join()
		p2.join()
		p3.join()
		p4.join()


################################################################################################################

	### this needs to be adjusted to handle flexible read length sizes
	def avggene_riboshift_plot(self, alignposition, ribosome_site, fiveorthreeprime='5', normalizaiton='uneq', threshold='0'):
		"""
		plotting of avegene's for riboshifts
		"""
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # should be 5' mapped at this point
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

		thresh = threshold

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		# avggene_csv_path = '%s/avggene%s_ORFraw_%sshift_%s%s150' % (fp_assign_path, alignpos, assignpos, mappingpos, norm) 
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, self.norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now

		outpath = '%s/metaplots' % (fp_assign_path)
		if not os.path.exists(outpath):   os.makedirs(outpath)


		### get size ranges used for fl, eA, and eE
		fl_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eA_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eE_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)

		fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
		fl_ribo_shift.columns = ['riboshift']
		eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
		eA_ribo_shift.columns = ['riboshift']
		eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
		eE_ribo_shift.columns = ['riboshift']

		# print eE_ribo_shift.index.tolist()
		fl_sizes = fl_ribo_shift.index.tolist()
		eA_sizes = eA_ribo_shift.index.tolist()
		eE_sizes = eE_ribo_shift.index.tolist()

		minlen_fl = fl_sizes[0]
		maxlen_fl = fl_sizes[-1]
		minlen_eA = eA_sizes[0]
		maxlen_eA = eA_sizes[-1]
		minlen_eE = eE_sizes[0]
		maxlen_eE = eE_sizes[-1]
		###

		## get paths to stored csv average gene files
		fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s.csv' % (avggene_csv_path, self.samp, thresh, minlen_fl, maxlen_fl, alignpos)
		eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s.csv' % (avggene_csv_path, self.samp, thresh, minlen_eA, maxlen_eA, alignpos)
		eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s.csv' % (avggene_csv_path, self.samp, thresh, minlen_eE, maxlen_eE, alignpos)
		aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s.csv' % (avggene_csv_path, self.samp, thresh, alignpos)

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header=0)


		# fl plot
		fl_plot_outfile = "%s/%s_fl_%s_rpkmThresh%s_%s_%sshift_%s.pdf" % (outpath, self.samp, alignpos, thresh, self.norm_type, ribosome_shift, norm)
		fl_avggene_df.plot.line(y='avg', use_index = True)
		plt.savefig(fl_plot_outfile, format = 'pdf')
		plt.close()

		# eA plot
		eA_plot_outfile = "%s/%s_eA_%s_rpkmThresh%s_%s_%sshift_%s.pdf" % (outpath, self.samp, alignpos, thresh, self.norm_type, ribosome_shift, norm)
		eA_avggene_df.plot.line(y='avg', use_index = True)
		plt.savefig(eA_plot_outfile, format = 'pdf')
		plt.close()
		
		# eE plot 
		eE_plot_outfile = "%s/%s_eE_%s_rpkmThresh%s_%s_%sshift_%s.pdf" % (outpath, self.samp, alignpos, thresh, self.norm_type, ribosome_shift, norm)
		eE_avggene_df.plot.line(y='avg', use_index = True)
		plt.savefig(eE_plot_outfile, format = 'pdf')
		plt.close()

		# aL plot
		aL_plot_outfile = "%s/%s_aL_%s_rpkmThresh%s_%s_%sshift_%s.pdf" % (outpath, self.samp, alignpos, thresh, self.norm_type, ribosome_shift, norm)
		aL_avggene_df.plot.line(y='avg', use_index = True)
		plt.savefig(aL_plot_outfile, format = 'pdf')
		plt.close()


		# print fl_avggene_df
		


################################################################################################################

	#### Individual Metagene heatmaps

	def metagene_heatmap(self, alignposition, fiveorthreeprime, normalizaiton, ribosome_site='0'):
		"""
		Build a metagene_heatmap plot using R
		"""
		alignpos = alignposition	# '1' for start '2' for stop
		assignment = fiveorthreeprime  	# '5' for fiveprime '3' for threeprime
		norm = normalizaiton 	# 'uneq' or 'eq'
		ribosome_shift = ribosome_site # A P or E

		if alignpos == '1':
			pos = 'START'
		if alignpos == '2':
			pos = 'STOP' # numbering on STOP plot needs some help at the moment

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		
		plot_outpath = '%s/metaplots' % (fp_assign_path)
		if not os.path.exists(plot_outpath):   os.makedirs(plot_outpath)

		wd = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # the last 150 value for cdslength is hardcoded at the moment
		#print wd
		Rcode1 = '%s/riboseq/%s_heatmap.R' % (rootDir,pos)
		Rcode2 = '%s/riboseq/%s_heatmap_local.R' % (rootDir,pos)
		outfile1 = '%s/%s_%s_%sp_%sshift_heatmap_%s_reads.pdf' % (plot_outpath, self.samp, pos, assignment, ribosome_shift, self.norm_type)
		outfile2 = '%s/%s_%s_%sp_%sshift_heatmap_%s_reads_local.pdf' % (plot_outpath, self.samp, pos, assignment, ribosome_shift, self.norm_type)

		heatmap_wide = 'Rscript %s %s %s %s' % (Rcode1, self.samp, wd, outfile1)
		heatmap_local = 'Rscript %s %s %s %s' % (Rcode2, self.samp, wd, outfile2)
		subprocess.Popen(heatmap_wide, shell= True).wait()
		subprocess.Popen(heatmap_local, shell= True).wait()	



	def get_riboshift(self, shiftsite, alignposition="1", fiveorthreeprime="5", normalizaiton="uneq", ribosome_site="0"):
		"""
		Function that calculates A offsets for a library
		    Offsets are from 5' aligned reads
		    Mapping is relative to start codon
		    Dependent on running avggene_multi(alignposition="1", fiveorthreeprime="5", normalizaiton="uneq",ribosome_site="0")
		"""
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

		shiftpos = shiftsite
		if shiftpos == "P":
			shift_adjust = 0
		if shiftpos == "A":
			shift_adjust = 3 # A site is shifted 3 nt's farther than P site from 5' end
		if shiftpos == "E":
			shift_adjust = -3


		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		avggene_paths = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, self.norm_type, ribosome_shift, assignment, norm) # must use 5' start mapping for this

		frame = pd.DataFrame()
		list_ = []

		for ft in self.ftsize:
			infile = "%s/%s_%sf_avg_1.csv" % (avggene_paths, self.samp, ft)
			df = pd.read_csv(infile, index_col="position", header=0)
			df.columns = ["%d" % (ft)]
			list_.append(df)

		frame = pd.concat(list_, axis=1)

		# Get the max (or modal read number per position) of footprints relative to start codon
		# Adding a +1 to account of +/- 1nt from digestion heterogeneity
		# Only look between 11 to 16 nts upstream of start for 5' aligned reads
		# Define full length reads from 31-34 nt's in length
		    # retrieve these using 31-15 == 16, and 34-15 == 19+1, thus 16:20
		offset_fl = 50+1+shift_adjust-frame.iloc[34:39,self.flmin-15:self.flmax+1-15].idxmax()
		# Define emptyA reads from 21-24 nt's in length
		    # retrieve these using 21-15 == 6, and 24-15 == 9+1, thus 6:10
		offset_emptyA = 50+1+shift_adjust-frame.iloc[34:39,self.eAmin-15:self.eAmax+1-15].idxmax()
		# Define emptyE reads from 18-19 nt's in length - need to cut close to avoid empty A
			# retrieve these using 18-15 == 3, and 19-15 == 4, thus 3:5
			# also need to only look around the start codon for these reads since the 5' end is directly over the 'A' of the AUG
		offset_emptyE = 50+1+shift_adjust-frame.iloc[49:52,self.eEmin-15:self.eEmax+1-15].idxmax()

		outpath = "%s/FPassignment/%s/%s/%s/riboshift" % (self.rootpath, self.genome_name, self.experiment, self.samp)
		if not os.path.exists(outpath):     os.makedirs(outpath)

		outfile_fl = "%s/full_length_%soffsets.csv" % (outpath, shiftpos)
		outfile_emptyA = "%s/emptyA_%soffsets.csv" % (outpath, shiftpos)
		outfile_emptyE = "%s/emptyE_%soffsets.csv" % (outpath, shiftpos)
		
		offset_fl.to_csv(outfile_fl, index = True, header = False)
		offset_emptyA.to_csv(outfile_emptyA, index = True, header = False)
		offset_emptyE.to_csv(outfile_emptyE, index = True, header = False)
		

	def plot_riboshift_metagenes(self, shiftsite, alignposition="1", fiveorthreeprime="5", normalizaiton="uneq", ribosome_site="0"):
		"""
		Plot average gene plots using A site shifts from get_riboshift_Asites()
		"""
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

		shiftpos = shiftsite # should be "P" or "A"

		# I believe I could also load ALL_start_pops.csv file here intead of regenerating the dataFrame with some tweaking
		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		avggene_paths = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, self.norm_type, ribosome_shift, assignment, norm) # must use 5' start mapping for this

		frame = pd.DataFrame()
		list_ = []

		for ft in self.ftsize:
			infile = "%s/%s_%sf_avg_1.csv" % (avggene_paths, self.samp, ft)
			df = pd.read_csv(infile, index_col="position", header=0)
			df.columns = ["%d" % (ft)]
			list_.append(df)

		frame = pd.concat(list_, axis=1)

		frame_fl = frame.iloc[:,self.flmin-15:self.flmax+1-15] # account for large mounts of reads
		frame_emptyA = frame.iloc[:,self.eAmin-15:self.eAmax+1-15]
		frame_emptyE = frame.iloc[:,self.eEmin-15:self.eEmax+1-15]

		frame_fl['readpos'] = frame_fl.index
		frame_emptyA['readpos'] = frame_emptyA.index
		frame_emptyE['readpos'] = frame_emptyE.index   

		shift_path = "%s/FPassignment/%s/%s/%s/riboshift" % (self.rootpath, self.genome_name, self.experiment, self.samp)
		shift_fl = "%s/full_length_%soffsets.csv" % (shift_path, shiftpos)
		shift_eA = "%s/emptyA_%soffsets.csv" % (shift_path, shiftpos)
		shift_eE = "%s/emptyE_%soffsets.csv" % (shift_path, shiftpos)

		fl_shift_df = pd.read_csv(shift_fl, header=None, index_col=0)
		eA_shift_df = pd.read_csv(shift_eA, header=None, index_col=0) 
		eE_shift_df = pd.read_csv(shift_eE, header=None, index_col=0) 	    

		fl_dict = fl_shift_df.iloc[:,0].to_dict()
		eA_dict = eA_shift_df.iloc[:,0].to_dict()
		eE_dict = eE_shift_df.iloc[:,0].to_dict()

		df_fl_list = []
		df_eA_list = []
		df_eE_list = []

		for i in fl_shift_df.index: 
			df_name = 'fl_df_'+str(i)
			fl_df = frame_fl.loc[:,[str(i), 'readpos']]
			fl_df['readpos'] = fl_df['readpos']+(fl_dict[i])
			fl_df.name = df_name
			df_fl_list.append(fl_df)

		for i in eA_shift_df.index: 
			df_name = 'eA_df_'+str(i)
			eA_df = frame_emptyA.loc[:,[str(i), 'readpos']]
			eA_df['readpos'] = eA_df['readpos']+(eA_dict[i])
			eA_df.name = df_name
			df_eA_list.append(eA_df)

		for i in eE_shift_df.index: 
			df_name = 'eE_df_'+str(i)
			eE_df = frame_emptyE.loc[:,[str(i), 'readpos']]
			eE_df['readpos'] = eE_df['readpos']+(eE_dict[i])
			eE_df.name = df_name
			df_eE_list.append(eE_df)

		## eventually this could be changed to add after every call to subplot
			# I'm not sure how to draw lines on the individual plots at this point
		def plot_features():
			plt.axvline(x=50-0.5,color="#999999",dashes=[3,2],zorder=-1)
			plt.axvline(x=52+0.5,color="#999999",dashes=[3,2],zorder=-1)
			plt.axvline(x=55+0.5,color="#999999",dashes=[3,2],zorder=-1)
			plt.text(x=50+1, y=0, s="P")
			plt.text(x=53+1, y=0, s="A")

		### Plot full length reads:
		# need some for loop here in the futrue that is more flexible for variable length reads
		fl_plot_outfile = "%s/fullLength_%s_%sshift.pdf" % (shift_path, self.norm_type, shiftpos)
		fig, axes = plt.subplots(nrows=len(df_fl_list), ncols=1, sharex=True)
		for i in range(len(df_fl_list)):
		    df_fl_list[i].plot.line(x='readpos',y=df_fl_list[i].columns.values[0], ax=axes[i])
		plt.xlim(xmin = 35, xmax = 70)
		plt.axvline(x=50-0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.axvline(x=52+0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.axvline(x=55+0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.text(x=50+0.5, y=0, s="P")
		plt.text(x=53+0.5, y=0, s="A")
		plt.savefig(fl_plot_outfile, format='pdf')
		plt.close()

		### Now plot emptyA reads:
		eA_plot_outfile = "%s/emptyA_%s_%sshift.pdf" % (shift_path, self.norm_type, shiftpos)
		fig, axes = plt.subplots(nrows=len(df_eA_list), ncols=1, sharex=True)
		for i in range(len(df_eA_list)):
		    df_eA_list[i].plot.line(x='readpos',y=df_eA_list[i].columns.values[0], ax=axes[i])
		plt.xlim(xmin = 35, xmax = 70)
		plt.axvline(x=50-0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.axvline(x=52+0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.axvline(x=55+0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.text(x=50+0.5, y=0, s="P")
		plt.text(x=53+0.5, y=0, s="A")
		plt.savefig(eA_plot_outfile, format='pdf')
		plt.close()

		### Now plot emptyE reads:
		eE_plot_outfile = "%s/emptyE_%s_%sshift.pdf" % (shift_path, self.norm_type, shiftpos)
		fig, axes = plt.subplots(nrows=len(df_eE_list), ncols=1, sharex=True)
		for i in range(len(df_eE_list)):
		    df_eE_list[i].plot.line(x='readpos',y=df_eE_list[i].columns.values[0], ax=axes[i])
		df_eE_list[0].plot.line(x='readpos', y='18', ax=axes[0])
		df_eE_list[1].plot.line(x='readpos', y='19', ax=axes[1])
		# df_eE_list[2].plot.line(x='readpos', y='23', ax=axes[2])
		# df_eE_list[3].plot.line(x='readpos', y='24', ax=axes[3])
		plt.axvline(x=50-0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.axvline(x=52+0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.axvline(x=55+0.5,color="#999999",dashes=[3,2],zorder=-1)
		plt.text(x=50+0.5, y=0, s="P")
		plt.text(x=53+0.5, y=0, s="A")
		plt.xlim(xmin = 35, xmax = 70)
		plt.savefig(eE_plot_outfile, format='pdf')
		plt.close()


	### building count tables:
	def build_count_tables(self, alignposition="1", fiveorthreeprime="5", normalizaiton="uneq", ribosome_site="A"):
		"""
		wrapper to direct building of density tables
		count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
		inset_choice contains adjustments to avoid start and stop codon peaks in data tables
		output is a csv file with densities for each region and 3'UTR occupancies
		"""
		# mappingpos = "5"
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

		flminlen = str(self.flmin)
		flmaxlen = str(self.flmax)
		eAminlen = str(self.eAmin)
		eAmaxlen = str(self.eAmax)
		eEminlen = str(self.eEmin)
		eEmaxlen = str(self.eEmax)


		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store count tables
		if not os.path.exists(outfolder):   os.makedirs(outfolder)
		
		# if self.norm_type == "raw":
		# 	densitystring = "DensityUnnormalized"
		# elif self.norm_type == "rpm":
		# 	densitystring = "Density_rpm"
		# else:
		# 	print "Normalization is not set for densebuilder!"
		# 	sys.exit()

		### require Density_rpm for densityfilestring
		densitystring = "Density_rpm"
		totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.samp)
		totreadcountf = open(totreads_countfile, "r")
		totreads = int(totreadcountf.read())
		totreadcountf.close()


		def fl_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, flminlen, flmaxlen, self.samp, flminlen, flmaxlen, self.samp, flminlen, flmaxlen)
			outfilestring = "%s/%s_fl_%s_%sto%s_countTable_rpkm" % (outfolder, self.samp, self.norm_type, flminlen, flmaxlen)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads)

			print commandstring
			os.system(commandstring)

		def eA_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, eAminlen, eAmaxlen, self.samp, eAminlen, eAmaxlen, self.samp, eAminlen, eAmaxlen)
			outfilestring = "%s/%s_eA_%s_%sto%s_countTable_rpkm" % (outfolder, self.samp, self.norm_type, eAminlen, eAmaxlen)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads)
			print commandstring
			os.system(commandstring)

		def eE_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, eEminlen, eEmaxlen, self.samp, eEminlen, eEmaxlen, self.samp, eEminlen, eEmaxlen)
			outfilestring = "%s/%s_eE_%s_%sto%s_countTable_rpkm" % (outfolder, self.samp, self.norm_type, eEminlen, eEmaxlen)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads)
			print commandstring
			os.system(commandstring)

		def aL_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_aL/%s_aLf/%s_aLf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, self.samp, self.samp,)
			outfilestring = "%s/%s_aL_%s_countTable_rpkm" % (outfolder, self.samp, self.norm_type)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads)
			print commandstring
			os.system(commandstring)


		p1 = multiprocessing.Process(target=fl_countTable)
		p1.start()
		p2 = multiprocessing.Process(target=eA_countTable)
		p2.start()
		p3 = multiprocessing.Process(target=eE_countTable)
		p3.start()
		p4 = multiprocessing.Process(target=aL_countTable)
		p4.start()
		p1.join()
		p2.join()
		p3.join()
		p4.join()
		# def eA_countTable():

		# def eE_countTable():

		# def aL_countTable():


	def build_count_tables_utr3adj(self, alignposition="1", fiveorthreeprime="5", normalizaiton="uneq", ribosome_site="A"):
		"""
		wrapper to direct building of density tables
		count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
		inset_choice contains adjustments to avoid start and stop codon peaks in data tables
		output is a csv file with densities for each region and 3'UTR occupancies
		"""
		# mappingpos = "5"
		stopcodons = self.stopcodons
		alignpos = alignposition # '1' for start, '2' for stop
		assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
		norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

		flminlen = str(self.flmin)
		flmaxlen = str(self.flmax)
		eAminlen = str(self.eAmin)
		eAmaxlen = str(self.eAmax)
		eEminlen = str(self.eEmin)
		eEmaxlen = str(self.eEmax)


		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store count tables
		if not os.path.exists(outfolder):   os.makedirs(outfolder)
		
		# if self.norm_type == "raw":
		# 	densitystring = "DensityUnnormalized"
		# elif self.norm_type == "rpm":
		# 	densitystring = "Density_rpm"
		# else:
		# 	print "Normalization is not set for densebuilder!"
		# 	sys.exit()

		### require Density_rpm for densityfilestring
		densitystring = "Density_rpm"
		totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.samp)
		totreadcountf = open(totreads_countfile, "r")
		totreads = int(totreadcountf.read())
		totreadcountf.close()


		def fl_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, flminlen, flmaxlen, self.samp, flminlen, flmaxlen, self.samp, flminlen, flmaxlen)
			outfilestring = "%s/%s_fl_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, self.samp, self.norm_type, flminlen, flmaxlen)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm_utr3adj.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads, stopcodons)

			print commandstring
			os.system(commandstring)

		def eA_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, eAminlen, eAmaxlen, self.samp, eAminlen, eAmaxlen, self.samp, eAminlen, eAmaxlen)
			outfilestring = "%s/%s_eA_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, self.samp, self.norm_type, eAminlen, eAmaxlen)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm_utr3adj.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads, stopcodons)
			print commandstring
			os.system(commandstring)

		def eE_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, eEminlen, eEmaxlen, self.samp, eEminlen, eEmaxlen, self.samp, eEminlen, eEmaxlen)
			outfilestring = "%s/%s_eE_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, self.samp, self.norm_type, eEminlen, eEmaxlen)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm_utr3adj.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads, stopcodons)
			print commandstring
			os.system(commandstring)

		def aL_countTable():
			trspdictfilestring = '%s/%s/density%sp_%sshift_aL/%s_aLf/%s_aLf_' %(
					fp_assign_path, densitystring, assignment, ribosome_shift, self.samp, self.samp,)
			outfilestring = "%s/%s_aL_%s_countTable_rpkm_utr3adj" % (outfolder, self.samp, self.norm_type)
			print ""
			# need to input trspdictfilestring, UTRfilesting, cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice

			commandstring = 'python2 %s/riboseq/riboseq_buildDenseTables_rpkm_utr3adj.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (
				rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring, totreads, stopcodons)
			print commandstring
			os.system(commandstring)


		p1 = multiprocessing.Process(target=fl_countTable)
		p1.start()
		p2 = multiprocessing.Process(target=eA_countTable)
		p2.start()
		p3 = multiprocessing.Process(target=eE_countTable)
		p3.start()
		p4 = multiprocessing.Process(target=aL_countTable)
		p4.start()
		p1.join()
		p2.join()
		p3.join()
		p4.join()

	# def 



#####################################
#		END of RPexperiment 		#
#####################################



### metagene heatmap plotting:
class seaborn_metagene_heatmaps(object):
	'''
	use seaborn to generate automated plot of all metagene heatmaps in the dataset
	create metagene heatmaps around the start and stop codon
		Use both RPM values scalled within each dataset and robust settings to allow visualization of majority of data set
			Robust just reduces vmax and increases vmin, helps with high abundance start peaks
	'''
	def __init__(self, rootpath, genome_name, experiment, samples, ftsize):
		self.rootpath = rootpath
		self.genome_name = genome_name
		self.experiment = experiment
		self.samples = samples
		self.ftsize = ftsize

	
	def build_tidy_df(self, position, folder, file):
		"""
		build a tidy dataframe containing all samples in a single frame
		position must be 'start' or 'stop'
		takes each average gene file for each read length, creates a data frame, 
			and then combines these all into a single df for each sample
		"""
		dframe = pd.DataFrame()
		df_list = []
		avgpath = "%s/FPassignment/%s/%s/%s" % (self.rootpath, self.genome_name, folder, file)
		if position == "start":
			pospath = 'avggene1_ORFrpm_0shift_5uneq150'
			avggennumb = '1'
		if position == "stop":
			pospath = 'avggene2_ORFrpm_0shift_5uneq150'
			avggennumb = '2'

		for i in ftsize[::-1]: # reverse order of ftsize to start with 40 and work backward to 15, possibly simplifies plotting
			readlen = str(i)
			df = pd.read_csv("%s/%s/%s_%sf_avg_%s.csv" % (avgpath, pospath, file, readlen, avggennumb))
			df['sample'] = file # create a new column with the sample name
			df['readlen'] = readlen # create another new column for the readlength being measured
			df_list.append(df) # add this to the empty placeholder list
		dframe = pd.concat(df_list, axis=0) # combine all dataframes into a single dataframe, combining by the rows - keeping same columns
		return dframe

	def concat_tidy_dfs(self, position):
		"""
		combine all tidy dataframes for each sample, into one master dataframe for all samples
		also, this is where local values are set
		"""
		for folder in self.samples:
			all_samples = [] # an empty list that will eventually contain all entries for all samples
			for file in self.samples[folder]:
				out_df = self.build_tidy_df(position, folder, file) # build a tidy data frame for each sample
				all_samples.append(out_df)
			df_all = pd.concat(all_samples, axis=0, ignore_index=True)
			if position == 'start':
				df_all['position'] = df_all['position']-50
			if position == 'stop':
				df_all['position'] = df_all['position']-150
		df_local = df_all.copy()
		if position == 'start':
			local_min = -20
			local_max = 20
		if position == 'stop':
			local_min = -30
			local_max = 20
		df_local = df_local.loc[(df_local['position']>= local_min) & (df_local['position']<=local_max)]
		return df_all, df_local

	def draw_heatmap(self, *args, **kwargs):
		"""
		Define the heatmap function to be used for plotting the data
		All of the pyplot variables are set here
		"""
		data = kwargs.pop('data')
		### I'm not entirely sure why we need to use pivot... something to do with how fg.map_dataframe feeds in data to this function
		d = data.pivot(index=args[1], columns=args[0], values=args[2]) 
		ax = sns.heatmap(d, cmap="binary", **kwargs)
		ax.invert_yaxis()
		
		### set x ticks
		xtix = range(0,len(d.columns))
		xtix = [x+0.5 for x in xtix]
		xtixminor = xtix[::2]
		xtixmajor = xtix[::10]

		xtixlabels = range(data['position'].min(), data['position'].max()+1)
		xtixlabels = xtixlabels[::10]

		ax.set_xticklabels(xtixlabels, horizontalalignment='center')
		ax.set_xticks(xtixminor, minor=True)
		ax.set_xticks(xtixmajor, minor=False)

		### set y ticks
		ytix = range(0,25)
		ytix = [y+0.5 for y in ytix]
		ytixminor = ytix[::1]
		ytixmajor = ytix[::5]

		ytixlabels = range(15,41)
		ytixlabels = ytixlabels[::5]

		ax.set_yticks(ytixminor, minor=True)
		ax.set_yticks(ytixmajor, minor=False)
		ax.set_yticklabels(ytixlabels, verticalalignment='center')

	def meta_heatmap(self, position, df, meta_window):
		"""
		use a facetgrid to output a heatmap from each library to the same figure
		map the 'draw_heatmap' function defined above
		"""
		plotoutdir = "%s/FPassignment/%s/%s/analysis/heatmaps" % (self.rootpath, self.genome_name, self.experiment)
		if not os.path.exists(plotoutdir):	os.makedirs(plotoutdir)
		if meta_window == "all":
			asp = 4
		if meta_window == "local":
			asp = 3
		fg = sns.FacetGrid(df, row="sample", aspect=asp) # aspect sets the aspect ratio to be what we want
		fg.map_dataframe(self.draw_heatmap, 'position', 'readlen', 'avg', robust=False)
		if position == 'start':
			metapos = 'start'
		if position == 'stop':
			metapos = 'stop'

		plotoutfig_pdf = "%s/meta_heatmap_%s_%s_%s.pdf" % (plotoutdir, self.experiment, metapos, meta_window)
		plt.savefig(plotoutfig_pdf)

	def meta_heatmap_robust(self, position, df, meta_window):
		plotoutdir = "%s/FPassignment/%s/%s/analysis/heatmaps" % (self.rootpath, self.genome_name, self.experiment)
		if not os.path.exists(plotoutdir):	os.makedirs(plotoutdir)
		fg = sns.FacetGrid(df, row="sample", aspect=4) # aspect sets the aspect ratio to be what we want
		fg.map_dataframe(self.draw_heatmap, 'position', 'readlen', 'avg', robust=True)
		if position == 'start':
			metapos = 'start'
		if position == 'stop':
			metapos = 'stop'
		plotoutfig_pdf = "%s/meta_heatmap_robust_%s_%s_%s.pdf" % (plotoutdir, self.experiment, metapos, meta_window)
		plt.savefig(plotoutfig_pdf)

	def metaplotter(self):
		"""
		handler function that runs runs all of the above functions in the proper order
		set to build start and stop heatmaps by default
		"""
		poslist = ['start', 'stop']
		for pos in poslist:
			df_all, df_local = self.concat_tidy_dfs(pos)
			self.meta_heatmap(pos, df_all, meta_window="all")
			self.meta_heatmap(pos, df_local, meta_window="local")
			self.meta_heatmap_robust(pos, df_all, meta_window="all")
			self.meta_heatmap_robust(pos, df_local, meta_window="local")

	######


def main():
	for experiment in samples:
		filenumb = 0 # add a counter to keep track of which sample is being analyzed
		for samp in samples[experiment]:
			print "******  starting sample: ", samples[experiment][filenumb], "******"
			
			### adjust for variable linkers used in library prep

			trim5Pval = trim5Plist[filenumb]

			if trim5Pval == 0:
				linker = linker1
				fivePrimeNTstoTrim = "0"
				deduplicate_reads = False
				print "trimming 5prime: ", fivePrimeNTstoTrim, "and using linker: ", linker, "dedupe == ", deduplicate_reads
			elif trim5Pval == 1:
				linker = linker2
				fivePrimeNTstoTrim = "4"
				deduplicate_reads = True
				print "trimming 5prime: ", fivePrimeNTstoTrim, "and using linker: ", linker, "dedupe == ", deduplicate_reads

			RPexp = RPexperiment(samp, experiment, rootpath, fastqPath, fivePrimeNTstoTrim, threadNumb, 
								genome_name, ncRNAstarGenome, starGenome, linker, GTFfile, 
								twobitfile, threshold, UTRfilestring, ftsize, flmin, flmax, eAmin, eAmax, eEmin, eEmax,
								cdsDenThresh, raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice, deduplicate_reads,
								read_proc_dir, align_dir, stopcodons)
			
			## Read Processing
			RPexp.deduplicate_input_reads()
			RPexp.cut_5p_end()
			RPexp.remove_adapter()
			RPexp.ncRNA_subtract_STAR()
			# RPexp.STAR_align() ### unique alignments
			RPexp.STAR_align_Single_Random_Mulitmapper() ### take single best multimapper


			### Density builders
			RPexp.densebuilder_multi_raw("5")
			RPexp.densebuilder_multi_norm("5")		
			# RPexp.densebuilder_multi("3")

			### Avggene builders
			RPexp.avggene_multi(alignposition="1", fiveorthreeprime="5", normalizaiton="uneq", ribosome_site="0") 
			RPexp.avggene_multi(alignposition="2", fiveorthreeprime="5", normalizaiton="uneq", ribosome_site="0")
			# RPexp.avggene_multi(alignposition="1", fiveorthreeprime="3", normalizaiton="uneq", ribosome_site="0")
			# RPexp.avggene_multi(alignposition="2", fiveorthreeprime="3", normalizaiton="uneq", ribosome_site="0")
			
			### individual metaplots			
			RPexp.metagene_heatmap("1", "5", "uneq") 
			RPexp.metagene_heatmap("2", "5", "uneq")
			# RPexp.metagene_heatmap("1", "3", "uneq")
			# RPexp.metagene_heatmap("2", "3", "uneq")

			### read size distribution
			RPexp.read_size_dist_single()
			
			### QC
			# RPexp.read_count_tracking()
			# RPexp.ncRNA_distribution()

			### Riboshift offsets:
			RPexp.get_riboshift("P")
			RPexp.plot_riboshift_metagenes("P")
			RPexp.get_riboshift("A")	
			RPexp.plot_riboshift_metagenes("A") 

			### density building with riboshifted reads:
			RPexp.densebuilder_riboshift("A")
			RPexp.avggene_riboshift("1","A")
			RPexp.avggene_riboshift("2","A")
			RPexp.avggene_riboshift_plot("1", "A")
			RPexp.avggene_riboshift_plot("2", "A")

			RPexp.avggene_riboshift("1","A", normalizaiton="eq", threshold='10')
			RPexp.avggene_riboshift("2","A", normalizaiton="eq", threshold='10')
			RPexp.avggene_riboshift_plot("1", "A", normalizaiton="eq", threshold='10')
			RPexp.avggene_riboshift_plot("2", "A", normalizaiton="eq", threshold='10')

			### build count tables:
			RPexp.build_count_tables()
			RPexp.build_count_tables_utr3adj()

			filenumb+=1 

	### Seaborn Metagene Heatmaps
	snsHeatmap = seaborn_metagene_heatmaps(rootpath, genome_name, experiment, samples, ftsize)
	snsHeatmap.metaplotter()
	

if __name__ == '__main__':
	# execute only if run as a script
	main()