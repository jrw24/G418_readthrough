__author__ = "jamie wangen"

"""
main file for handling RNAseq reads here
"""

## start with SR1 datasets:


### Import dependencies
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import argparse
import subprocess
import multiprocessing
import collections
import pysam
import pandas as pd 
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
import importlib


### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/RNAseq" % args.rootDir)
sys.path.append("%s/RNAseq/libsettings" % args.rootDir)

rootDir = args.rootDir
libsetName = args.libSetFile

print libsetName
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)

threadNumb = str(args.threadNumb)


class RNAseqExperiment(object):
	"""
	Workflow: trim adapter, ncRNA subtraction, STAR align -> check the alignments

	"""
	def __init__(self, samp, experiment, rootpath, fastqPath, threadNumb, genome_name, 
					ncRNAstarGenome, starGenome, linker, GTFfile, twobitfile,
					UTRfilestring, raw_dense_thresh, rpm_dense_thresh, norm_type, 
					inset_choice, paired_end_reads, threshold, read_proc_dir, align_dir):

		self.samp = samp 
		self.experiment = experiment
		self.rootpath = rootpath
		self.fastqPath = fastqPath
		self.threadNumb = threadNumb
		self.genome_name = genome_name
		self.ncRNAstarGenome = ncRNAstarGenome
		self.starGenome = starGenome
		self.linker = linker
		self.GTFfile = GTFfile
		self.twobitfile = twobitfile
		self.UTRfilestring = UTRfilestring
		self.raw_dense_thresh = raw_dense_thresh
		self.rpm_dense_thresh = rpm_dense_thresh
		self.norm_type = norm_type
		self.inset_choice = inset_choice
		self.paired_end_reads = paired_end_reads
		self.threshold = threshold
		self.cdsDenThresh = cdsDenThresh
		self.read_proc_dir = read_proc_dir
		self.align_dir = align_dir



	"""


	"""

	def remove_adapter(self):
		'''
		-x specifies the 3' adapter on on the read
		-l specifies the minimum post-trimming read length
		-L specifies the maximum post-trimming read length
		-Q specifies the minimum quality score 

		For TruSeq Kits:
		Illumina Read 1 Adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
		Illumina Read 2 Adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
		'''
		# input fastq files:
		fastqPath = self.fastqPath
		
		if self.paired_end_reads == True:
			fastq1 = '%s/%s*1.fastq.gz' % (fastqPath, self.samp)
			fastq2 = '%s/%s*2.fastq.gz' % (fastqPath, self.samp)

			trimming_out_path = '%s/%s' % (self.rootpath, self.read_proc_dir)
			if not os.path.exists(trimming_out_path): os.makedirs(trimming_out_path)
			trimming_out_file = '%s/%s' % (trimming_out_path, self.samp)

			trimming_out_filegz_1 = '%s/%s-trimmed-pair1.fastq.gz' % (trimming_out_path, self.samp)
			trimming_out_filegz_2 = '%s/%s-trimmed-pair2.fastq.gz' % (trimming_out_path, self.samp)
			
			if os.path.isfile(trimming_out_filegz_1) and os.path.isfile(trimming_out_filegz_2):
				print "adapter previously removed"
			else: 
				# leaving this without a defined adapter for now...
				# skewer is using -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC, -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
				read_trimming = 'skewer -m pe -Q 10 -l 15 --quiet -o %s --threads %s %s %s' % (
								trimming_out_file, self.threadNumb, fastq1, fastq2)
				print read_trimming
				subprocess.Popen(read_trimming, shell=True).wait()
				
				compression_command_1 = 'pigz -p %s %s-trimmed-pair1.fastq' % (self.threadNumb, trimming_out_file)
				compression_command_2 = 'pigz -p %s %s-trimmed-pair2.fastq' % (self.threadNumb, trimming_out_file)
				print compression_command_1
				subprocess.Popen(compression_command_1, shell=True).wait()
				print compression_command_2
				subprocess.Popen(compression_command_2, shell=True).wait()
				print "adapter trimming finished"

		# not yet working for unpaired reads...
		if self.paired_end_reads == False:

			fastq_name = '%s/%s*fastq.gz' % (fastqPath, self.samp)

			trimming_out_path = '%s/%s' % (self.rootpath, self.read_proc_dir)
			if not os.path.exists(trimming_out_path): os.makedirs(trimming_out_path)
			
			trimming_out_file = '%s/%s' % (trimming_out_path, self.samp)
			trimming_out_filegz = '%s-trimmed.fastq.gz' % (trimming_out_file)
			
			if os.path.isfile(trimming_out_filegz):
				print "adapter previously removed"
			else:
				read_trimming = 'skewer -x %s -Q 10 -l 15 --quiet -o %s --threads %s %s' % (
					self.linker, trimming_out_file, self.threadNumb, fastq_name)
				print read_trimming
				subprocess.Popen(read_trimming, shell=True).wait()
				
				compression_command = 'pigz -p %s %s-trimmed.fastq' % (self.threadNumb, trimming_out_file)
				print compression_command
				subprocess.Popen(compression_command, shell=True).wait()




	def ncRNA_subtract_STAR(self):
		'''
		using --alignIntronMax 1 to prevent splicing
		set to 1 mismatch using --outFilterMismatchNmax 1
		'''
		ncRNA_subSTAR_outpath = '%s/%s' % (self.rootpath, self.read_proc_dir)
		ncRNA_subSTAR_outfile = '%s/%s' % (ncRNA_subSTAR_outpath, self.samp)
		if not os.path.exists(ncRNA_subSTAR_outpath): os.makedirs(ncRNA_subSTAR_outpath)

		if self.paired_end_reads == True:
			ncRNA_sub_star = 'STAR \
							  --runThreadN %s \
							  --genomeDir %s \
							  --readFilesIn %s-trimmed-pair1.fastq.gz %s-trimmed-pair2.fastq.gz \
							  --readFilesCommand gunzip -c \
							  --outReadsUnmapped Fastx \
							  --outSAMtype None \
							  --outFilterMismatchNmax 1 \
							  --outFilterMultimapNmax 100 \
							  --alignIntronMax 1 \
							  --outFileNamePrefix %s' % (
				self.threadNumb, self.ncRNAstarGenome, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
			print ncRNA_sub_star
			subprocess.Popen(ncRNA_sub_star, shell=True).wait()
			compression_command_1 = 'pigz -p %s -c %sUnmapped.out.mate1 > %s_no_ncRNA_1.fastq.gz' % (self.threadNumb, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
			compression_command_2 = 'pigz -p %s -c %sUnmapped.out.mate2 > %s_no_ncRNA_2.fastq.gz' % (self.threadNumb, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
			print compression_command_1
			subprocess.Popen(compression_command_1, shell=True).wait()
			subprocess.Popen('rm %sUnmapped.out.mate1' % (ncRNA_subSTAR_outfile), shell=True).wait() 
			print compression_command_2
			subprocess.Popen(compression_command_2, shell=True).wait()
			subprocess.Popen('rm %sUnmapped.out.mate2' % (ncRNA_subSTAR_outfile), shell=True).wait()
			print "Non-coding RNA successfully removed \n" 

		# not yet validated for unpaired libraries ...
		if self.paired_end_reads == False:
			ncRNA_sub_star = 'STAR --runThreadN %s --genomeDir %s --readFilesIn %s-trimmed.fastq.gz --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outSAMtype None --outFilterMismatchNmax 1 --alignIntronMax 1 --outFileNamePrefix %s' % (
				self.threadNumb, self.ncRNAstarGenome, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
			print ncRNA_sub_star
			subprocess.Popen(ncRNA_sub_star, shell=True).wait()
			compression_command = 'pigz -p %s -c %sUnmapped.out.mate1 > %s_no_ncRNA.fastq.gz' % (self.threadNumb, ncRNA_subSTAR_outfile, ncRNA_subSTAR_outfile)
			subprocess.Popen(compression_command, shell=True).wait()
			subprocess.Popen('rm %sUnmapped.out.mate1' % (ncRNA_subSTAR_outfile), shell=True).wait() 


	def STAR_align(self):
		"""
		Align reads to the genome using STAR
		Make sure to allow for soft-clipping to permit alignment of reads with A or T added by superscript III
		45 seems to be the max number of supported threads for alignment
		"""
		if self.paired_end_reads == True:
			star_infile_1 = '%s/%s/%s_no_ncRNA_1.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
			star_infile_2 = '%s/%s/%s_no_ncRNA_2.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
			star_out_path = '%s/%s/%s/STAR/%s' % (self.rootpath, self.align_dir, self.genome_name, self.samp)
			star_out_file = '%s/%s' % (star_out_path, self.samp)
			if not os.path.exists(star_out_path): os.makedirs(star_out_path)
			
			star_command = 'STAR \
							--runThreadN %s \
							--genomeDir %s \
							--readFilesIn %s %s \
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
				self.threadNumb, self.starGenome, star_infile_1, star_infile_2, star_out_file)
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
			print "STAR alignment has finished successfully \n"

		
		# not yet validated...
		if self.paired_end_reads == False:
			star_infile = '%s/%s/%s_no_ncRNA.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
			star_out_path = '%s/%s/%s/STAR/%s' % (self.rootpath, self.align_dir, self.genome_name, self.samp)
			star_out_file = '%s/%s' % (star_out_path, self.samp)
			if not os.path.exists(star_out_path): os.makedirs(star_out_path)
			
			star_command = 'STAR --runThreadN %s --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterType BySJout --outFilterMultimapNmax 1 --outWigType wiggle read1_5p --outFileNamePrefix %s --quantMode TranscriptomeSAM --outReadsUnmapped Fastx' % (self.threadNumb, self.starGenome, star_infile, star_out_file)
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


	def STAR_align_rMult(self):
		"""
		Align reads to the genome using STAR, allow multimapping and take the best match
		Make sure to allow for soft-clipping to permit alignment of reads with A or T added by superscript III
		45 seems to be the max number of supported threads for alignment
		"""
		if self.paired_end_reads == True:
			star_infile_1 = '%s/%s/%s_no_ncRNA_1.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
			star_infile_2 = '%s/%s/%s_no_ncRNA_2.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
			star_out_path = '%s/%s/%s/STAR/%s' % (self.rootpath, self.align_dir, self.genome_name, self.samp)
			star_out_file = '%s/%s' % (star_out_path, self.samp)
			if not os.path.exists(star_out_path): os.makedirs(star_out_path)
			
			star_command = 'STAR \
							--runThreadN %s \
							--genomeDir %s \
							--readFilesIn %s %s \
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
				self.threadNumb, self.starGenome, star_infile_1, star_infile_2, star_out_file)
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
			print "STAR alignment has finished successfully \n"

		
		# # not currently validated...
		# if self.paired_end_reads == False:
		# 	star_infile = '%s/%s/%s_no_ncRNA.fastq.gz' % (self.rootpath, self.read_proc_dir, self.samp)
		# 	star_out_path = '%s/%s/%s/STAR/%s' % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		# 	star_out_file = '%s/%s' % (star_out_path, self.samp)
		# 	if not os.path.exists(star_out_path): os.makedirs(star_out_path)
			
		# 	star_command = 'STAR --runThreadN %s --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterType BySJout --outFilterMultimapNmax 1 --outWigType wiggle read1_5p --outFileNamePrefix %s --quantMode TranscriptomeSAM --outReadsUnmapped Fastx' % (self.threadNumb, self.starGenome, star_infile, star_out_file)
		# 	print star_command
		# 	subprocess.Popen(star_command, shell=True).wait()
			
		# 	## sort and index bam outputs
		# 	## genome
		# 	index_genome = 'samtools index %sAligned.sortedByCoord.out.bam' % (star_out_file)
		# 	print index_genome
		# 	subprocess.Popen(index_genome, shell=True).wait()
		# 	## transcriptome
		# 	sort_transcriptome = 'samtools sort -@ %s %sAligned.toTranscriptome.out.bam %sAligned.toTranscriptome.out.sorted' % (self.threadNumb, star_out_file, star_out_file)
		# 	index_transcriptome = 'samtools index %sAligned.toTranscriptome.out.sorted.bam' % (star_out_file)
		# 	print sort_transcriptome
		# 	subprocess.Popen(sort_transcriptome, shell=True).wait()
		# 	print index_transcriptome
		# 	subprocess.Popen(index_transcriptome, shell=True).wait()
		# 	subprocess.Popen('rm %sAligned.toTranscriptome.out.bam' % (star_out_file), shell = True).wait()	


	def densebuilder_pairedreads(self, normalization_type):
		print "***** starting densebuilder run *****"

		assignment = "5" # "5" or "3" for fiveprime or threeprime alignment
		bam_inpath = "%s/%s/%s/STAR/%s" % (self.rootpath, self.align_dir, self.genome_name, self.samp)
		bamfileinput = "%s/%sAligned.sortedByCoord.out.bam" % (bam_inpath, self.samp)
		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')

		### testing broken LUC...

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
		if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

		# read_length_list = self.ftsize
		### for calculating rpm values 
		# if self.norm_type == "raw":
		if normalization_type == "raw":
			densitystring = "DensityUnnormalized"
			totreads = 1E6 # this effectively devides readcounts at every position by 1 so that data is not normalized
		# elif self.norm_type == "rpm":
		elif normalization_type == "rpm":
			densitystring = "Density_rpm"
			totreads = (bamfile.count()/2) # count the total number of reads in the bamfile and divide by this number (in millions) to get rpm
											# this is diveded by 2 because of paired end reads. 
		else:
			print "Normalization is not set for densebuilder!"
			sys.exit()
	

		## create a dictionary for input read size with no offset
		# riboshiftdict = "{"+str(readsize)+':[0]'+"}" # created in file instead
		bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.samp)
		bamfileoutput = '%s/%s_match.bam' % (bamfileoutpath, self.samp)
		
		## choose where to write density files:
		densityfilepath = '%s/%s/density5p_rnaseq/%s' %(fp_assign_path, densitystring, self.samp)
		densityfileout = '%s/%sf_' % (densityfilepath, self.samp)

		if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
		# removing riboshiftdict here
		commandstring = 'python2 %s/RNAseq/densebuilder_main_RNAseq.py --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (
			rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, self.threshold, totreads, densityfileout, bamfileoutput)
		print commandstring
		os.system(commandstring)

		sort_bamfile = 'samtools sort -@ %s %s %s/%s_match.sorted' % (self.threadNumb, bamfileoutput, bamfileoutpath, self.samp)
		print sort_bamfile
		subprocess.Popen(sort_bamfile, shell=True).wait()
		index_bamfile = 'samtools index %s/%s_match.sorted.bam' % (bamfileoutpath, self.samp)
		print index_bamfile
		subprocess.Popen(index_bamfile, shell=True).wait()
		clean_bamfile = 'rm %s' % (bamfileoutput)
		print clean_bamfile
		subprocess.Popen(clean_bamfile, shell=True).wait()

	def count_total_reads(self ):
		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		bamfileinput = '%s/%s_star_default/%s_match.sorted.bam' % (fp_assign_path, self.samp, self.samp)


		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')
		readnumber = bamfile.count()
		
		readcount_outfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.samp)
		readcountout = open(readcount_outfile, "w")
		readcountout.write(str(readnumber))
		readcountout.close()

	def build_countTables(self, fiveorthreeprime="5"):
		"""
		wrapper to direct building of density tables
		count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
		inset_choice contains adjustments to avoid start and stop codon peaks in data tables
		output is a csv file with densities for each region and 3'UTR occupancies
		"""
		
		assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.experiment, self.samp)
		outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store count tables
		if not os.path.exists(outfolder):   os.makedirs(outfolder)
		
		if self.norm_type == "raw":
			densitystring = "DensityUnnormalized"
		elif self.norm_type == "rpm":
			densitystring = "Density_rpm"
		else:
			print "Normalization is not set for densebuilder!"
			sys.exit()

		trspdictfilestring = '%s/%s/density%sp_rnaseq/%s/%sf_' %(
			fp_assign_path, densitystring, assignment, self.samp, self.samp)
		outfilestring = "%s/%s_%sp_rnaseq_%s_countTable" % (outfolder, self.samp, assignment, self.norm_type)
		print ""

		commandstring = 'python2 %s/RNAseq/RNAseq_buildDenseTables.py --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s' % (
			rootDir,trspdictfilestring,self.UTRfilestring,self.cdsDenThresh,self.norm_type,self.raw_dense_thresh,self.rpm_dense_thresh,self.inset_choice, outfilestring)
		print commandstring
		os.system(commandstring)



def main():
	for experiment in samples:
		for samp in samples[experiment]:
			RNAexp = RNAseqExperiment(samp, experiment, rootpath, fastqPath, threadNumb, genome_name, 
									ncRNAstarGenome, starGenome, linker, GTFfile, twobitfile,
									UTRfilestring, raw_dense_thresh, rpm_dense_thresh, norm_type, 
									inset_choice, paired_end_reads, threshold, read_proc_dir, align_dir)
			print RNAexp.cdsDenThresh
			RNAexp.remove_adapter()
			RNAexp.ncRNA_subtract_STAR()
			if multiMap == False:
				RNAexp.STAR_align()
			if multiMap == True:
				RNAexp.STAR_align_rMult()
			RNAexp.densebuilder_pairedreads(normalization_type = 'raw')
			RNAexp.densebuilder_pairedreads(normalization_type = 'rpm')
			RNAexp.count_total_reads()
			RNAexp.build_countTables()

if __name__ == '__main__':
	# execute only if run as a script
	main()
