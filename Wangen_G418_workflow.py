__author__ = 'Jamie Wangen'

"""
Workflow for generating figures for G418 paper: https://doi.org/10.1101/798579


### Dependencies:

Python:
matplotlib	2.2.2
numpy 	1.14.0
pandas	0.22.0
pysam	0.13
scipy	1.0.1
seaborn	0.9.0
twobitreader	3.1.5
argparse	1.1
multiprocessing	0.70a1
pathos	0.2.1


"""



import sys
import os
import importlib
import argparse
import subprocess
from ftplib import FTP

rootDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
threadNumb = 40



##### 1) Generate Genome Files:
	
	### download hg38 from genecode:

class generateGenomes(object):

	### building all files need for downstream analysis

	def __init__(self, rootDir, threadNumb):
		self.rootDir = rootDir
		self.threadNumb = threadNumb

	def download_hg38(self):
		"""

		"""
		outfile = "%s/genomes/GRCh38.primary_assembly.genome.fa.gz" % self.rootDir
		hg38file = "GRCh38.primary_assembly.genome.fa.gz"

		if not os.path.exists("%s/genomes" % rootDir):	os.makedirs("%s/genomes" % self.rootDir)

		ftp = FTP('ftp.ebi.ac.uk')
		ftp.login()
		ftp.cwd('pub/databases/gencode/Gencode_human/release_30/')

		with open(outfile, 'wb') as fp:
			ftp.retrbinary('RETR %s' % hg38file, fp.write)

		ftp.quit()


	def clean_hg38(self):

		GRCh38 = "%s/genomes/GRCh38.primary_assembly.genome.fa.gz" % self.rootDir

		# checkCmnd = "zcat %s | grep -in '>' " % GRCh38
		# print checkCmnd
		# subprocess.check_output(checkCmnd, shell=True)

		gunzipCmnd = "gunzip %s" % GRCh38
		subprocess.Popen(gunzipCmnd, shell= True).wait()

		### 51471479 last line of valid chromosomes
		cleanCmnd = "sed -n '1,51471479p' %s/genomes/GRCh38.primary_assembly.genome.fa > %s/genomes/hg38.fa" % (self.rootDir, self.rootDir)
		subprocess.Popen(cleanCmnd, shell= True).wait()

		indexCmnd = "samtools faidx %s/genomes/hg38.fa" % (self.rootDir)
		subprocess.Popen(indexCmnd, shell= True).wait()

		TwoBitCmnd = "faToTwoBit %s/genomes/hg38.fa %s/genomes/hg38.2bit" % (self.rootDir, self.rootDir)
		subprocess.Popen(TwoBitCmnd, shell= True).wait()

		removeGRCh38 = "rm %s/genomes/GRCh38.primary_assembly.genome.fa" % (self.rootDir)
		subprocess.Popen(removeGRCh38, shell= True).wait()


	def download_Gencode_annotation(self):

		outfile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)
		annotationFile = "gencode.v30.annotation.gtf.gz"

		outfile_tRNA = "%s/genomes/gencode.v30.tRNAs.gtf.gz" % (self.rootDir)
		annotationFiletRNA = "gencode.v30.tRNAs.gtf.gz"

		ftp = FTP('ftp.ebi.ac.uk')
		ftp.login()
		ftp.cwd('pub/databases/gencode/Gencode_human/release_30/')

		with open(outfile, 'wb') as fp:
			ftp.retrbinary('RETR %s' % annotationFile, fp.write)

		with open(outfile_tRNA, 'wb') as fp:
			ftp.retrbinary('RETR %s' % annotationFiletRNA, fp.write)		

		ftp.quit()

	def parse_GTF_file(self):

		gtfInFile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)

		parseCmnd = "python2 %s/utils/GTF_hg38_protCode_termStop_validUTR.py --gtfInFile %s --rootDir %s" % (
			self.rootDir, gtfInFile, self.rootDir)

		subprocess.Popen(parseCmnd, shell= True).wait()

	def build_annotation_files(self):

		gtfInFilePrefix = "%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs" % (self.rootDir)
		twoBitGenome = "%s/genomes/hg38.2bit" % (self.rootDir)

		mRNA_cmnd = "python2 %s/utils/mRNA_sequences_from_gtf.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(mRNA_cmnd, shell= True).wait()

		prot_cmnd = "python2 %s/utils/protein_sequences_from_gtf.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(prot_cmnd, shell= True).wait()

		utrTable_cmnd = "python2 %s/utils/makeUTRtable.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(utrTable_cmnd, shell= True).wait()

		stopFinder_cmnd = "python2 %s/utils/stopcodon_finder.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(stopFinder_cmnd, shell= True).wait()

		stopPositions_cmnd = "python2 %s/utils/utr3_stop_positions.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(stopPositions_cmnd, shell= True).wait()

		uORFs_cmnd = "python2 %s/utils/uORF_finder.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(uORFs_cmnd, shell= True).wait()

		codons_cmnd = "python2 %s/utils/codonFinder.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s --threadNumb %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome, self.threadNumb)
		subprocess.Popen(codons_cmnd, shell= True).wait()

	def parse_GTF_allTr(self):

		gtfInFile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)
		parseCmnd = "python2 %s/utils/GTF_hg38_all_tr.py --gtfInFile %s --rootDir %s" % (
			self.rootDir, gtfInFile, self.rootDir)

		subprocess.Popen(parseCmnd, shell= True).wait()

	def build_annotation_allTR(self):

		gtfInFilePrefix = "%s/genomes/gencodeV30_all_tr" % (self.rootDir)
		twoBitGenome = "%s/genomes/hg38.2bit" % (self.rootDir)

		utrTable_cmnd = "python2 %s/utils/makeUTRtable.py --gtfInFilePrefix %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, twoBitGenome)
		subprocess.Popen(utrTable_cmnd, shell= True).wait()



	def build_ncRNA_depletion(self):

		gtfInFile = "%s/genomes/gencode.v30.annotation.gtf.gz" % (self.rootDir)
		gtfInFileTrna = "%s/genomes/gencode.v30.tRNAs.gtf.gz" % (self.rootDir)
		rRNAncbi = "%s/genomes/refseq_human_rRNA_all.fa" % (self.rootDir)
		twoBitGenome = "%s/genomes/hg38.2bit" % (self.rootDir)

		ncRNA_depletion_cmnd = "python2 %s/utils/rRNA_depletion_hg38.py --gtfInFile %s --gtfInFileTrna %s --rRNAncbi %s --rootDir %s --twoBitGenome %s" % (
			self.rootDir, gtfInFile, gtfInFileTrna, rRNAncbi, self.rootDir, twoBitGenome)
		subprocess.Popen(ncRNA_depletion_cmnd, shell= True).wait()


	def build_STAR_indexes(self):

		### ncRNA
		ncSparsity = 1
		ncGenomeDir = "%s/genomes/star_hg38_ncRNA" % self.rootDir
		ncGenomeFasta = "%s/genomes/gencodeV30_ncRNA_all.fa" % self.rootDir
		ncSjdbGTF = "0"
		ncSAindexNbases = 9

		ncRNA_build_STAR_index = "python2 %s/utils/buildStarIndex.py --rootDir %s --threadNumb %s --STARsparsity %s --genomeDir %s --genomeFastaFiles %s --sjdbGTF %s --SAindexNbases %s" % (
			self.rootDir, self.rootDir, self.threadNumb, ncSparsity, ncGenomeDir, ncGenomeFasta, ncSjdbGTF, ncSAindexNbases)
		subprocess.Popen(ncRNA_build_STAR_index, shell=True).wait()

		### hg38
		## need to unzip GTF file to work with STAR
		unzip_cmnd = "gunzip %s/genomes/gencode.v30.annotation.gtf.gz" % self.rootDir
		subprocess.Popen(unzip_cmnd, shell=True).wait()


		hgSparsity = 1
		hgGenomeDir = "%s/genomes/star_gtf_gencodeV30annotation" % self.rootDir
		hgGenomeFasta = "%s/genomes/hg38.fa" % self.rootDir
		hgSjdbGTF = "%s/genomes/gencode.v30.annotation.gtf" % self.rootDir
		hgSAindexNbases = 14

		hg38_build_STAR_index = "python2 %s/utils/buildStarIndex.py --rootDir %s --threadNumb %s --STARsparsity %s --genomeDir %s --genomeFastaFiles %s --sjdbGTF %s --SAindexNbases %s" % (
			self.rootDir, self.rootDir, self.threadNumb, hgSparsity, hgGenomeDir, hgGenomeFasta, hgSjdbGTF, hgSAindexNbases)
		subprocess.Popen(hg38_build_STAR_index, shell=True).wait()

		gzip_cmnd = "gzip %s/genomes/gencode.v30.annotation.gtf" % self.rootDir
		subprocess.Popen(gzip_cmnd, shell=True).wait()



##### 2) Retrieve RawData

class RawData(object):

	### building all files need for downstream analysis

	def __init__(self, rootDir, threadNumb):
		self.rootDir = rootDir
		self.threadNumb = threadNumb


	def FASTQ_dump_sequences(self):
		"""
		Download all FASTQ_files in appropriate directories for downstream analysis
		"""

		fq_dump_cmnd = "python2 %s/utils/fastqDump_G418.py --rootDir %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.threadNumb)
		subprocess.Popen(fq_dump_cmnd, shell=True).wait()


	def merge_allAG_experiment(self):
		"""
		merge fastq files from allAG profiling experiment
		"""

		inputDir = "%s/Data/RPF/FASTQ" % self.rootDir
		outputDir = "%s/Data/RPF/FASTQ/allAGmerge" % self.rootDir

		if not os.path.exists(outputDir):	os.makedirs(outputDir)

		mergeFQ_cmnd = "python2 %s/utils/merge_fastq_allAG.py --inputDir %s --outputDir %s" % (
			self.rootDir, inputDir, outputDir)
		subprocess.Popen(mergeFQ_cmnd, shell=True).wait()

##### 3) Run main analysis pipeline on ribosome profiling data

class RibosomeProfiling_workflow(object):

	def __init__(self, rootDir, threadNumb, libSetFile, libSetFileAllTr):
		self.rootDir = rootDir
		self.threadNumb = threadNumb
		self.libSetFile = libSetFile
		self.libSetFileAllTr = libSetFileAllTr


	def RPexp(self):
		"""
		run the main ribosome profiling workflow for all samples here
		"""

		riboseq_cmnd = "python2 %s/riboseq/riboseq_main.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)

		subprocess.Popen(riboseq_cmnd, shell=True).wait()

	def RP_raw_countTables(self):
		"""
		build raw count tables for RNA seq files
		"""
		raw_countTables_cmnd = "python2 %s/riboseq/riboseq_build_exp_RAW_countTables.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(raw_countTables_cmnd, shell=True).wait()

	def RP_avgene_cdsNorm_start(self):

		avgene_cmnd = "python2 %s/riboseq/riboseq_avggene_cdsNorm_start.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(avgene_cmnd, shell=True).wait()

	def RP_avgene_cdsNorm_stop(self):

		avgene_cmnd = "python2 %s/riboseq/riboseq_avggene_cdsNorm_stop.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(avgene_cmnd, shell=True).wait()

	def RP_codon_occ(self):

		gtfInFilePrefix = "gencodeV30_protCode_TermStopCodon_validUTRs" 

		codon_occ_cmnd = "python2 %s/riboseq/riboseq_codon_occ_workflow.py --gtfInFilePrefix %s --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, gtfInFilePrefix, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(codon_occ_cmnd, shell=True).wait()

	def densebuild_allTr(self):

		dballTr_cmnd = "python2 %s/riboseq/riboseq_allTr_wf.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFileAllTr, self.threadNumb) 
		subprocess.Popen(dballTr_cmnd, shell=True).wait()



class RNAseq_workflow(object):

	def __init__(self, rootDir, threadNumb, libSetFile):
		self.rootDir = rootDir
		self.threadNumb = threadNumb
		self.libSetFile = libSetFile

	def RNAexp(self):

		rnaseq_cmnd = "python2 %s/RNAseq/RNAseq_main.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(rnaseq_cmnd, shell=True).wait()

	def RNA_raw_countTables(self):
		raw_countTables_cmnd = "python2 %s/RNAseq/RNAseq_build_exp_RAW_countTables.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			self.rootDir, self.rootDir, self.libSetFile, self.threadNumb)
		subprocess.Popen(raw_countTables_cmnd, shell=True).wait()


### plot figures:

class plot_figures(object):

	def __init__(self, rootDir, threadNumb, 
			libSet_RP_allG418, libSet_RP_allAGmerge,
			libSet_RNA_HEK293T, libSet_RNA_Calu6,
			libSet_RP_allG418_allTr):
		self.rootDir = rootDir
		self.threadNumb = threadNumb
		self.libSet_RP_allG418 = libSet_RP_allG418
		self.libSet_RP_allAGmerge = libSet_RP_allAGmerge
		self.libSet_RNA_HEK293T = libSet_RNA_HEK293T
		self.libSet_RNA_Calu6 = libSet_RNA_Calu6
		self.libSet_RP_allG418_allTr = libSet_RP_allG418_allTr


	def plot_figure_1(self):

		fig_cmnd = "python2 %s/figures/figscripts/plot_figure1.py --rootDir %s --threadNumb %s" % (
			rootDir, rootDir, threadNumb)
		subprocess.Popen(fig_cmnd, shell=True).wait()

	def plot_figure_2(self):

		fig_cmnd_AB = "python2 %s/figures/figscripts/plot_figure2AB.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allAGmerge, threadNumb)
		subprocess.Popen(fig_cmnd_AB, shell=True).wait()

		fig_cmnd_C = "python2 %s/figures/figscripts/plot_figure2C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_C, shell=True).wait()

		fig_cmnd_D = "python2 %s/figures/figscripts/plot_figure2D.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allAGmerge, threadNumb)
		subprocess.Popen(fig_cmnd_D, shell=True).wait()

		fig_cmnd_S1 = "python2 %s/figures/figscripts/plot_figure2S1.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1, shell=True).wait()

		fig_cmnd_S2A = "python2 %s/figures/figscripts/plot_figure2S2A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S2A, shell=True).wait()

		fig_cmnd_S2C = "python2 %s/figures/figscripts/plot_figure2S2C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S2C, shell=True).wait()

		fig_cmnd_S3A = "python2 %s/figures/figscripts/plot_figure2S3A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allAGmerge, threadNumb)
		subprocess.Popen(fig_cmnd_S3A, shell=True).wait()

		fig_cmnd_S3B = "python2 %s/figures/figscripts/plot_figure2S3B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S3B, shell=True).wait()


	def plot_figure_3(self):

		fig_cmnd_A = "python2 %s/figures/figscripts/plot_figure3A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_A, shell=True).wait()

		fig_cmnd_BC = "python2 %s/figures/figscripts/plot_figure3BC.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_BC, shell=True).wait()

		fig_cmnd_D = "python2 %s/figures/figscripts/plot_figure3D.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_D, shell=True).wait()

		fig_cmnd_E = "python2 %s/figures/figscripts/plot_figure3E.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_E, shell=True).wait()

		fig_cmnd_S1A = "python2 %s/figures/figscripts/plot_figure3S1A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1A, shell=True).wait()

		fig_cmnd_S1B = "python2 %s/figures/figscripts/plot_figure3S1B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1B, shell=True).wait()

		fig_cmnd_S1C = "python2 %s/figures/figscripts/plot_figure3S1C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1C, shell=True).wait()

		fig_cmnd_S2 = "python2 %s/figures/figscripts/plot_figure3S2.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S2, shell=True).wait()


	def plot_figure_4(self):

		fig_cmnd_A = "python2 %s/figures/figscripts/plot_figure4A.py --rootDir %s --threadNumb %s" % (
			rootDir, rootDir, threadNumb)
		subprocess.Popen(fig_cmnd_A, shell=True).wait()

		fig_cmnd_B = "python2 %s/figures/figscripts/plot_figure4B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_B, shell=True).wait()

		fig_cmnd_C = "python2 %s/figures/figscripts/plot_figure4C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_C, shell=True).wait()

		fig_cmnd_D = "python2 %s/figures/figscripts/plot_figure4D.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_D, shell=True).wait()

		fig_cmnd_S1A = "python2 %s/figures/figscripts/plot_figure4S1A.py --rootDir %s --threadNumb %s" % (
			rootDir, rootDir, threadNumb)
		subprocess.Popen(fig_cmnd_S1A, shell=True).wait()

		fig_cmnd_S1B = "python2 %s/figures/figscripts/plot_figure4S1B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1B, shell=True).wait()

	def plot_figure_5(self):

		fig_cmnd_A = "python2 %s/figures/figscripts/plot_figure5A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_A, shell=True).wait()

		fig_cmnd_B = "python2 %s/figures/figscripts/plot_figure5B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_B, shell=True).wait()

		fig_cmnd_C = "python2 %s/figures/figscripts/plot_figure5C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_C, shell=True).wait()

		fig_cmnd_D = "python2 %s/figures/figscripts/plot_figure5D.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_D, shell=True).wait()

		fig_cmnd_S1A = "python2 %s/figures/figscripts/plot_figure5S1A.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1A, shell=True).wait()

		fig_cmnd_S1B = "python2 %s/figures/figscripts/plot_figure5S1B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1B, shell=True).wait()

		fig_cmnd_S2 = "python2 %s/figures/figscripts/plot_figure5S2.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S2, shell=True).wait()

	def plot_figure_6(self):

		fig_cmnd_B = "python2 %s/figures/figscripts/plot_figure6B.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RNA_Calu6, threadNumb)
		subprocess.Popen(fig_cmnd_B, shell=True).wait()

		fig_cmnd_C = "python2 %s/figures/figscripts/plot_figure6C.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418_allTr, threadNumb)
		subprocess.Popen(fig_cmnd_C, shell=True).wait()

		fig_cmnd_S1 = "python2 %s/figures/figscripts/plot_figure6S1.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_S1, shell=True).wait()

	def plot_figure_7(self):

		libSetPathForR = "%s/RNAseq/libsettings/%s.py" % (rootDir, self.libSet_RNA_HEK293T)
		# print libSetPathForR

		fig_cmnd_AB_1 = "Rscript %s/figures/figscripts/plot_figure7AB.R %s %s" % (
			rootDir, libSetPathForR, rootDir)
		subprocess.Popen(fig_cmnd_AB_1, shell=True).wait()

		fig_cmnd_AB_2 = "python2 %s/figures/figscripts/plot_figure7AB.py --rootDir %s --threadNumb %s" % (
			rootDir, rootDir, threadNumb)
		subprocess.Popen(fig_cmnd_AB_2, shell=True).wait()

		fig_cmnd_Cl = "python2 %s/figures/figscripts/plot_figure7C_left.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allAGmerge, threadNumb)
		subprocess.Popen(fig_cmnd_Cl, shell=True).wait()

		fig_cmnd_Cr = "python2 %s/figures/figscripts/plot_figure7C_right.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allAGmerge, threadNumb)
		subprocess.Popen(fig_cmnd_Cr, shell=True).wait()

		fig_cmnd_Dl = "python2 %s/figures/figscripts/plot_figure7D_left.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418_allTr, threadNumb)
		subprocess.Popen(fig_cmnd_Dl, shell=True).wait()

		fig_cmnd_Dr = "python2 %s/figures/figscripts/plot_figure7D_right.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418_allTr, threadNumb)
		subprocess.Popen(fig_cmnd_Dr, shell=True).wait()

		fig_cmnd_El = "python2 %s/figures/figscripts/plot_figure7E_left.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_El, shell=True).wait()

		fig_cmnd_Er = "python2 %s/figures/figscripts/plot_figure7E_right.py --rootDir %s --libSetFile %s --threadNumb %s" % (
			rootDir, rootDir, self.libSet_RP_allG418, threadNumb)
		subprocess.Popen(fig_cmnd_Er, shell=True).wait()

		fig_cmnd_S1B = "python2 %s/figures/figscripts/plot_figure7S1B.py --rootDir %s --threadNumb %s" % (
			rootDir, rootDir, threadNumb)
		subprocess.Popen(fig_cmnd_S1B, shell=True).wait()



##### MAIN #####

def main():

	### 1) Generate Genome Files
	genomeGen = generateGenomes(rootDir, threadNumb)
	genomeGen.download_hg38()
	genomeGen.clean_hg38()
	genomeGen.download_Gencode_annotation()
	genomeGen.parse_GTF_file() ## this takes awhile
	genomeGen.build_annotation_files()
	genomeGen.parse_GTF_allTr() ## this takes a very long time
	genomeGen.build_annotation_allTR()
	genomeGen.build_ncRNA_depletion()
	genomeGen.build_STAR_indexes()

	### 2) Process Raw Data
	rawData = RawData(rootDir, threadNumb)
	rawData.FASTQ_dump_sequences()
	rawData.merge_allAG_experiment()

	### 3) Run Ribosome Profiling analysis pipeline
	RP = RibosomeProfiling_workflow(rootDir, threadNumb, 
		libSetFile="riboseq_libsettings_allG418",
		libSetFileAllTr="riboseq_libsettings_allG418_allTr")
	RP.RPexp()
	RP.RP_raw_countTables()
	RP.RP_avgene_cdsNorm_stop()
	RP.RP_codon_occ()
	RP.densebuild_allTr()

	RP2 = RibosomeProfiling_workflow(rootDir, threadNumb, 
		libSetFile="riboseq_libsettings_allAGmerge",
		libSetFileAllTr="riboseq_libsettings_allAGmerge_allTr")
	RP2.RPexp()
	RP2.RP_raw_countTables()
	RP2.RP_avgene_cdsNorm_start()
	RP2.RP_avgene_cdsNorm_stop()

	### 4) Run RNAseq analysis pipeline
	RNA = RNAseq_workflow(rootDir, threadNumb, libSetFile='RNAseq_libsettings_HEK293T')
	RNA.RNAexp()
	RNA.RNA_raw_countTables()

	RNA2 = RNAseq_workflow(rootDir, threadNumb, libSetFile='RNAseq_libsettings_Calu6')
	RNA2.RNAexp()
	RNA2.RNA_raw_countTables()

	### 5) Plot Figures
	pltFig = plot_figures(rootDir, threadNumb,
		libSet_RP_allG418 = 'riboseq_libsettings_allG418',
		libSet_RP_allAGmerge = 'riboseq_libsettings_allAGmerge',
		libSet_RNA_HEK293T = 'RNAseq_libsettings_HEK293T',
		libSet_RNA_Calu6 = 'RNAseq_libsettings_Calu6',
		libSet_RP_allG418_allTr = 'riboseq_libsettings_allG418_allTr')
	pltFig.plot_figure_1()
	pltFig.plot_figure_2()
	pltFig.plot_figure_3()
	pltFig.plot_figure_4()
	pltFig.plot_figure_5()
	pltFig.plot_figure_6()
	pltFig.plot_figure_7()



if __name__ == '__main__':
	main()





