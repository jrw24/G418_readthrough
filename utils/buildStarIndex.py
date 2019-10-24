__author__ = "jamie_wangen"

"""
build index for star algnments

"""

import sys, os
import argparse

### inputs
parser = argparse.ArgumentParser()
parser.add_argument('--rootDir', help = 'root directory')
parser.add_argument('--threadNumb', help = 'number of threads to run')
parser.add_argument('--STARsparsity', help = 'sparsity setting for star, default is 1 but set to higher value if RAM is an issue')
parser.add_argument('--genomeDir', help = 'directory to build STAR index')
parser.add_argument('--genomeFastaFiles', help = 'path to input fasta files for the genome')
parser.add_argument('--sjdbGTF', help = 'gtf file for splice junction annotations, set to "0" if there is none')
parser.add_argument('--SAindexNbases', help= 'length of SA pre-indexing string')

args = parser.parse_args()


threadNumb = args.threadNumb
sparsity = args.STARsparsity ### set to higher number if RAM is an issue
genomeDir = args.genomeDir # output directory to store index
genomeFastaFiles = args.genomeFastaFiles ## fasta file with the genome
sjdbGTFfile = args.sjdbGTF # annotation file in gtf format
SAindexNbases = args.SAindexNbases # default is 14, scale down to 9 for small genomes like ncRNA
# sjdbOverhang = 100 # default value


if not os.path.exists(genomeDir): os.makedirs(genomeDir) # make output directory if it doesnt exist

### ncRNA has 5537 entries, total sequence length of 580,856
### SAindexNbases:  log2(580856)/2 - 1 == 8.57

def build_star_index():
	"""

	"""

	if sjdbGTFfile == "0":
		### for ncRNA
		CMMD= 'STAR \
			--runThreadN %s \
			--runMode genomeGenerate \
			--genomeSAsparseD %s \
			--genomeDir %s \
			--genomeFastaFiles %s \
			--genomeSAindexNbases %s \ ' % (
				threadNumb, sparsity, genomeDir, genomeFastaFiles, SAindexNbases)
	else:
		### default
		CMMD= 'STAR \
			--runThreadN %s \
			--runMode genomeGenerate \
			--genomeSAsparseD %s \
			--genomeDir %s \
			--genomeFastaFiles %s \
			--genomeSAindexNbases %s \
			--sjdbGTFfile %s' % (
				threadNumb, sparsity, genomeDir, genomeFastaFiles, SAindexNbases, sjdbGTFfile)

	
	print CMMD
	os.system(CMMD)

build_star_index()