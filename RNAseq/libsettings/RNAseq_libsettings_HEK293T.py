"""
Inputs for RNAseq sample processing

"""

import os
os.path.join(os.path.dirname(__file__))

libsetDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
rootDir = libsetDir.replace('/RNAseq/libsettings','')


rootpath = "%s/Data/RNA" % rootDir ## path to data from root dir e.g.: /home/jamie/scripts/Wangen_G418/Data/RPF
fastqPath = rootpath+'/FASTQ' #Add path to fastq files, ideally in a directory named 'FASTQ'


#####

# Library input Settings:


experiment = 'HEK_G418'
samplelist = [
	"1_allAG_1_Untr_A",
	"2_allAG_2_Untr_B",
	"3_allAG_3_G418_500_A",
	"4_allAG_4_G418_500_B",
	"5_allAG_5_G418_2000_A",
	"6_allAG_6_G418_2000_B",
	"7_allAG_17_G418_10min_A",
	"8_allAG_18_G418_10min_B"
]

samples = {experiment:samplelist}

### Annotation Files: 
genome_name = 'hg38_protCode' 	# add this for storing alignment files to correct directory
ncRNAstarGenome = '%s/genomes/star_hg38_ncRNA' % rootDir # not sure if last backslash is necessary
starGenome = '%s/genomes/star_gtf_gencodeV30annotation' % rootDir
twobitfile = '%s/genomes/hg38.2bit' % rootDir
GTFfile = '%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs.gtf' % rootDir
UTRfilestring = '%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_UTRs.csv' % rootDir


### Alignment variables
paired_end_reads = True
multiMap = True
read_proc_dir = 'read_processing_hg38_multiRandom'
align_dir = 'alignments_hg38_multiRandom'

# set threads
threadNumb = '40' # number of threads for processing

linker = "NNNNNNCACTCGGGCACCAAGGAC" ### trim using skewer Illumina settings


### Footprint assignment
norm_type = "rpm" # either 'raw' for no normalization or 'rpm' for normalization to all mapped reads in alignment BAM file
threshold = '0' 
totreads = '1E6' 

### Count tables:
cdsDenThresh = False  # set to false for countTables
raw_dense_thresh = 0.01
rpm_dense_thresh = 0.001

# select insets to use to avoid start and stop peaks
inset_choice = 'default' # defaultInsets, zeroInsets, or customInsets



