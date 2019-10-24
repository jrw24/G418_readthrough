"""
Input settings for ribosome profiling 

"""

import os

#####
### Sample inputs:
# rootpath = '/home/jamie/Data/G418/all/RPF'

libsetDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
rootDir = libsetDir.replace('/riboseq/libsettings','')


rootpath = "%s/Data/RPF" % rootDir ## path to data from root dir e.g.: /home/jamie/scripts/Wangen_G418/Data/RPF
fastqPath = rootpath+'/FASTQ/allAGmerge' #Add path to fastq files, ideally in a directory named 'FASTQ'

"""
experiments: 
	rpJRW002 - ata - 3
	rpJRW003 - cal - 2
	rpJRW009 - rep - 4
	rpJRW014 - ag3 - 2 ## omitting TIG treated samples
	rpJRW015 - ag4 - 8 ## including 10min timepoint

### need to change read processing settings for ata and cal libs


"""

experiment = 'allAGmerge'
samplelist = [
	"1_Untr",
	"2_G418_500",
	"3_G418_2k",
	"4_Genta",
	"5_Parom",
	"6_Neo",
	"7_Tobra",
	"8_Amik",
	"9_G418_500_10min"
]

trim5Plist = [1,1,1,1,1,1,1,1,1]

samples = {experiment:samplelist}

### Annotation Files: 
genome_name = 'hg38_protCode' 	# add this for storing alignment files to correct directory
ncRNAstarGenome = '%s/genomes/star_hg38_ncRNA' % rootDir # not sure if last backslash is necessary
starGenome = '%s/genomes/star_gtf_gencodeV30annotation' % rootDir
twobitfile = '%s/genomes/hg38.2bit' % rootDir
GTFfile = '%s/genomes/gencodeV30_all_tr.gtf' % rootDir
UTRfilestring = '%s/gencodeV30_all_tr_UTRs.csv' % rootDir
# stopcodons = '%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_stopcodons.csv' % rootDir


read_proc_dir = 'read_processing_hg38_multiRandom'
align_dir = 'alignments_hg38_multiRandom'
deduplicate_reads = True
### Alignment variables
# set threads
threadNumb = '40' # number of threads for processing
### Choose linker 1, 2, or 3
linker1 = "CTGTAGGCACCATCAAT" ## this is linker-1
linker2 = "NNNNNNCACTCGGGCACCAAGGAC" ## this is linker-2 (obz-407)

# linker = "CTGTAGGCACCATCAAT" ## this is linker-1
# linker = "NNNNNNCACTCGGGCACCAAGGAC" ## this is linker-2 (obz-407)
#linker = "NBAAAAAAAAAAAAAAAAAAAA" # linker 3 for polyA priming
# fivePrimeNTstoTrim = "4" # number of 5' nucleotides added by RT primer selection
# 						 # Trim 0 for old RT primer (ni-ni-9?), and 4 for new RT primer (obz-408)

### Footprint assignment
norm_type = "rpm" # either 'raw' for no normalization or 'rpm' for normalization to all mapped reads in alignment BAM file
threshold = '0' 

# ftsize is the iterable for pool.map(), it is good to keep this as a list
ftsize= [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]

### Ribosome populations
# Define readsizes to count for 'full-length' (fl), empty-Asite (eA), and empty-Esite (eE) ribosomes
# These vary emperically be RNase conditions
flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19

### Count tables:
cdsDenThresh = True
raw_dense_thresh = 0.01
rpm_dense_thresh = 0.001

# select insets to use to avoid start and stop peaks
inset_choice = 'default' # defaultInsets, zeroInsets, or customInsets


