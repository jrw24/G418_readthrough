

import sys
import os
import numpy as np
import argparse

from pathos.multiprocessing import ProcessingPool as Pool


### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--gtfInFilePrefix', help= 'the input annotation downloaded from gencode without ".gtf"')
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')

args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)

rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)

threadNumb = str(args.threadNumb)
GenAnPath = '%s/genomes' % args.rootDir
genAnName = "%s" % args.gtfInFilePrefix


### function inputs:

motiffilerootpath = '%s/%s_codons/%s_' % (GenAnPath, genAnName, genAnName)
# motiffilerootpath= '/home/jamie/genomes/human/hg38/gencodeV30_protCode_TermStopCodon_validUTRs_codons/gencodeV30_protCode_TermStopCodon_validUTRs_'

motifs= """ "['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGG','AGT','AGC','GGT','GGC','GGA','GGG','AGA','AGG','CGT','CGC','CGA','CGG','TAA','TAG','TGA']" """
cds5trim= '15'
cds3trim= '15'
seqwin= """ "['0','3']" """  # For density files shifted to the first nt of A sites, should be multiple of 3
exclusionfiles=  """ "['0']" """# check this


## define handler function:

def get_codon_occ_fl(file):
	"""
	function to get codon occupancies
	footprint_pop: fl, eA, eE, or aL
	ribosome_site: A, P, E, or 0 
	"""

	footprint_pop = 'fl'
	ribosome_site = 'A'
	fiveorthreeprime = '5'

	### set normalization type:
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	else:
		print "Normalization is not set!"
		sys.exit()

	### set footprint population:
	if footprint_pop == 'fl':
		readsizemin = str(flmin)
		readsizemax = str(flmax)
	elif footprint_pop == 'eA':
		readsizemin = str(eAmin)
		readsizemax = str(eAmax)
	elif footprint_pop == 'eE':
		readsizemin = str(eEmin)
		readsizemax = str(eEmax)
	elif footprint_pop == 'aL':
		readsizemin = 'aL'

	### density file paths:
	
	trspdictfilestring = '%s/FPassignment/%s/%s/%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (
		rootpath, genome_name, experiment, file, densitystring, fiveorthreeprime, ribosome_site, readsizemin, 
		readsizemax, file, readsizemin, readsizemax,file,readsizemin,readsizemax)
	
	sample_name = file

	outlistfilepath = '%s/FPassignment/%s/%s/%s/codon' % (rootpath, genome_name, experiment, file)
	if not os.path.exists(outlistfilepath): os.makedirs(outlistfilepath)
	outlistfile = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy.csv' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)
	outfileparams = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy_params.txt' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)

	commandstring= 'python2 %s/riboseq/codon_occupancy_main.py \
						--motiffilerootpath %s \
						--motifs %s \
						--trspdictfilestring %s \
						--sample_name %s \
						--UTRfilestring %s \
						--cds5trim %s \
						--cds3trim %s \
						--seqwin %s \
						--exclusionfiles %s \
						--outfileparams %s \
						--outlistfile %s' % (

						rootDir,
						motiffilerootpath, 
						motifs, 
						trspdictfilestring, 
						sample_name, 
						UTRfilestring, 
						cds5trim, 
						cds3trim, 
						seqwin, 
						exclusionfiles, 
						outfileparams, 
						outlistfile)
	print commandstring
	os.system(commandstring)


def get_codon_occ_eA(file):
	"""
	function to get codon occupancies
	footprint_pop: fl, eA, eE, or aL
	ribosome_site: A, P, E, or 0 although zero is not very informative here
	"""
	print "starting sample %s" % (file)

	footprint_pop = 'eA'
	ribosome_site = 'A'
	fiveorthreeprime = '5'

	### set normalization type:
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	else:
		print "Normalization is not set!"
		sys.exit()

	### set footprint population:
	if footprint_pop == 'fl':
		readsizemin = str(flmin)
		readsizemax = str(flmax)
	elif footprint_pop == 'eA':
		readsizemin = str(eAmin)
		readsizemax = str(eAmax)
	elif footprint_pop == 'eE':
		readsizemin = str(eEmin)
		readsizemax = str(eEmax)
	elif footprint_pop == 'aL':
		readsizemin = 'aL'

	# what is codonocc reading these in as?
	### density file paths:
	
	trspdictfilestring = '%s/FPassignment/%s/%s/%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (
		rootpath, genome_name, experiment, file, densitystring, fiveorthreeprime, ribosome_site, readsizemin, 
		readsizemax, file, readsizemin, readsizemax,file,readsizemin,readsizemax)
	
	sample_name = file

	outlistfilepath = '%s/FPassignment/%s/%s/%s/codon' % (rootpath, genome_name, experiment, file)
	if not os.path.exists(outlistfilepath): os.makedirs(outlistfilepath)
	outlistfile = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy.csv' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)
	outfileparams = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy_params.txt' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)

	commandstring= 'python2 %s/riboseq/codon_occupancy_main.py \
						--motiffilerootpath %s \
						--motifs %s \
						--trspdictfilestring %s \
						--sample_name %s \
						--UTRfilestring %s \
						--cds5trim %s \
						--cds3trim %s \
						--seqwin %s \
						--exclusionfiles %s \
						--outfileparams %s \
						--outlistfile %s' % (

						rootDir,
						motiffilerootpath, 
						motifs, 
						trspdictfilestring, 
						sample_name, 
						UTRfilestring, 
						cds5trim, 
						cds3trim, 
						seqwin, 
						exclusionfiles, 
						outfileparams, 
						outlistfile)

	print commandstring
	os.system(commandstring)


def main():

	p = Pool(nodes = 40)
	p.map(get_codon_occ_fl, samplelist)
	p2 = Pool(nodes = 40)
	p2.map(get_codon_occ_eA, samplelist)


if __name__ == '__main__':
	main()


