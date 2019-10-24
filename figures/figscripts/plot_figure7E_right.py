

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines


import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import struct
import ast
import pandas as pd
import argparse
import importlib

parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help= 'riboseq libsettings file to run riboseq_main')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)
import rphelper as rph

rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)


### Inputs:

stopPositionsPath = '%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_utr3StopPositions.csv' % rootDir
mRNAdf = pd.read_csv('%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_mRNAseqs.csv' % rootDir,
						sep=',', index_col=0)
AAdf = pd.read_csv('%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_Protseqs.csv' % rootDir,
						sep=',', index_col=0)

###

samples = [samplelist[18]]+[samplelist[20]]
samples_plotted = '_vs_'.join(samples)

colorlist = ['black','orange'] #, 'green']
# colorlist = ['orange', 'black'] #, 'green']
linewidthsplt = [1.0, 0.75]
# print colorlist

genelist = ['AMD1']

pop = "fl"
assignment = "5"
ribosome_shift = "A"

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
else:
	print "oops"
	
	
### Set normalization used: default is Density_rpm
if norm_type == "raw":
	densitystring = "DensityUnnormalized"
elif norm_type == "rpm":
	densitystring = "Density_rpm"
else:
	print "Normalization is not set for densebuilder!"
	sys.exit()



def load_stops_df(stopPositionsPath):
	stops = pd.read_csv(stopPositionsPath, index_col=0)
	return stops

def get_trspname(genename):
	"""
	look up the trspname given a genename
	"""
	# kgXref = pd.read_csv('/home/jamie/genomes/human/hg19/gene_annotation/hg19_kgXref.txt', sep = '\t', dtype='string')

	UTRfile = pd.read_csv(UTRfilestring, sep=',')

	trsp = UTRfile.loc[UTRfile['gene_name'] == genename]['#transcript'].item()
	return trsp


def open_single_dense_file (infilestring):
	transcript_reads = []
	with open(infilestring, 'rb') as f:
		nextval = f.read(4)
		while nextval != "":
			transcript_reads.append(struct.unpack("f",nextval)[0])
			nextval = f.read(4)
	return transcript_reads

def load_genomes(UTRfilestring, twobitfile):
	"""
	make this a separate function so that these only need to be loaded a single time
	"""
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	genome= twobitreader.TwoBitFile(twobitfile) # do we actually need to load this in here?
	return UTRdict, genome


### retrieve gene density file for a given sample 
def retrieve_gene_density(exp, samp, UTRdict, trspname):


	## set path to density files
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, exp, samp)
	trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_%s' %(
		fp_assign_path, densitystring, assignment, ribosome_shift, 
		minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen, trspname)

	if UTRdict.has_key(trspname):
		exonsplicedcounts = open_single_dense_file(trspfilestring)
	else:
		print "transcript %s is not in UTR file" % trspname
# 		sys.exit()

	gene_df = pd.DataFrame(exonsplicedcounts)
	gene_df.columns = [samp]
	return gene_df


def plot_genemodel_firstStop_Utr3(df_all, UTRdict, colorlist, trspname, stops, genename):


	plotoutfile = "%s/figures/Fig7E_right.pdf" % rootDir
	
	mrna_tr = mRNAdf.loc[trspname]['mRNAseqs']
	
	mrnalen = int(UTRdict[trspname][3])
	cdslen = int(UTRdict[trspname][4])
	utr5len = int(UTRdict[trspname][5])
	utr3len = int(UTRdict[trspname][6])

	cdsstart = utr5len
	cdsend = mrnalen-utr3len
	utr3seq = mrna_tr[cdsend:].upper()

	sc_peak_height = df_all.iloc[cdsend-12:].max().max()
	plot_max = sc_peak_height

	sc_zero = ast.literal_eval(stops.loc[trspname]['frameZeroStopPosMRNA'])
	sc_minusOne = ast.literal_eval(stops.loc[trspname]['frameMinusOneStopPosMRNA'])
	sc_plusOne = ast.literal_eval(stops.loc[trspname]['framePlusOneStopPosMRNA'])

	fig, axes = plt.subplots(figsize = (6,3))	

	
	for i in range(len(df_all.columns)):

		plot_height = plot_max+(plot_max/10) # 110% of the max rpm value
		plot_height = plot_max+(plot_max*0.10) # 110% of the max rpm value
		plot_min = (-plot_max*0.125)
		gene_length = len(df_all)

		utrmax = plot_max*0.02
		cdsmin = (plot_min*0.5)+(plot_min*0.25)
		cdsmax = plot_max*0.05

		utrmin = cdsmin*0.75
		utrmax = -cdsmin*0.5

		utrwidth = (-utrmin+utrmax)
		utrwidAdj = utrwidth*0.1

		min_scalor = plot_max*0.015


		upstream_bound = cdsend - 75
		downstream_bound = cdsend + 425
		
		if downstream_bound > gene_length:
			downstream_bound = gene_length

		dfToPlot = df_all.iloc[:,i]
		dfToPlot.plot.line(y='rpm', use_index=True, ax = axes, c=colorlist[i], linewidth = 0.75)


		rect_utr5 = matplotlib.patches.Rectangle((0,utrmin),utr5len,utrmax,edgecolor='black')
		rect_cds = matplotlib.patches.Rectangle((utr5len,cdsmin),cdslen,-cdsmin,edgecolor='black')
		rect_utr3 = matplotlib.patches.Rectangle((utr5len+cdslen,utrmin),utr3len,utrmax,edgecolor='black')

		axes.axis([upstream_bound, downstream_bound, plot_min, 6.0])
		axes.add_patch(rect_utr5)
		axes.add_patch(rect_cds)
		axes.add_patch(rect_utr3)

		### draw stop codons in each frame
		# for stop in sc_minusOne:
		# 	s0 = matplotlib.patches.Rectangle((stop,utrmin),3,utrmax,edgecolor='black', fc = 'gray',)
		# 	axes.add_patch(s0)

		# for stop in sc_plusOne:
		# 	s0 = matplotlib.patches.Rectangle((stop,utrmin),3,utrmax,edgecolor='black', fc = 'gray',)
		# 	axes.add_patch(s0)

		for stop in sc_zero:
			s0 = matplotlib.patches.Rectangle((stop,utrmin),3,utrmax,edgecolor='black', fc = 'r',)
			axes.add_patch(s0)

		axes.set_title("%s, %s" % (genename, trspname))
	plt.savefig(plotoutfile, format='pdf')



def main():
	for gene in genelist:

		# global trspname
		# try:
		# 	trspname
		# except NameError:
		# 	trspname_exists = False
		# else:
		# 	trspname_exists = True

		# if trspname_exists == False:
		# 	trspname = get_trspname(gene)
		trspname = get_trspname(gene)
		stops = load_stops_df(stopPositionsPath)

		UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)
		if not trspname in UTRdict:

			print "%s not present in UTRdict" % trspname
			# continue
		dflist= []
		for samp in samples:

			dfsamp = retrieve_gene_density(experiment, samp, UTRdict, trspname)
			dflist.append(dfsamp)
		df_all = pd.concat(dflist, axis=1, ignore_index = False)

		plot_genemodel_firstStop_Utr3(df_all, UTRdict, colorlist, trspname, stops, genename = gene)



if __name__ == '__main__':
	main()



