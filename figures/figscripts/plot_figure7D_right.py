#####


import matplotlib
matplotlib.use('Agg') ## set backend here
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



###

samplelist = [samplelist[0]]+ [samplelist[1]]
samples_plotted = '_vs_'.join(samples)

treatList = ['Untr', 'G418']

colorlist = ['black', 'orange'] #, 'green']
# colorlist = ['orange']

### trspname must be specified here,
	## multiple transcripts allowed per gene

genename = "GPX4"
trspname = "ENST00000611653.4"

# genename = "SELENOH"
# trspname = "ENST00000534355.6"
# # trspname = "ENST00000622257.4"


# genename = "TP53"
# trspname = "ENST00000445888.6"

# genename = 'AMD1'
# trspname = 'uc001gjk.3'

setYmax = True
yplotmax = 5.0 # must be a float...

# mRNAdf = pd.read_csv('%s/genomes/gencodeV30_all_tr_mRNAseqs.csv',
# 						sep=',', index_col=0)

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
	

densitystring = "Density_allTrDf"
	

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
	print 'retrieving density'

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

def retrieve_gene_density_DF(exp, samp, UTRdict, trspname):
	print "retreiving gene density from data.frame"
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, exp, samp)
	trDFinpath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_.csv.gz' %(
		fp_assign_path, densitystring, assignment, ribosome_shift, 
		minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen)

	### read in dataframe here:
	trDF = pd.read_csv(trDFinpath, compression='gzip', index_col=0)

	if UTRdict.has_key(trspname):
		exonsplicedcounts = ast.literal_eval(trDF.loc[trspname][0])
	else:
		print "transcript %s is not in UTR file" % trspname

	# print exonsplicedcounts

	gene_df = pd.DataFrame(exonsplicedcounts)
	gene_df.columns = [samp]
	return gene_df


def plot_genemodels_combined(df_all, UTRdict, colorlist, trspname):
	print 'starting plots'
	plotoutfolder = '%s/FPassignment/%s/%s/analysis/genemodels/%s' % (rootpath, genome_name, experiment, densitystring)
	# plotoutfolder = '/home/jamie/misc/histoneTest'
	if not os.path.exists(plotoutfolder):   os.makedirs(plotoutfolder)
	# plotoutfile = "%s/%s_%s_genemodel_overlay_combined.pdf" % (plotoutfolder, genename, pop)
	plotoutfile = "%s/figures/Fig7D_right.pdf" % (rootDir)
	
	
	mrnalen = int(UTRdict[trspname][3])
	cdslen = int(UTRdict[trspname][4])
	utr5len = int(UTRdict[trspname][5])
	utr3len = int(UTRdict[trspname][6])

	cdsstart = utr5len
	cdsend = mrnalen-utr3len
	
	# sc_zero = ast.literal_eval(stops.loc[trspname]['frameZeroStopPosMRNA'])

	# fig, axes = plt.subplots(figsize = (10,9))	
	fig, axes = plt.subplots(figsize = (6,4))
	# fig, axes = plt.subplots(figsize = (3.5,1.5))	
	
	for i in range(len(df_all.columns)):
		### set plot boundaries
		
		if setYmax == True:
			plot_max = yplotmax
		else:
			plot_max = df_all.max().max() # max of each col, then max of all cols

		## original
		# plot_height = plot_max+(plot_max/10) # 110% of the max rpm value
		# plot_min = (-plot_max/10)
		# gene_length = len(df_all)

		# utrmin = (plot_min/2)+(plot_min/10)
		# utrmax = plot_max/50
		# cdsmin = (plot_min/2)+(plot_min/4)
		# cdsmax = plot_max/20

		# ### alt gene model
		# plot_height = plot_max+(plot_max/10) # 110% of the max rpm value
		# plot_height = plot_max+(plot_max*0.10) # 110% of the max rpm value
		# # plot_min = (-plot_max*0.125) ## smaller gene model
		# plot_min = (-plot_max*0.2) ## bigger gene model
		# gene_length = len(df_all)

		# # utrmin = (plot_min/2)+(plot_min/10)
		# utrmax = plot_max*0.02
		# cdsmin = (plot_min*0.5)+(plot_min*0.25)
		# cdsmax = plot_max*0.05

		# utrmin = cdsmin*0.75
		# utrmax = -cdsmin*0.5

		# utrwidth = (-utrmin+utrmax)
		# utrwidAdj = utrwidth*0.1
		# min_scalor = plot_max*0.015

		### best gene model
		plot_height = plot_max+(plot_max/10) # 110% of the max rpm value
		plot_height = plot_max+(plot_max*0.10) # 110% of the max rpm value
		# plot_min = (-plot_max*0.125) ## smaller gene model
		plot_min = (-plot_max*0.2) ## bigger gene model
		gene_length = len(df_all)

		## set CDS
		cdsmin = (plot_min*0.75)
		cdsmax = -plot_min*0.5

		## set UTR
		utrmin = plot_min*0.625
		utrmax = -plot_min*0.25 # half height of cds

		upstream_bound = cdsstart-25
		downstream_bound = cdsend+50
		
		# print df_all.head()
		dfToPlot = df_all.iloc[:,i]

		# print dfToPlot.index
		dfToPlot.plot.line(y='rpm', use_index=True, ax = axes, c=colorlist[i], linewidth = 0.75)

		### for fill_between option:
		# dfToPlot.plot.line(y='rpm', use_index=True, ax = axes, c=colorlist[i], alpha=1, linewidth=0.5)
		# axes.fill_between(dfToPlot.index, dfToPlot, 0, color=colorlist[i], alpha=0.25)


		### draw the plot
		# df_all.iloc[:,i].plot.line(y='rpm', use_index=True, ax = axes, c=colorlist[i])
		rect_utr5 = matplotlib.patches.Rectangle((0,utrmin),utr5len,utrmax,edgecolor='black')
		rect_cds = matplotlib.patches.Rectangle((utr5len,cdsmin),cdslen,cdsmax,edgecolor='black')
		rect_utr3 = matplotlib.patches.Rectangle((utr5len+cdslen,utrmin),utr3len,utrmax,edgecolor='black')
	#         ptc = matplotlib.patches.Rectangle((593,cdsmin),3,cdsmax,edgecolor='red', fc = 'r') #PTC for R154X Nanoluc

		# ### larger rectangles
		# rect_utr5 = matplotlib.patches.Rectangle((0,utrmin),utr5len,utrmax,edgecolor='black')
		# rect_cds = matplotlib.patches.Rectangle((utr5len,cdsmin),cdslen,-cdsmin,edgecolor='black')
		# rect_utr3 = matplotlib.patches.Rectangle((utr5len+cdslen,utrmin),utr3len,utrmax,edgecolor='black')

		axes.axis([0, gene_length, plot_min, plot_height])
		# axes.axis([upstream_bound, downstream_bound, plot_min, plot_height])
		# axes.axis([0, 1500, plot_min, 2])
		# axes.axis([1250, 1750, plot_min, 10])
		axes.add_patch(rect_utr5)
		axes.add_patch(rect_cds)
		axes.add_patch(rect_utr3)
		# for stop in sc_zero:
		# 	s0 = matplotlib.patches.Rectangle((stop,utrmin),3,utrmax,edgecolor='red', fc = 'r',)
		# 	axes[0].add_patch(s0)
		if setYmax == True:
			axes.set_ylim(plot_min, yplotmax)
		axes.set_title("%s, %s" % (genename, trspname))
		# axes.axvline(722, linestyle="--", color='red')

	plt.savefig(plotoutfile, format='pdf')


def main():
	
	UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)
	if not trspname in UTRdict:

		print "%s not present in UTRdict" % trspname
		# continue
	dflist= []
	for samp in samplelist:
		print samp
		dfsamp = retrieve_gene_density_DF(experiment, samp, UTRdict, trspname)
		dflist.append(dfsamp)
	df_all = pd.concat(dflist, axis=1, ignore_index = False)
	# print df_all
	plot_genemodels_combined(df_all, UTRdict, colorlist, trspname)

if __name__ == '__main__':
	main()

