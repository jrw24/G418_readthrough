

"""
overlay mutliple genes ect...

"""


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

rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)
import rphelper as rph

### Inputs:

uORFinfile = "%s/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_uORFtable.csv" % rootDir

samples = [samplelist[18]]+ [samplelist[20]]
print samples
samples_plotted = '_vs_'.join(samples)
print samples_plotted

colorlist = ['black', 'orange'] #, 'green']
linewidthsplt = [1.0,0.75]


genename = "AMD1"

setYmax = False
# yplotmax = 8.0 # must be a float...



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

def get_trspname(genename):
	"""
	look up the trspname given a genename
	"""

	UTRfile = pd.read_csv(UTRfilestring, sep=',')

	trsp = UTRfile.loc[UTRfile['gene_name'] == genename]['#transcript'].item()
	return trsp


def find_trsp_uORFs(uORFinfile, trspname):
	"""
	look up the uORF entries for a given transcripts
	return a dataframe of all annotated uORFs here
	if there are no uORFs, return a df with no rows
	"""
	

	uORFdb = pd.read_csv(uORFinfile, index_col = 0)
	uORFdb.set_index('trxname', inplace=True)
	
	if trspname in uORFdb.index:
		uORFtr = uORFdb.loc[[trspname]]
	else:
		uORFtr = pd.DataFrame(columns=uORFdb.columns)
#     uORFtr
	return uORFtr


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

def plot_genemodels_combined(df_all, UTRdict, colorlist, trspname, uORFtr, genename):
	print 'starting plots'

	# plotoutfile = "%s/%s_%s_genemodel_overlay_combined.pdf" % (plotoutfolder, genename, pop)
	# plotoutfile = "%s/GMoverlay_uORFs_%s_%s_%s_%s.pdf" % (plotoutfolder, pop, genename, trspname, samples_plotted)
	plotoutfile = '%s/figures/Fig7E_left.pdf' % rootDir
	
	
	mrnalen = int(UTRdict[trspname][3])
	cdslen = int(UTRdict[trspname][4])
	utr5len = int(UTRdict[trspname][5])
	utr3len = int(UTRdict[trspname][6])

	cdsstart = utr5len
	cdsend = mrnalen-utr3len
	
	# sc_zero = ast.literal_eval(stops.loc[trspname]['frameZeroStopPosMRNA'])

	fig, axes = plt.subplots(figsize = (8,4))	
	# fig, axes = plt.subplots(figsize = (10,9))	
	# fig, axes = plt.subplots(nrows=1, ncols=1, sharex = False, figsize = (10,9))
	
	for i in range(len(df_all.columns)):
		### set plot boundaries
		# plot_max = df_all.max().max() # max of each col, then max of all cols
		# plot_height = plot_max+(plot_max/10) # 110% of the max rpm value
		# plot_min = (-plot_max/10)
		# gene_length = len(df_all)

		# utrmin = (plot_min/2)+(plot_min/10)
		# utrmax = plot_max/50
		# cdsmin = (plot_min/2)+(plot_min/4)
		# cdsmax = plot_max/20
		# utrwidth = (-utrmin+utrmax)
		# utrwidAdj = utrwidth/10


		plot_max = df_all.max().max() # max of each col, then max of all cols
		plot_height = plot_max+(plot_max*0.10) # 110% of the max rpm value
		plot_min = (-plot_max*0.125)
		gene_length = len(df_all)

		# utrmin = (plot_min/2)+(plot_min/10)
		utrmax = plot_max*0.02
		cdsmin = (plot_min*0.5)+(plot_min*0.25)
		cdsmax = plot_max*0.05

		utrmin = cdsmin*0.75
		utrmax = -cdsmin*0.5

		utrwidth = (-utrmin+utrmax)
		utrwidAdj = utrwidth*0.1

		min_scalor = plot_max*0.015
#         print utrwidth
		# print utrwidAdj
		
		cdsEndextend = utr5len + 50 
		if len(uORFtr) > 0:
			maxuORFextend = uORFtr.iloc[-1].startPosition + uORFtr.iloc[-1].uORFlen + 20
	#         print maxuORFextend, "this is maximal extension" 
		else:
			maxuORFextend = 0
		
		if maxuORFextend > cdsEndextend:
			maxXval = maxuORFextend
		else:
			maxXval = cdsEndextend
		
		# maxValCustom = 150
		# if maxValCustom != 0:
		#     maxXval = maxValCustom

		### draw the plot
		df_all.iloc[:,i].plot.line(y='rpm', use_index=True, ax = axes, c=colorlist[i], linewidth =linewidthsplt[i])
		rect_utr5 = matplotlib.patches.Rectangle((0,utrmin - min_scalor),utr5len,utrmax,edgecolor='black')
		rect_cds = matplotlib.patches.Rectangle((utr5len,cdsmin - min_scalor),cdslen,-cdsmin,edgecolor='black')
		rect_utr3 = matplotlib.patches.Rectangle((utr5len+cdslen,utrmin - min_scalor),utr3len,utrmax,edgecolor='black')
#         ptc = matplotlib.patches.Rectangle((593,cdsmin),3,cdsmax,edgecolor='red', fc = 'r') #PTC for R154X Nanoluc
		# axes.axis([0, gene_length, plot_min, plot_height])
		axes.axis([0, maxXval, plot_min, plot_height])
		axes.add_patch(rect_utr5)
		axes.add_patch(rect_cds)
		axes.add_patch(rect_utr3)


		if len(uORFtr) > 0:
			for i in range(len(uORFtr)):
				upORF = uORFtr.iloc[i] # save each row to a new object
			#     print upORF
#             (utrmin-utrwidAdj)
				rectuORFtotal = matplotlib.patches.Rectangle((upORF.startPosition, utrmin-utrmax/4 - min_scalor), upORF.uORFlen, utrmax+utrmax/2,
															  facecolor="orange", edgecolor="black") #, alpha=0.25)
#                 rectuORFtotal = matplotlib.patches.Rectangle((upORF.startPosition, utrmin-plot_min/10, upORF.uORFlen, (utrmax+utrmax/10),
#                                                               facecolor="orange", edgecolor="black")
				rectuORFstart = matplotlib.patches.Rectangle((upORF.startPosition, utrmin-utrmax/4 - min_scalor), 3, utrmax+utrmax/2, 
															  facecolor="green", edgecolor="black")
				rectuORFstop = matplotlib.patches.Rectangle((upORF.startPosition+upORF.uORFlen-3, utrmin-utrmax/4 -min_scalor), 3, utrmax+utrmax/2, 
															  facecolor="red", edgecolor="black")

				axes.add_patch(rectuORFtotal)
				axes.add_patch(rectuORFstart)
				axes.add_patch(rectuORFstop)
		axes.set_title("%s, %s" % (genename, trspname))

	plt.savefig(plotoutfile, format='pdf', bbox_inches="tight")
	plt.close()


def main():

	trspname = get_trspname(genename)

	UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)
	if not trspname in UTRdict:
		print "%s not present in UTRdict" % trspname
	uORFtr = find_trsp_uORFs(uORFinfile, trspname)
	print uORFtr
	dflist= []
	for samp in samples:
		print samp
		dfsamp = retrieve_gene_density(experiment, samp, UTRdict, trspname)
		dflist.append(dfsamp)
	df_all = pd.concat(dflist, axis=1, ignore_index = False)
	# print df_all
	plot_genemodels_combined(df_all, UTRdict, colorlist, trspname, uORFtr, genename)

if __name__ == '__main__':
	main()





