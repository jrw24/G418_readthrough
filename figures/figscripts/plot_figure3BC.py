
## plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import pandas as pd
import argparse
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
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)
sys.path.append('%s/riboseq' % rootDir)
import rphelper as rph




## inputs:

alignposition = "2" 	# set 1 for start codon and 2 for stop codon
pop = 'fl'
ribosome_site = "A" 	# A, P, E, or 0
normalizaiton = 'eq' # 'uneq' for rpm and 'eq' for rpkm

if normalizaiton == 'eq':
	threshold = '10'
else:
	threshold = '0'

folder = experiment
samples = samplelist[23:25]

figList = ['B', 'C']
treatList = ['Untr', 'G418']
# samples = samplelist
# samples = samplelist[1:3] 
sample_plot_names = samples

black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'
colorList = [black, orange, red, green, blue, yellow, pink]
col_list = colorList


if pop == "fl":
    minlen = str(flmin)
    maxlen = str(flmax)
elif pop == "eA":
    minlen = str(eAmin)
    maxlen = str(eAmax)
elif pop == "eE":
    minlen = str(eEmin)
    maxlen = str(eEmax)
elif pop == "custom":
	minlen = customSize
	maxlen = customSize
else:
    print "whoops, something went wrong here, horribly horribly wrong!"
    # pop = ftsize

samples_plotted = '_vs_'.join(samples)

### for already shifted samples, for now
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, folder, normalizaiton, fiveorthreeprime='5', threshold='0'):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

	# print "threshold ==", threshold

	sample_numb = len(samples)
	sample_index = range(0,len(samples))
	sample_dict = {}
	for i in sample_index:
		sample_dict[samples[i]] = sample_index[i]

	df_fl_list = []
	df_eA_list = []
	df_eE_list = []
	df_aL_list = []

	df_custom_list = []


	outpath = '%s/figures' % rootDir
	if not os.path.exists(outpath):   os.makedirs(outpath)

	print "using %s for normalization" % (norm)
	
	for file in samples:

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, folder, file) 
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now

	###

	## get paths to stored csv average gene files
		if pop == 'custom':
			# if norm == 'uneq':
			custom_csv = '%s/%s_%s_shiftCustom_rpkmThresh%s_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, threshold, customSize, customSize, alignpos)


		if norm == 'uneq':
			fl_avggene_csv = '%s/%s_fl_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)
			eA_avggene_csv = '%s/%s_eA_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, eAmin, eAmax, alignpos)
			eE_avggene_csv = '%s/%s_eE_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, eEmin, eEmax, alignpos)
			aL_avggene_csv = '%s/%s_rpkmThresh0_aL_aLf_avg_%s.csv' % (avggene_csv_path, file, alignpos)

		elif norm == 'eq':
			fl_avggene_csv = '%s/%s_fl_rpkmThresh10_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)
			eA_avggene_csv = '%s/%s_eA_rpkmThresh10_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, eAmin, eAmax, alignpos)
			eE_avggene_csv = '%s/%s_eE_rpkmThresh10_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, eEmin, eEmax, alignpos)
			aL_avggene_csv = '%s/%s_rpkmThresh10_aL_aLf_avg_%s.csv' % (avggene_csv_path, file, alignpos)

		else:
			print "no norm selected"

		

		# fl_avggene_csv = '%s/%s_fl_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)
		# eA_avggene_csv = '%s/%s_eA_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, eAmin, eAmax, alignpos)
		# eE_avggene_csv = '%s/%s_eE_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, eEmin, eEmax, alignpos)
		# aL_avggene_csv = '%s/%s_aL_aLf_avg_%s.csv' % (avggene_csv_path, file, alignpos)


		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		if pop == 'custom':
			custom_avggene_df = pd.read_csv(custom_csv, index_col = 0, header = 0)
			df_custom_list.append(custom_avggene_df)

		df_fl_list.append(fl_avggene_df)
		df_eA_list.append(eA_avggene_df)
		df_eE_list.append(eE_avggene_df)
		df_aL_list.append(aL_avggene_df)


	#######
	## reading frame analysis, for default -150, +100 reading frame
	counter = 0
	# for avdf in df_fl_list:
	if pop == 'fl':
		dfFrame = df_fl_list
	elif pop == 'eA':
		dfFrame = df_eA_list
	elif pop == 'eE':
		dfFrame = df_eE_list
	elif pop == 'custom':
		dfFrame = df_custom_list
	else:
		print "pop not set, no plot was made"

	for avdf in dfFrame:
		# print avdf

		### CDS - count regions from nt 0 to 130 
		### UTR3 - major stop codon peak is at postion 148, 
			### "0" is actually the position after the stop codon... need to doublecheck Colins avggene

		frame_minusOne = 0
		frame_zero = 0
		frame_plusOne = 0

		rpfCounts = avdf['avg']
		cdsCounts = rpfCounts[0:132]

		utr3Counts = rpfCounts[153:-1] ## do not include final position

		cds_frame_minusOne_counts = sum(cdsCounts[0::3])
		cds_frame_zero_counts = sum(cdsCounts[1::3])
		cds_frame_plusOne_counts = sum(cdsCounts[2::3])

		cdstotal = cds_frame_plusOne_counts+cds_frame_zero_counts+cds_frame_minusOne_counts

		utr3_frame_minusOne_counts = sum(utr3Counts[0::3])
		utr3_frame_zero_counts = sum(utr3Counts[1::3])
		utr3_frame_plusOne_counts = sum(utr3Counts[2::3])

		utr3total = utr3_frame_plusOne_counts+utr3_frame_zero_counts+utr3_frame_minusOne_counts

		cdsFrame = [
			(cds_frame_minusOne_counts/cdstotal)*100,
			(cds_frame_zero_counts/cdstotal)*100,
			(cds_frame_plusOne_counts/cdstotal)*100
		]

		utr3Frame = [ 
			(utr3_frame_minusOne_counts/utr3total)*100,
			(utr3_frame_zero_counts/utr3total)*100,
			(utr3_frame_plusOne_counts/utr3total)*100
		]


		### plotting separated by region, and not frame

		frames = [cdsFrame, utr3Frame]
		avgRegions = ['CDS', 'utr3']

		dfplt = pd.DataFrame(frames, index=avgRegions)

		figout = "%s/figures/Fig3%s.pdf" % (rootDir, figList[counter])
		fig,ax = plt.subplots(figsize=(5,5))
		ax.set_ylim(0,60)
		dfplt.plot.bar(ax=ax)
		ax.axhline(y=33.333, linestyle = "--")
		ax.set_title(treatList[counter])
		fig.savefig(figout, format='pdf')
		plt.close()

		counter +=1





##########


avggene_riboshift_plot_overlay(alignposition, ribosome_site, folder, normalizaiton, threshold=threshold)
