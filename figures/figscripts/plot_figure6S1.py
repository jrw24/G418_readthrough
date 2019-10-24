"""
This script is for taking input density files and plotting average genes around the start codon
	can be used to overlay multiple samples on the same plot

"""
## plotting
import matplotlib
matplotlib.use('Agg') ## set backend here
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


### inputs:
alignposition = "2" 	# set 1 for start codon and 2 for stop codon
pop = "fl" 	# fl, eA, or eE for density files
customSize = 30
threshold = 10
ribosome_site = "A" 	# A, P, E, or 0


samples = samplelist[27:29]


sample_plot_names = samples
samples_plotted = '_vs_'.join(samples)

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
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "whoops, something went wrong here, horribly horribly wrong!"
	# pop = ftsize



### for already shifted samples, for now
def avggene_riboshift_plot_overlay_zoom(alignposition, ribosome_site, experiment, fiveorthreeprime='5', normalizaiton='uneq', threshold ="0"):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalizaiton # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

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


	for file in samples:
		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, file)
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


		fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
		# eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
		# eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
		# aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

		if pop == "custom":
			custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
			print custom_avggene_csv
			custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
			df_custom_list.append(custom_avggene_df)

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		# eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		# eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		# aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		df_fl_list.append(fl_avggene_df)
		# df_eA_list.append(eA_avggene_df)
		# df_eE_list.append(eE_avggene_df)
		# df_aL_list.append(aL_avggene_df)

	plot_outfile = "%s/figures/Fig6S1.pdf" % (rootDir)


	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))


	ax.set_ylim(0, 0.3)
	ax.set_xlim(0,100)

	counter = 0 ## for making offsets
	if pop == 'fl':
		for i in range(len(df_fl_list)):
			df_fl_list[i].plot.line(x=df_fl_list[i].index+counter, y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
			# counter +=6
			# df_fl_list[i].plot.scatter(y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], s=2, use_index=True)
	if pop == 'eA':
		for i in range(len(df_eA_list)):
			df_eA_list[i].plot.line(y=df_eA_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
	if pop == 'custom':
		for i in range(len(df_fl_list)):
			df_custom_list[i].plot.line(x=df_custom_list[i].index+counter, y=df_custom_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])

	plt.legend(loc=1, prop={'size': 6})
	plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
	plt.close()


##########

def main():
	avggene_riboshift_plot_overlay_zoom(alignposition, ribosome_site, experiment, normalizaiton='eq', threshold=threshold)

if __name__ == '__main__':
	main()
