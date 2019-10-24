"""
script is based off of riboseq_readingFrame_utr3meta
"""
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
pd.set_option('display.max_columns', 50)
pd.set_option('display.max_rows', 300)
import seaborn as sns
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
import rphelper as rph

rootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)


### inputs:
alignposition = "2" 	# set 1 for start codon and 2 for stop codon
pop = 'fl'
ribosome_site = "A" 	# A, P, E, or 0
normalization = 'eq' # 'uneq' for rpm and 'eq' for rpkm

folder = experiment
samples = samplelist[0:4]+samplelist[6:16]
sample_plot_names = samples

### set colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'

colorList = [black, orange, red, green, blue, yellow, pink]


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


samples_plotted = '_vs_'.join(samples)

### for already shifted samples, for now
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, folder, normalization, fiveorthreeprime='5'):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalization # 'uneq' for no normalization; 'eq' to give equal weight to all genes
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

	# outpath = '%s/FPassignment/%s/%s/analysis/metaplots' % (rootpath, genome_name, folder)
	outpath = '%s/figures' % rootDir
	if not os.path.exists(outpath):   os.makedirs(outpath)
	
	for file in samples:

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, folder, file)
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now

	###

	## get paths to stored csv average gene files
		if pop == 'custom':
			if norm == 'uneq':
				custom_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, customSize, customSize, alignpos)

		if norm == 'uneq':
			fl_avggene_csv = '%s/%s_fl_rpkmThresh0_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)
			# eA_avggene_csv = '%s/%s_eA_rpkmThresh0_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, eAmin, eAmax, alignpos)
			# eE_avggene_csv = '%s/%s_eE_rpkmThresh0_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, eEmin, eEmax, alignpos)
			# aL_avggene_csv = '%s/%s_rpkmThresh0_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, alignpos)

		elif norm == 'eq':
			fl_avggene_csv = '%s/%s_fl_rpkmThresh10_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, flmin, flmax, alignpos)
			# eA_avggene_csv = '%s/%s_eA_rpkmThresh10_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, eAmin, eAmax, alignpos)
			# eE_avggene_csv = '%s/%s_eE_rpkmThresh10_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, eEmin, eEmax, alignpos)
			# aL_avggene_csv = '%s/%s_rpkmThresh10_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, alignpos)

		else:
			print "no norm selected"



		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		# eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		# eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		# aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		# if pop == 'custom':
		# 	custom_avggene_df = pd.read_csv(custom_csv, index_col = 0, header = 0)
		# 	df_custom_list.append(custom_avggene_df)

		df_fl_list.append(fl_avggene_df)
		# df_eA_list.append(eA_avggene_df)
		# df_eE_list.append(eE_avggene_df)
		# df_aL_list.append(aL_avggene_df)


	#######

	counter = 0
	# for avdf in df_fl_list:
	relDenList = []

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

		# print avdf.loc[-148:-17] ## not zero based with loc
		# print len(avdf.loc[-148:-17])

		# print avdf.loc[5:100, 'avg'] ## not zero based with loc
		# print len(avdf.loc[5:100])

		### CDS - count regions from nt -148 to -17 
		### UTR3 - countr regiona from 5 to 100, inclusive


		rpfCounts = avdf['avg']

		cdsCounts = avdf.loc[-148:-17, 'avg']
		cdsDense = (sum(cdsCounts)/len(cdsCounts))

		utr3Counts = avdf.loc[5:100, 'avg'] ## do not include final position
		utr3Dense = (sum(utr3Counts)/len(utr3Counts))

		relDense = (utr3Dense/cdsDense)*100
		relDenList.append(relDense)

		counter +=1


	figout = "%s/figures/Fig2C.pdf" % (rootDir)
	fig,ax = plt.subplots(figsize=(5,5))
	dfplt = pd.DataFrame(relDenList, index=samples, columns=['utrPer'])
	dfplt['treatments'] = ['untr', 'untr', 'g418','g418', 'gent','gent', 'par','par', 'neo','neo', 'tobra','tobra', 'amik', 'amik']
	dfplt['colorVal'] = [black, black, orange, orange, red, red, green, green, blue, blue, yellow, yellow, pink, pink]


	new_ind = range(len(dfplt))[::2]

	av_vals = []
	for i in new_ind:

		i = int(i)

		v1 = dfplt.iloc[i,0]
		v2 = dfplt.iloc[i+1,0]

		av =(v1+v2)/2

		av_vals.append(av)



	dfav = pd.DataFrame.from_dict({
			"tr":dfplt['treatments'][::2],
			'av':av_vals
		})


	### seaborn catplots:
	fig, ax = plt.subplots(figsize=(6,6))
	sns.boxplot(data=dfav, x='tr', y='av', showbox=False , width = 0.5, 
					 showcaps=False, color = "black")
	sns.swarmplot(data=dfplt, x='treatments', y='utrPer', size=8, ax=ax, palette=colorList)
	ax=plt.gca()

	for item in ax.get_xticklabels():
		item.set_rotation(90)
	ax.set_ylim(0, 25)
	plt.savefig(figout, format='pdf', bbox_inches = "tight")


##########

def main():
	avggene_riboshift_plot_overlay(alignposition, ribosome_site, folder, normalization)

if __name__ == '__main__':
	main()

