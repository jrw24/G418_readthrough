"""

"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines
##
import sys, os
import pandas as pd
import numpy as np
import seaborn as sns
import pysam
import argparse
import importlib
from collections import OrderedDict

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

# print samplelist
samplelist = [samplelist[22]]
# samplelist = samplelist[22]
# samplelist 
# samp = samplelist[22]

# plotoutdir = "%s/FPassignment/%s/%s/analysis/heatmaps" % (rootpath, genome_name, experiment)
# if not os.path.exists(plotoutdir):	os.makedirs(plotoutdir)


def build_tidy_df(position, folder, file):
    """
	build a tidy dataframe containing all samples in a single frame
	position must be 'start' or 'stop'
    """
    dframe = pd.DataFrame()
    df_list = []
    avgpath = "%s/FPassignment/%s/%s/%s" % (rootpath, genome_name, folder, file)
    if position == "start":
    	pospath = 'avggene1_ORFrpm_0shift_5uneq150'
    	avggennumb = '1'
    if position == "stop":
    	pospath = 'avggene2_ORFrpm_0shift_5uneq150'
    	avggennumb = '2'

    for i in ftsize[::-1]: # reverse order of ftsize to start with 40 and work backward to 15, possibly simplifies plotting
        readlen = str(i)
        df = pd.read_csv("%s/%s/%s_%sf_avg_%s.csv" % (avgpath, pospath, file, readlen, avggennumb))
        df['sample'] = file # create a new column with the sample name
        df['readlen'] = readlen # create another new column for the readlength being measured
        df_list.append(df) # add this to the empty placeholder list
    dframe = pd.concat(df_list, axis=0) # combine all dataframes into a single dataframe, combining by the rows - keeping same columns
    return dframe

def concat_tidy_dfs(position):
	for folder in samples:
		all_samples = [] # an empty list that will eventually contain all entries for all samples
		for file in samples[folder]:
			out_df = build_tidy_df(position, folder, file) # build a tidy data frame for each sample
			all_samples.append(out_df)
		df_all = pd.concat(all_samples, axis=0, ignore_index=True)
		if position == 'start':
			df_all['position'] = df_all['position']-50
		if position == 'stop':
			df_all['position'] = df_all['position']-150
	return df_all

def draw_heatmap(*args, **kwargs):
	data = kwargs.pop('data')
	d = data.pivot(index=args[1], columns=args[0], values=args[2]) # I'm not entirely sure why we need t
	ax = sns.heatmap(d, cmap="Blues", **kwargs)
	ax.invert_yaxis()
	major_xticks=ax.get_xticks()[::5]
	minor_xticks=ax.get_xticks()
	major_tick_labels=ax.get_xticklabels()[::5]
	ax.set_xticklabels(major_tick_labels)
	ax.set_xticks(minor_xticks, minor=True)
	ax.set_xticks(major_xticks, minor=False)
	major_yticks=ax.get_yticks()[::2]
	minor_yticks=ax.get_yticks()
	major_ytick_labels=ax.get_yticklabels()[::2]
	ax.set_yticklabels(major_ytick_labels)
	ax.set_yticks(minor_yticks, minor=True)
	ax.set_yticks(major_yticks, minor=False)


def meta_heatmap(position, df_all):
	fg = sns.FacetGrid(df_all, row="sample", aspect=4) # aspect sets the aspect ratio to be what we want
	fg.map_dataframe(draw_heatmap, 'position', 'readlen', 'avg', robust=False)
	if position == 'start':
		metapos = 'start'
	if position == 'stop':
		metapos = 'stop'
	
	plotoutfig_pdf = "%s/meta_heatmap_%s.pdf" % (plotoutdir, metapos)
	plt.savefig(plotoutfig_pdf)

def meta_heatmap_robust(position, df_all):
	fg = sns.FacetGrid(df_all, row="sample", aspect=4) # aspect sets the aspect ratio to be what we want
	fg.map_dataframe(draw_heatmap, 'position', 'readlen', 'avg', robust=True)
	if position == 'start':
		metapos = 'start'
	if position == 'stop':
		metapos = 'stop'
	plotoutfig_pdf = "%s/meta_heatmap_robust_%s.pdf" % (plotoutdir, metapos)
	plt.savefig(plotoutfig_pdf)



def get_single_size_distribution(samp):
	"""
	This is a script to build a data frame based on the total number of read counts for each read lengthi in ftsize
	Output is a data frame with number of mapped reads per read length
	"""

	### set footprint assignment path and the path to bamfiles used to count each read length
		## included read lengths are defined by ftsize
	fp_assign_path = '%s/FPassignment/%s/%s/%s' % (rootpath, genome_name, experiment, samp)
	bamfile_prefix = '%s/%s_star_default/%s' % (fp_assign_path, samp, samp)
	df = pd.DataFrame(columns = ['count'])
	# print df

	for read_size in ftsize:
		inbamfile = '%s_%s_match.sorted.bam' % (bamfile_prefix, read_size)
		bamfile = pysam.AlignmentFile(inbamfile, 'rb')
		read_count = bamfile.count()
		df.loc[read_size] = [read_count]

	total_reads = df['count'].sum()
	colname = "%s_percent" % (samp)
	df[colname] = (df['count'] / total_reads)*100
	return df

def main():
	poslist = ['stop']

	for pos in poslist:
		df_all = concat_tidy_dfs(pos)
		
		scWindow = range(-19,-16) ## the stop codon...
		windf = df_all[df_all['position'].isin(scWindow)]

		dl = []
		for samp in samplelist:
			tempdf = windf[windf['sample'] == samp]
			dl.append(tempdf)

		rldict = OrderedDict()
		totdense = 0
		for ft in ftsize:
			rldf = dl[0][dl[0]['readlen'] == str(ft)]
			totdense += rldf['avg'].sum()
			rldict[ft] = rldf['avg'].sum()

		perdict = OrderedDict()
		for rd in rldict:
			perc = rldict[rd]/totdense*100
			perdict[rd] = perc 

		perdf = pd.DataFrame.from_dict(perdict, orient = "index")

		plotoutfile = '%s/figures/Fig4C.pdf' % rootDir

		fullDfl = []
		for samp in samplelist:
			print samp
		full_df = get_single_size_distribution(samp)
		fullDfl.append(full_df)
		

		newdf = perdf.join(full_df) ### join merges on index
		newdf.drop('count',axis=1,inplace=True)

		print newdf

		fig, ax = plt.subplots(figsize=(3,6))

		newdf = newdf.loc[25:35,]
		print newdf
		print newdf.columns
		print newdf.index

		mycolors = ['red', 'black']
		# ax.set_xlim(xmin=25, xmax=35)
		for i in range(len(newdf.columns)):
		# readsize_df.plot.line(y=readsize_df.columns[i], style='.-', ax=ax, legend=True)
			newdf.plot.line(y=newdf.columns[i], style='-', ax=ax, color=mycolors[i], legend=True)
		# plt.xlim(xmin = min(ftsize), xmax = max(ftsize))
		# plt.xlim(xmin = 25, xmax = 35)
		# ax.set_xlim(xmin=, xmax=35)
		plt.xticks(newdf.index, fontsize=6)
		plt.ylabel('Percent of Reads')
		plt.xlabel('Read lengths')
		L=plt.legend(loc=2, prop={'size': 6})
		L.get_texts()[0].set_text('Stop Codon')
		L.get_texts()[1].set_text('All Reads')
		plt.savefig(plotoutfile, format='pdf', bbox_inches="tight")
		plt.close()



if __name__ == '__main__':
	# execute only if run as a script
	main()


