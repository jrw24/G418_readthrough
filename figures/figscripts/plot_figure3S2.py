

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
import pandas as pd
import pysam
import numpy as np
import csv
from datetime import datetime
import argparse
from pathos.multiprocessing import ProcessingPool as Pool
pd.set_option('display.max_columns', 40)
import argparse
import importlib
from lifelines import KaplanMeierFitter

parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--libSetFile', help = 'libset file')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()
sys.path.append("%s/riboseq" % args.rootDir)
sys.path.append("%s/riboseq/libsettings" % args.rootDir)
import rphelper as rph

ootDir = args.rootDir
libsetName = args.libSetFile
libset = importlib.import_module("%s" % libsetName)
for attr in dir(libset):
	if not attr.startswith("_"):
		globals()[attr] = getattr(libset, attr)
threadNumb = str(args.threadNumb)



scdf = pd.read_csv(stopcodons)

def plot_inFrame_stops(scdf):

	scdict = {0:0,
		  1:0,
		  2:0,
		  3:0,
		  4:0,
		  5:0,
		  6:0,
		  7:0,
		  8:0,
		  9:0,
		  10:0,
		  11:0,
		  12:0,
		  13:0,
		  14:0,
		  15:0,
		  16:0,
		  17:0,
		  18:0,
		  19:0,
		  20:0,
		  21:0,
		  22:0,
		  23:0,
		  24:0,
		  25:0,
		  26:0,
		  27:0,
		  28:0,
		  29:0,
		  30:0,
		  31:0,
		  32:0,
		  33:0,
		  34:0,
		  35:0,
		  36:0,
		  37:0,
		  38:0,
		  39:0,
		  40:0
		 }

	for i in scdf.index:

		if scdf.loc[i, 'frameZeroStopCount'] >= 40:
			scdict[40] +=1
		else:
			scdict[scdf.loc[i, 'frameZeroStopCount']] +=1
			
	fig, ax = plt.subplots(figsize=(4,3))
	plt.bar(range(len(scdict)), list(scdict.values()), align='center')
	plt.savefig('%s/figures/Fig3S2A.pdf' % rootDir, format='pdf', bbox_inches="tight")

def plot_survival_function(scdf):
	
	dfl = scdf.copy()
	dfl['sc3'] = dfl['frameZeroUtr3LenAdj'] +3

	for i in dfl.index:
		if dfl.loc[i,'sc3'] < 101:
			dfl.loc[i,'kill'] = 1
		else:
			dfl.loc[i,'kill'] = 0

	kfm = KaplanMeierFitter()
	T = dfl['sc3']
	E = dfl['kill']

	kfm.fit(T, event_observed=E)
	kfm.survival_function_.plot()

	ax = plt.gca()

	ax.set_ylim(0,1)
	ax.set_xlim(0,100)

	figout = "%s/figures/Fig3S2B.pdf" % rootDir
	plt.savefig(figout, format ='pdf', bbox_inches="tight")


def main():

	plot_inFrame_stops(scdf)
	plot_survival_function(scdf)

if __name__ == '__main__':
	main()

		

		