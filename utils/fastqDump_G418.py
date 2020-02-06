### python script for downloading fastq files

import os, sys
import subprocess
from collections import OrderedDict
from pathos.multiprocessing import ProcessingPool as Pool
import argparse

### import libsettings file and add to global namespace
parser= argparse.ArgumentParser()
parser.add_argument('--rootDir', help= 'the root directory containing data and scripts')
parser.add_argument('--threadNumb', help= 'number of threads')
args = parser.parse_args()

rootDir = args.rootDir
threadNumb = str(args.threadNumb)


### 29 ribosome profiling samples to download and 12x2 RNA-seq samples


outdir_RP = '%s/Data/RPF/FASTQ' % rootDir
outdir_RNA = '%s/Data/RNA/FASTQ' % rootDir

if not os.path.exists(outdir_RP):	os.makedirs(outdir_RP)
if not os.path.exists(outdir_RNA):	os.makedirs(outdir_RNA)


srrList_RP = [
	'SRR10251273',
	'SRR10251274',
	'SRR10251275',
	'SRR10251276',
	'SRR10251277',
	'SRR10251278',
	'SRR10251279',
	'SRR10251280',
	'SRR10251281',
	'SRR10251282',
	'SRR10251283',
	'SRR10251284',
	'SRR10251285',
	'SRR10251286',
	'SRR10251287',
	'SRR10251288',
	'SRR10251289',
	'SRR10251290',
	'SRR10251299',
	'SRR10251300',
	'SRR10251301',
	'SRR10251302',
	'SRR10251303',
	'SRR10251304',
	'SRR10251305',
	'SRR10251306',
	'SRR10251307',
	'SRR10251308',
	'SRR10251309'
]

samplelist_RP = [
	"1_allAG_1_Untr_A",
	"2_allAG_2_Untr_B",
	"3_allAG_3_G418_500_A",
	"4_allAG_4_G418_500_B",
	"5_allAG_5_G418_2000_A",
	"6_allAG_6_G418_2000_B",
	"7_allAG_7_Genta_A",
	"8_allAG_8_Genta_B",
	"9_allAG_9_Parom_A",
	"10_allAG_10_Parom_B",
	"11_allAG_11_Neo_A",
	"12_allAG_12_Neo_B",
	"13_allAG_13_Tobra_A",
	"14_allAG_14_Tobra_B",
	"15_allAG_15_Amik_A",
	"16_allAG_16_Amik_B",
	"17_allAG_17_G418_500_10min_A",
	"18_allAG_18_G418_500_10min_B",
	"19_G418rep1_1_Untr_A",
	"20_G418rep1_2_Untr_B",
	"21_G418rep1_3_G418_A",
	"22_G418rep1_4_G418_B",
	"23_G418rep2_1_Untr_A", 
	"24_G418rep2_2_Untr_B", 
	"25_G418rep2_3_G418_A", 
	"26_G418rep3_1_Untr_A",
	"27_G418rep3_2_G418_A",
	"28_Calu6_1_Untr_A", 
	"29_Calu6_2_G418_A" 
]

RP_dict = OrderedDict(zip(srrList_RP, samplelist_RP))

srrList_RNA = [
	'SRR10251291',
	'SRR10251292',
	'SRR10251293',
	'SRR10251294',
	'SRR10251295',
	'SRR10251296',
	'SRR10251297',
	'SRR10251298',
	'SRR10251310',
	'SRR10251311',
	'SRR10251312',
	'SRR10251313'
]

samplelist_RNA = [
	"1_allAG_1_Untr_A",
	"2_allAG_2_Untr_B",
	"3_allAG_3_G418_500_A",
	"4_allAG_4_G418_500_B",
	"5_allAG_5_G418_2000_A",
	"6_allAG_6_G418_2000_B",
	"7_allAG_17_G418_10min_A",
	"8_allAG_18_G418_10min_B",
	"9_Calu6_1_Untr_A",
	"10_Calu6_2_Untr_B",
	"11_Calu6_3_G418_A",
	"12_Calu6_4_G418_B"
]

RNA_dict = OrderedDict(zip(srrList_RNA, samplelist_RNA))


def get_FASTQ_RP(srr):

	command_to_run = "fastq-dump --readids --split-files --skip-technical --gzip %s --outdir %s" % (srr, outdir_RP)
	print command_to_run
	os.system(command_to_run)

def get_FASTQ_RNA(srr):

	command_to_run = "fastq-dump --readids --split-files --skip-technical --gzip %s --outdir %s" % (srr, outdir_RNA)
	print command_to_run
	os.system(command_to_run)


def rename_RP_fastq():
	for srr in RP_dict:
		rename_cmnd = "mv %s/%s_1.fastq.gz %s/%s_RP.fastq.gz" % (
			outdir_RP, srr, outdir_RP, RP_dict[srr])
		os.system(rename_cmnd)

def rename_RNA_fastq():
	for srr in RNA_dict:
		rename_cmnd_1 = "mv %s/%s_1.fastq.gz %s/%s_RNA_1.fastq.gz" % (
			outdir_RNA, srr, outdir_RNA, RNA_dict[srr])
		os.system(rename_cmnd_1)

		rename_cmnd_2 = "mv %s/%s_2.fastq.gz %s/%s_RNA_2.fastq.gz" % (
			outdir_RNA, srr, outdir_RNA, RNA_dict[srr])
		os.system(rename_cmnd_2)

def main():

	print "starting fastq-dump"

	p = Pool(nodes=40)
	p.map(get_FASTQ_RP, srrList_RP)
	p.map(get_FASTQ_RNA, srrList_RNA)
	rename_RP_fastq()
	rename_RNA_fastq()


if __name__ == '__main__':
	main()






###