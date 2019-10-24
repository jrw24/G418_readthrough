### Script for merging fastq files from seperate experiments

import sys
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inputDir', help= 'directory with fastq files')
parser.add_argument('--outputDir', help = 'directory to send output')

args = parser.parse_args()

inpath = args.inputDir
outpath = args.outputDir

fq1 = [
	"1_allAG_1_Untr_A",
	"3_allAG_3_G418_500_A",
	"5_allAG_5_G418_2000_A",
	"7_allAG_7_Genta_A",
	"9_allAG_9_Parom_A",
	"11_allAG_11_Neo_A",
	"13_allAG_13_Tobra_A",
	"15_allAG_15_Amik_A",
	"17_allAG_17_G418_500_10min_A",
]

fq2 = [
	"2_allAG_2_Untr_B",
	"4_allAG_4_G418_500_B",
	"6_allAG_6_G418_2000_B",
	"8_allAG_8_Genta_B",
	"10_allAG_10_Parom_B",
	"12_allAG_12_Neo_B",
	"14_allAG_14_Tobra_B",
	"16_allAG_16_Amik_B",
	"18_allAG_18_G418_500_10min_B",
]

fq_merged = [

	"1_Untr",
	"2_G418_500",
	"3_G418_2k",
	"4_Genta",
	"5_Parom",
	"6_Neo",
	"7_Tobra",
	"8_Amik",
	"9_G418_500_10min"
]


# if not os.path.exists(FASTQpath): os.makedirs(FASTQpath)

def mergeFastQ(fq1Input, fq2Input, fqOutput):
	fq1 = '%s/%s*.fastq.gz' % (inpath, fq1Input)
	fq2 = '%s/%s*.fastq.gz' % (inpath, fq2Input)

	fqOut = '%s/%s.fastq.gz' % (outpath, fqOutput)

	merge_command = 'cat %s %s > %s' % (fq1, fq2, fqOut)
	print merge_command
	os.system(merge_command)

for sample in range(len(fq_merged)):
	mergeFastQ(fq1[sample], fq2[sample], fq_merged[sample])
