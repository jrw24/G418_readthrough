__author__ = 'Jamie Wangen'

"""
workflow for processing a custom readsize that falls outside of fl, eA, eE, or aL designations
	fl == 'full length'; typically 28-35
	eA == 'empty Asite'; typically 21-24
	eE == 'empty Esite'; 18-19 - hypothesized reads around start codon
	aL == 'all Lengths'; combination of all reads from above after respecitive Asite shifts

This can be run after read processing > genome mapping > densebuilder unshifted > riboshift
	riboshift is necessary to get offset for single custom readlength

Next, build density for A shift of readlength

"""

import sys, os 
import pysam
import time
import pandas as pd
import collections
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
from collections import OrderedDict
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




print samples
customSize = 30 ## necessary, but not used unless specifying a custom read length



### make this a class
class RP_custom_readlength(object):

	### instantiate the required inputs for the experiment

	def __init__(self,file,folder,rootpath,align_dir,genome_name,customSize,
				 flmin,flmax,eAmin,eAmax,eEmin,eEmax,
				 GTFfile,twobitfile,threshold,UTRfilestring,norm_type):

		self.file = file
		self.folder = folder
		self.rootpath = rootpath
		self.align_dir = align_dir
		self.genome_name = genome_name
		self.customSize = customSize
		self.flmin = flmin
		self.flmax = flmax
		self.eAmin = eAmin
		self.eAmax = eAmax
		self.eEmin = eEmin
		self.eEmax = eEmax
		self.GTFfile = GTFfile
		self.twobitfile = twobitfile
		self.threshold = threshold
		self.UTRfilestring = UTRfilestring
		self.norm_type = norm_type


	def densebuilder_riboshift(self, ribosome_site):
		"""
		Function to map densities with no shift and combine all read lengths
		could be taken out of this class for use in multiprocessing possibly
		"""

		assignment = "5" # should be five prime assignment to use shift dicts
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
		
		bam_inpath = "%s/%s/%s/STAR/%s" % (self.rootpath, self.align_dir, self.genome_name, self.file)
		bamfileinput = "%s/%sAligned.sortedByCoord.out.bam" % (bam_inpath, self.file)
		# bamfile = pysam.AlignmentFile(bamfileinput, 'rb')

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.folder, self.file)
		bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.file)
		if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

		densitystring = "Density_allTrDf"
		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')
		totreads = bamfile.count() ## simply count the number of reads in the bamfile to use for normalization


		# Set paths to riboshift dicts for each species of interest
		fl_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eA_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		eE_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)


		fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
		fl_ribo_shift.columns = ['riboshift']
		eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
		eA_ribo_shift.columns = ['riboshift']
		eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
		eE_ribo_shift.columns = ['riboshift']

		# print eE_ribo_shift.index.tolist()
		fl_sizes = fl_ribo_shift.index.tolist()
		eA_sizes = eA_ribo_shift.index.tolist()
		eE_sizes = eE_ribo_shift.index.tolist()

		fl_shiftdict = collections.OrderedDict()
		eA_shiftdict = collections.OrderedDict()
		eE_shiftdict = collections.OrderedDict()

		for i in fl_sizes:
			fl_shiftdict.update( {str(i):'['+str(fl_ribo_shift.loc[i,'riboshift'])+']'})
		for i in eA_sizes:
			eA_shiftdict.update( {str(i):'['+str(eA_ribo_shift.loc[i,'riboshift'])+']'})
		for i in eE_sizes:
			eE_shiftdict.update( {str(i):'['+str(eE_ribo_shift.loc[i,'riboshift'])+']'})

		allLenghts_shiftdict = collections.OrderedDict()
		allLenghts_shiftdict.update(eE_shiftdict)
		allLenghts_shiftdict.update(eA_shiftdict)
		allLenghts_shiftdict.update(fl_shiftdict)

		def build_fl_dense(fldict):

			readsize_dict = fldict
			pop = "fl"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# fl_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			bamfileoutpath = '%s/%s_star_default_allTr' % (fp_assign_path, self.file)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.file, pop)
			if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.file, minlen, maxlen)
			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, self.file, minlen, maxlen)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main_allTR.py \
				--bamfileinput %s \
				--GTFfile %s \
				--twobitfile %s \
				--assignment %s \
				--riboshiftdict %s \
				--threshold %s \
				--totreads %s \
				--outputdata %s \
				--bamfileoutput %s' % (
					rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)
			
			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.file, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.file, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

		def build_eA_dense(eAdict):

			readsize_dict = eAdict
			pop = "eA"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# fl_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			bamfileoutpath = '%s/%s_star_default_allTr' % (fp_assign_path, self.file)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.file, pop)
			if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.file, minlen, maxlen)
			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, self.file, minlen, maxlen)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main_allTR.py \
				--bamfileinput %s \
				--GTFfile %s \
				--twobitfile %s \
				--assignment %s \
				--riboshiftdict %s \
				--threshold %s \
				--totreads %s \
				--outputdata %s \
				--bamfileoutput %s' % (
					rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)

			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.file, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.file, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

		def build_eE_dense(eEdict):

			readsize_dict = eEdict
			pop = "eE"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# fl_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			bamfileoutpath = '%s/%s_star_default_allTr' % (fp_assign_path, self.file)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.file, pop)
			if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)
			
			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, self.file, minlen, maxlen)
			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, self.file, minlen, maxlen)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main_allTR.py \
				--bamfileinput %s \
				--GTFfile %s \
				--twobitfile %s \
				--assignment %s \
				--riboshiftdict %s \
				--threshold %s \
				--totreads %s \
				--outputdata %s \
				--bamfileoutput %s' % (
					rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)

			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.file, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.file, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

		def build_allLengths_dense(alldict):
			
			readsize_dict = alldict
			pop = "aL"
			minlen = str(readsize_dict.keys()[0])
			maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# fl_shift:
			riboshiftdict = str(dict_conv).replace(" ","")
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			
			bamfileoutpath = '%s/%s_star_default_allTr' % (fp_assign_path, self.file)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.file, pop)
			if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

			densityfilepath = '%s/%s/density%sp_%sshift_%s/%s_%sf' %(fp_assign_path, densitystring, assignment, ribosome_shift, pop, self.file, pop)
			densityfileout = '%s/%s_%sf_' % (densityfilepath, self.file, pop)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main_allTR.py \
				--bamfileinput %s \
				--GTFfile %s \
				--twobitfile %s \
				--assignment %s \
				--riboshiftdict %s \
				--threshold %s \
				--totreads %s \
				--outputdata %s \
				--bamfileoutput %s' % (
					rootDir, bamfileinput, self.GTFfile, self.twobitfile, assignment, riboshiftdict, self.threshold, totreads, densityfileout, bamfileoutput)
			print commandstring
			os.system(commandstring)

			## sort and index bamfile
			sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.file, pop)
			os.system(sort_bamfile)
			index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.file, pop)
			os.system(index_bamfile)
			clean_bamfile = 'rm %s' % (bamfileoutput)
			os.system(clean_bamfile)

			### function to combine valid lengths into one density file
			## for example: 31-34, 21-24, and 18-19 with appropriate shifts

		# p1 = multiprocessing.Process(target=build_fl_dense, args=(fl_shiftdict,))
		# p1.start()
		# p2 = multiprocessing.Process(target=build_eA_dense, args=(eA_shiftdict,))
		# p2.start()
		# p3 = multiprocessing.Process(target=build_eE_dense, args=(eE_shiftdict,))
		# p3.start()
		# p4 = multiprocessing.Process(target=build_allLengths_dense, args=(allLenghts_shiftdict,))
		# p4.start()
		# p1.join()
		# p2.join()
		# p3.join()
		# p4.join()
		build_fl_dense(fl_shiftdict)




	def densebuilder_riboshift_custom(self, ribosome_site):
		"""
		Function to map densities with no shift and combine all read lengths
		could be taken out of this class for use in multiprocessing possibly
		Normalization is set to RPM by default
		"""
		print "cust size db", self.customSize

		assignment = "5" # should be five prime assignment to use shift dicts
		ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
		
		bam_inpath = "%s/%s/%s/STAR/%s" % (self.rootpath, self.align_dir, self.genome_name, self.file)
		bamfileinput = "%s/%sAligned.sortedByCoord.out.bam" % (bam_inpath, self.file)
		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')

		fp_assign_path = '%s/FPassignment/%s/%s/%s' % (self.rootpath, self.genome_name, self.folder, self.file)
		bamfileoutpath = '%s/%s_star_default_allTr' % (fp_assign_path, self.file)
		if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)

		densitystring = "Density_allTrDf"
		bamfile = pysam.AlignmentFile(bamfileinput, 'rb')
		totreads = bamfile.count() ## simply count the number of reads in the bamfile to use for normalization

		# totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, self.file)
		# totreadcountf = open(totreads_countfile, "r")
		# totreads = int(totreadcountf.read())
		# totreadcountf.close()

		# if self.norm_type == "raw":
		# 	densitystring = "DensityUnnormalized"
		# 	totreads = 1E6 # this effectively devides readcounts at every position by 1 so that data is not normalized
		# elif self.norm_type == "rpm":
		# 	densitystring = "Density_rpm"
		# 	totreads = bamfile.count() # count the total number of reads in the bamfile and divide by this number (in millions) to get rpm
		# else:
		# 	print "Normalization is not set for densebuilder!"
		# 	sys.exit()

		# Set paths to riboshift dicts for each species of interest
		##

		### define possible ranges for riboshifts
		flLens = range(flmin, flmax+1)
		eALens = range(eAmin, eAmax+1)
		eELens = range(eEmin, eEmax+1)

		print "This is cust size:", self.customSize

		### find which df contains customSize
		if self.customSize in flLens:
			riboshift_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		elif self.customSize in eALens:
			riboshift_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		elif self.customSize in eELens:
			riboshift_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)
		else:
			print "custom size is not in riboshiftdict!"
			sys.exit()


		riboshiftDF = pd.read_csv(riboshift_infile, index_col=0, header=None)
		riboshiftDF.columns = ['riboshift']

		print "custom size == ", self.customSize 

		print riboshiftDF
		print riboshiftDF.index.tolist()
		riboshift_sizes = riboshiftDF.index.tolist()
		# riboshiftdict = OrderedDict()
		riboshiftdict = dict()
		for i in riboshift_sizes:
			if i == self.customSize:
				# riboshiftdict.update( {str(i):'['+str(riboshiftDF.loc[i,'riboshift'])+']'})\
				print str(riboshiftDF.loc[i]['riboshift'])
				riboshiftdict[str(i)] = [str(riboshiftDF.loc[i]['riboshift'])]
		print riboshiftdict



		def build_custom_dense(riboshiftdict):

			readsize_dict = riboshiftdict
			# pop = "%s_shift_%snt" % (self.customSize, readsize_dict[str(self.customSize)])
			pop = "%s_shiftCustom" % (self.customSize)
			print "pop is:", pop 
			# minlen = str(readsize_dict.keys()[0])
			# maxlen = str(readsize_dict.keys()[-1])
			dict_conv = dict(readsize_dict)

			# fl_shift:
			riboshiftdict = str(dict_conv).replace(" ","") ### must remove the space between values for densebuilder input
			print riboshiftdict
			# riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
			bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, self.file)
			bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, self.file, pop)

			densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path, 
																			densitystring, 
																			assignment, ribosome_shift, self.customSize, self.customSize, 
																			self.file, self.customSize, self.customSize)

			densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, 
												self.file, self.customSize, self.customSize)
			if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)

			commandstring = 'python2 %s/riboseq/densebuilder_main_allTR.py \
				--bamfileinput %s \
				--GTFfile %s \
				--twobitfile %s \
				--assignment %s \
				--riboshiftdict %s \
				--threshold %s \
				--totreads %s \
				--outputdata %s \
				--bamfileoutput %s' % (
				
				rootDir,
				bamfileinput, 
				self.GTFfile, 
				self.twobitfile, 
				assignment, 
				riboshiftdict, 
				self.threshold, 
				totreads, 
				densityfileout, 
				bamfileoutput)

			print commandstring
			os.system(commandstring)
			
			# ## sort and index bamfile
			# sort_bamfile = 'samtools sort %s %s/%s_%s_match.sorted' % (bamfileoutput, bamfileoutpath, self.file, pop)
			# os.system(sort_bamfile)
			# index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath, self.file, pop)
			# os.system(index_bamfile)
			# clean_bamfile = 'rm %s' % (bamfileoutput)
			# os.system(clean_bamfile)

		build_custom_dense(riboshiftdict)


def main():
	print samples 
	for folder in samples:
		filenumb = 0
		for file in samples[folder]:
			print "***** starting sample: %s ******" % (samples[folder][filenumb])

			RPexp = RP_custom_readlength(file,folder,rootpath,align_dir,genome_name,customSize,
				 flmin,flmax,eAmin,eAmax,eEmin,eEmax,
				 GTFfile,twobitfile,threshold,UTRfilestring,norm_type)

			start1 = time.time()
			print "RPExp == ", RPexp.customSize
			# RPexp.densebuilder_riboshift_custom("A")
			RPexp.densebuilder_riboshift("A")

			# RPexp.avggene_riboshift_custom_size("1", "A")
			end1 = time.time()
			print(end1-start1, " : old time")

if __name__ == '__main__':
	main()



