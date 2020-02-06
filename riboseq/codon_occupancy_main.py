import sys
import os
import argparse
import struct
import math, csv
import numpy as np

riboseqDir = os.path.dirname(os.path.realpath(__file__)) ## add current directory as rootDir
sys.path.append(riboseqDir)
import rphelper as rph

# Workflow to get E, P, A sites occupancy on ribosomes
# density files are riboshifted already. 
# Exclusion module implemented


### this is developed by Colin Wu, with a few modifications 


class codonoccupancy(object):
	"""
	Instantiate the class

	"""
	def __init__(self,sample_name,motifs,seqwin,UTRfilestring,cds5trim,cds3trim,exclusiondict,outfilepath,outlistfile):
		self.sample_name = sample_name
		self.motifs= motifs
		self.seqwin= seqwin
		self.UTRfilestring= UTRfilestring
		self.cds5trim= cds5trim
		self.cds3trim= cds3trim
		self.exclusiondict= exclusiondict
		self.outfilepath= outfilepath
		self.outlistfile= outlistfile

	def writerows(self, intable, outfilestring):
		fc = open(outfilestring, "w")
		writer = csv.writer(fc)
		writer.writerows(intable)
		fc.close()

	def codonaverage(self):
		outlist,headers,motiffilelist= [],[],[] 
		headers.append("motif") 

		
		for motif in self.motifs: 
								
			motiffile= args.motiffilerootpath+ motif+ "_1.csv" 
			motiffilelist.append(motiffile)   
			headers.append(motif) 
		outlist.append(headers) 

		codon_occu = []
		codon_occu.append(self.sample_name)

		f_output= open(args.outfileparams, "w") 
		f_output.write("Density file is "+str(self.sample_name)+"\n") 
		f_output.write("cds5trim is "+str(args.cds5trim)+"\n")
		f_output.write("cds3trim is "+str(args.cds3trim)+"\n")
		f_output.write("Seqwin is "+str(args.seqwin)+"\n")
		f_output.write("Motiflist is "+str(motiffilelist)+"\n")

		readcountsdict= rph.readcountsf(args.trspdictfilestring) 
		exclusionmodule= exclusionfiles[0] 
		
		if exclusionmodule!= '0':
			exclusiondict= self.readindict(open(exclusionmodule,"rU"))
		else:
			exclusiondict= '0'
		print "Exclusion file is "+str(exclusionmodule)

		
		UTRdict= rph.readindict(open(args.UTRfilestring, "rU"))
		occupancy= self.occupancy(readcountsdict,motiffilelist,exclusiondict,codon_occu,UTRdict,f_output)
		outlist.append(codon_occu)
		f_output.close() 

		co= np.asarray(outlist) # convert outlist to a np.array
		output= co.T 
		# print "output: ", output
		# print "self.outlistfile: ", self.outlistfile
		self.writerows(output,self.outlistfile) # write these rows to a csv


	def occupancy(self,readcountsdict,motiffilelist,exclusiondict,codon_occu,UTRdict,f_output):
		"""
		This is the workhorse of this function
		It is handling the actual density values at the transcript level
		"""
		for genelist in motiffilelist: # motiffilelist is a list contaning all codons to be analyzed
			### genelist is csv file with all transcripts containing a specific codon
			
			f_output.write("Processing "+str(genelist)+"\n")
			f_csv= open(genelist,"rU")
			pausedict= self.readOrderedDict(f_csv) # read in the dictionary with all occurances of a codon in the genome
			count= 0 
			inputpos= 0 
			ratio_sum= 0

			for trsp in pausedict: # iterate through csv file of codon positions - this is pausedict - ex: hg19_AAA_1.csv
				if trsp== "headers":	continue # first line should be: 'headers,gene,chrom,trsp_num,cds_pos'
				else:	trspname= trsp.split("_",1)[0] # example line: 'uc004fop.3,PLCXD1,chrY,0,246' the first line for hg19_AAA_1.csv

				if exclusiondict!= '0':
					if exclusiondict.has_key(trspname):	continue # filter out transcripts present in exclusiondict

				if UTRdict.has_key(trspname):	# check that transcript is present in UTRdict and set utr5len and cdslen
					utr5len= int(UTRdict[trspname][5])
					cdslen= int(UTRdict[trspname][4])
				else:	continue 

				cdspos= int(pausedict[trsp][3])	
				mrnapos= cdspos+ utr5len 
				chrom= pausedict[trsp][1] 
				featnum= int(pausedict[trsp][2]) 
				inputpos+= 1 

				if readcountsdict.has_key(trspname):	cdscounts= readcountsdict[trspname][utr5len:utr5len+ cdslen]	
				else:	continue	

				tooclosetostartorstop= 6	# This is arbitrary, set boundaries around start and stop codon as function of seqwin
				# at a minimum, do not include any transcript within 2 codons

				if cdspos- int(self.seqwin[0])- tooclosetostartorstop < 0 or cdspos+ int(self.seqwin[1])+ tooclosetostartorstop> cdslen:	continue # throw out things too close to start or stop, within 2 codons minimum
				# set the location of density values to be queried
				# with default value of [0,3], window scores three nucleotides within the codon
				else:	loccounts= cdscounts[cdspos- int(self.seqwin[0]): cdspos+ int(self.seqwin[1])] # this is a 3nt window with default settings
					
				trimmed_cdscounts= cdscounts[int(self.cds5trim):cdslen- int(self.cds3trim)] # remove first N nt's and last nt's of cdscounts
				if len(trimmed_cdscounts)== 0 or len(loccounts)== 0:	continue

				cdsdensity= sum(trimmed_cdscounts)/len(trimmed_cdscounts) # density accross entire transcript
				locdensity= sum(loccounts)/len(loccounts) # density at querry location 

				# only throw away transcripts that have no reads in the cds
				if cdsdensity== 0: continue	# include ZERO motif density
				#if cdsdensity==0 or locdensity==0: continue	# exclude ZERO motif density
				
				# calculate codon density ratio relative to that of the transcript
				locratio= locdensity/cdsdensity
				count+= 1 # record a codon occurence as counted
				ratio_sum+= locratio	# sum accross all densities in transcriptome, keep running average, this is additive for all occurances of this codon			
		
			if count== 0: # for when nothing makes it through all forms of filtering
				#print "No genes to average."
				occupancy= "No genes to average"
				f_output.write("No genes to average"+"\n")	
			else:
				occupancy= ratio_sum/count # average density per occrance of this codon
				#print str(count)+" positions averaged. "	
				f_output.write(str(count)+" positions averaged. Occupancy= "+ str(occupancy)+ "\n")	
			codon_occu.append(occupancy) # record occupancy for this codon genome wide 		

	def readOrderedDict(self,f):
		from collections import OrderedDict
		previousgene= ""
		counter= 1
		filegen= csv.reader(f, delimiter=',')
		output= OrderedDict()
		for gene in filegen:
			if gene[0]== previousgene:
				modgenename= gene[0]+"_"+str(counter)
				counter+= 1
			else:
				modgenename= gene[0]
				counter= 1
			output[modgenename]= []
			for element in gene[1:]:
				output[modgenename].append(element)
			previousgene= gene[0]
		return output


	def readindict(self,f):
		previousgene=""
		counter=1
		filegen=csv.reader(f,delimiter=',')
		output = {}
		for gene in filegen:
			if gene[0]==previousgene:
				modgenename=gene[0]+"_"+str(counter)
				counter+=1
			else:
				modgenename=gene[0]
				counter=1
			output[modgenename]=[]
			for column in gene[1:]:
				output[modgenename].append(column)
			previousgene=gene[0]
		return output

if __name__== '__main__':

	parser= argparse.ArgumentParser()
	parser.add_argument('--motiffilerootpath', help= 'motiffilerootpath')
	parser.add_argument('--motifs', help= 'list of motifs')
	parser.add_argument('--trspdictfilestring', help= 'path to density files', required= True)
	parser.add_argument('--sample_name', help= 'name of sample being analyzed')
	parser.add_argument('--UTRfilestring', help= 'UTRfilestring', required= True)
	parser.add_argument('--cds5trim', help= 'trim from 5 end of cds')
	parser.add_argument('--cds3trim', help= 'trim from 3 end of cds')
	parser.add_argument('--seqwin', help= 'seqwin')
	parser.add_argument('--exclusionfiles', help= 'files to exclude transcripts')
	parser.add_argument('--outfileparams', help= 'output path for log file')
	parser.add_argument('--outlistfile', help= 'outputfile path and file name')
	args = parser.parse_args()


	print args.motifs

	import ast
	seqwin= ast.literal_eval(args.seqwin) # convert string to list object 
	# densityfiles= ast.literal_eval(args.densityfiles) # convert to list
	motifs= ast.literal_eval(args.motifs) # convert codons to list
	exclusionfiles= ast.literal_eval(args.exclusionfiles)


	codonavg= codonoccupancy(args.sample_name,motifs,seqwin,args.UTRfilestring,args.cds5trim,args.cds3trim,exclusionfiles,args.outfileparams,args.outlistfile)
	codonavg.codonaverage()

