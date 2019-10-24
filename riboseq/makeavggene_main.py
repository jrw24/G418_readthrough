import argparse
import rphelper as rph
import sys


class Avggene(object):
	# This is the workflow that creates the average (or "metagene") plots.
	# This function simply averages the counts for every gene.
	# alignpos is a variable: 1 means start codon, 2 means 1 bp past stop codon. 
	# regionlength5 and regionlength3 are the distances up and down of the alignment position that are considered in the average.
	# so average plot is length regionlength5+regionlength3, with alignent position the 1st in regionlength3.

	def __init__(self,regionlength5,regionlength3,trspdict,UTRdict,filterdict,exclusiondict,threshold,alignpos,equalweight,outfilebase):
		self.regionlength5= regionlength5
		self.regionlength3= regionlength3
		self.trspdict= trspdict
		self.UTRdict= UTRdict
		self.filterdict= filterdict
		self.exclusiondict= exclusiondict
		self.threshold= threshold
		self.alignpos= alignpos
		self.equalweight= equalweight
		self.outfilebase= outfilebase

	def totalavg(self):
		missedthresh= 0
		illegalgenes= 0
		genesinlist= 0
		tooshortlist= 0
		averagegene= [0 for num in range(0, int(self.regionlength5)+ int(self.regionlength3))]

		if self.equalweight== '1':	print "Equal weight"	# This is default
		elif self.equalweight== '0':	print "Unequal weight"
		else:	
			print "Err no normalization!!!"
			sys.exit()

		for trsp in self.trspdict:
			if self.UTRdict.has_key(trsp)!= True:	continue
			utr5len= int(self.UTRdict[trsp][5])
			utr3len= int(self.UTRdict[trsp][6])

			# Get proper extensions on gene. Check for cases here where UTR will be too short.
			if self.alignpos== '1':	# Align at start codon.
				if int(self.regionlength5)> utr5len:
					tooshortlist+=1
					continue
			elif self.alignpos== '2':	#align at stop codon.
				if int(self.regionlength3)> utr3len:
					tooshortlist+=1
					continue
			else:	print "No alignment position assigned!"

			exonsplicedcounts= self.trspdict[trsp]

			# Define ORF
			cdsstart= utr5len
			cdsend= len(exonsplicedcounts)- utr3len
			if cdsstart== cdsend:
				print "Error, gene length is 0 for transcript "+ trsp
				sys.exit()
				
			# Filter sequences for those of interest - Here is where the action is. i.e. specific sequences, or go term.
			if self.filterdict!= '0':
				if self.filterdict.has_key(trsp):
					print "Gene included: "+trsp
				else:	continue

			if self.exclusiondict!= '0':
				if self.exclusiondict.has_key(trsp):
					print "Gene excluded: "+trsp
					continue

			cdscounts= exonsplicedcounts[cdsstart:cdsend]
			cdsdensity= sum(cdscounts)/len(cdscounts)
			if cdsdensity*float(1000)< int(self.threshold):	# Threshold on cds density: (thresholding on "rpkm")
				missedthresh+= 1
				continue
		
			#Cut down to the size we want and make sure we don't run off the rails.
			if self.alignpos== '1':
				# Check for enough length on 3' end of START codons (5' end already checked above):
				if len(cdscounts)< int(self.regionlength3):
					tooshortlist+= 1
					continue

				countlist= exonsplicedcounts[cdsstart- int(self.regionlength5): (cdsstart+ int(self.regionlength3))]
				genesinlist+= 1

				if self.equalweight!= '1':	cdsdensity=1
				# Add to growing average.
				if cdsdensity!= 0:	
					for i in range(len(countlist)):	
						averagegene[i]+= countlist[i]/cdsdensity # Use CDS density for normalziation

			elif self.alignpos== '2':
				# Check for enough length on 5' end (3' end already checked above):
				if len(cdscounts)< int(self.regionlength5):
					tooshortlist+= 1
					continue

				countlist= exonsplicedcounts[cdsend- int(self.regionlength5): (cdsend+ int(self.regionlength3))]
				genesinlist+= 1

				if self.equalweight!= '1':	cdsdensity= 1
				# Add to growing average.
				if cdsdensity!= 0:	
					for i in range(len(countlist)):	averagegene[i]+= countlist[i]/cdsdensity # Use CDS density for normalziation
			else:	print "No alignment position assigned!"


		if self.equalweight== '1':
			for m in range(len(averagegene)):
				if genesinlist!= 0: averagegene[m]/= genesinlist
				#else: print "Error, no genes to average."	

		print "Genes under threshold: "+str(missedthresh)
		print "Genes removed for overlap or no UTR or dubious or non-major chrom or other: "+str(illegalgenes)
		print "Genes removed because feature(s) too short: "+str(tooshortlist)
		print "Genes in average: "+str(genesinlist)

		fc= open(self.outfilebase+"_"+str(self.alignpos)+"_output.txt","a")
		fc.write("\n")
		fc.write("Genes removed for overlap or no UTR or dubious or non-major chrom or other: "+str(illegalgenes)+"\n")
		fc.write("Genes removed because feature(s) too short: "+str(tooshortlist)+"\n")
		fc.write("Genes in average:"+str(genesinlist)+"\n")
		fc.close()
		
		# Write out avg to csv
		csvout= rph.avgcsvout(averagegene,self.outfilebase+"_avg_"+str(self.alignpos),[int(self.regionlength5),int(self.regionlength3)])
		

if __name__== '__main__':
	parser= argparse.ArgumentParser()
	parser.add_argument('--regionlength5', help= 'length upstream of START or STOP')
	parser.add_argument('--regionlength3', help= 'length downstream of START or STOP')
	parser.add_argument('--trspdictfilestring', help= 'input transcript density files')
	parser.add_argument('--UTRfilestring', help= 'UTRs file')
	parser.add_argument('--filtermodule', default= '0', help= 'motifs to include')
	parser.add_argument('--exclusionmodule', default= '0', help= 'motifs to exclude')
	parser.add_argument('--threshold', default= '0', help= 'thresholding on rpkm')
	parser.add_argument('--alignpos', default= '1', help= '1: START, 2: STOP')
	parser.add_argument('--equalweight', default= '1', help= 'equal or unequal weighting')
	parser.add_argument('--outfilebase', help= 'output file name')
	parser.add_argument
	args = parser.parse_args()

	# Write output file of comments.
	comments= "This is a test run of Avggene.\nIt makes an average of all genes.\n"
	comments+= "Threshold signifies minimal rpkm needed in coding region for gene to be in the average.\n"
	comments+= "alignpos =1 anchors average around the start codon and only includes 5'UTRs. alignpos =2 is the same for stop codon."
	fc= open(args.outfilebase+"_"+str(args.alignpos)+"_output.txt","w")
	fc.write(comments)
	fc.write("\n")
	fc.write("Avggene was called with parameters:\n")
	fc.write("transcripts= "+str(args.trspdictfilestring)+"\n")
	fc.write("filtermodule= "+str(args.filtermodule)+"\n")
	fc.write("exclusionmodule= "+str(args.exclusionmodule)+"\n")
	fc.write("threshold= "+str(args.threshold)+"\n")
	fc.write("regionlength5= "+str(args.regionlength5)+"\n")
	fc.write("regionlength3= "+str(args.regionlength3)+"\n")
	fc.write("equalweight= "+str(args.equalweight)+"\n")
	fc.write("alignpos= "+str(args.alignpos)+"\n")
	fc.close()

	if args.filtermodule!= '0':  filterdict= rph.readindict(open(args.filtermodule,"rU"))
	else:   filterdict= '0'
	if args.exclusionmodule!= '0':  exclusiondict= rph.readindict(open(args.exclusionmodule,"rU"))
	else:   exclusiondict= '0' 
	trspdict= rph.readcountsf(args.trspdictfilestring)
	UTRdict= rph.readindict(open(args.UTRfilestring, "rU"))
	metagene= Avggene(args.regionlength5,args.regionlength3,trspdict,UTRdict,filterdict,exclusiondict,args.threshold,args.alignpos,args.equalweight,args.outfilebase)
	metagene.totalavg()


