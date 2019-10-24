__author__ = 'Colin Wu'

"""
Helper functions for rib-pro
"""
import GFF
import urllib
import csv

def readcountsf(filestring):
    """ Function to read in counts files. Use for floats
    From bd.readcountsf(filestring)
    Called readcountsf(countsfilestring+"_plus_")
    Alternative: Read from wig instead? This would mean that the individual
                 binary files don't have to be created anymore
    returns dictionary resultlist[chrom] = []
    should be pretty easy to convert to resultlist
    """
    import struct
    keys=[]
    resultlist={}
    f2=open(filestring+"keys","r")
    for line in f2:
        keys.append(line.rstrip('\n'))  # rstrip takes off the end of line character.
    for chrom in keys:
        resultlist[chrom]=[]
        with open(filestring+chrom,"rb") as f:
            nextval = f.read(4) 
            while nextval != "": 
                resultlist[chrom].append(struct.unpack("f",nextval)[0])    
                # This returns a tuple, not an int, so need to take 1st item.
                nextval = f.read(4)
    f2.close()
    # Note that there is no file closing, probably not ideal, but everything seems to 
    # work fine.
    return resultlist
 
def makeGFFlist(GFFname):
    """ Tool for loading the entire yeast genome into memory
    From seqtools
    Called st.makeGFFlist(GFF.parse(codingGFF))
    Returns dictionary GFFlist[chr.id] = chr for chr in GFFgen
    Called for main coding GFF but not utr5GFF and utr3GFF -- generalize?
    Will this be affected it GFF is changed? Probably not, no parsing here, just storing
    """
    GFFlist={}
    for chr in GFF.parse(GFFname):
        GFFlist[chr.id]=chr
    return GFFlist 

def makeutrtable(utrgffgen):
    """ Table generator for next function.
        From genometools
        Called rphelper.makeutrtable(GFF.parse(utrGFF))
        Returns dictionary table[feature.id[:-5]] = [start, end]
        utrgffgen is the gffgen for 3'UTR from Nagalkshmi for R64.
        Called for utr5GFF and utr3GFF but not main coding GFF -- generalize?
        Will this be affected if GFF is changed? Probably not, no parsing here, just storing
    """
    table={}
    for chr in utrgffgen:
        for feature in chr.features:
            table[feature.id[:-5]]=[feature.location.start.position,feature.location.end.position]
    return table

def clean_feature_dict(featured):
    # delete 'chrMito' and '2-micron'
    # delete features whose type is not gene    
    # so that all neighbors are non-dubious
    for chrom in featured.keys():
        if chrom == 'chrMito' or chrom == '2-micron':
            del featured[chrom]
            continue
        temp = []
        for sub_feature in featured[chrom].features:
            if sub_feature.type == 'gene':
                if sub_feature.qualifiers.has_key("orf_classification"):
                    if sub_feature.qualifiers["orf_classification"][0] != "Dubious":
                        temp.append(sub_feature)
        featured[chrom].features = temp
    return featured

def parse_GFF(utrGFF):
    """ Very thin wrapper that tries to do makeutrtable(GFF.parse(utrGFF))
        but checks for IOError and returns an empty dictionary instead.
    """
    try:
        return makeutrtable(GFF.parse(utrGFF))
    except IOError:
        print "Warning! " + utrGFF + " couldn't be found."
        if raw_input("'c' to continue with empty dictionary\n") == 'c':
            return {}
        else:
            quit()

def qualifier_to_string(feature, qualkey):
    """ Another thin wrapper, this time to grab a tostring version of 
        whatever qualifier I am looking for.
    """
    try:
        return urllib.quote(";".join(feature.qualifiers[qualkey]))
    except KeyError:
        return 'NA'

def get_gene_list(csv):
    """ Parses a csv to get a list of genes.
    """
    pass

def check_GFF(GFFs):
    """ If GFFs is already a list, return (GFFlist, utrtable5, utrtable3)
        Otherwise return (GFFlist, {}, {})
        Maintains backward comptability
    """
    if type(GFFs)==list:        
        return (GFFs[0], GFFs[1], GFFs[2])
    else:
        print "Warning: UTRs of neighboring features not considered in overlap check."
        return (GFFs, {}, {})

def check_shift(bp):
    """ Checks the shift bp parameter
        If bp is an int, then return [bp, bp, 0]
        If bp is a list of length 2 then return [bp[0], bp[1], 0]
        Else return [bp[0], bp[1], bp[2]]
        Where the list of shift is [5' upstream, 3' downstream, ribosome shift]
        If bp provided is negative, raiseValueError
    """
    if type(bp) != int:
        if len(bp)==2:
            bp.append(0)
    else:
        bp = [bp, bp, 0]
    if bp[0] < 0 or bp[1] < 0:
        print "Error, bp is negative!"
        raise ValueError 
    return bp


def wigtocounts(wigfile,chrlen=None):
    """ Function to convert wigfiles back to a countsfile.
        This is pretty basic and handles 2 types of wigs (fixed and variable) 
        that I have created in the past. 
        Fixedstep assumed to be stepsize 1 and start value to be 1. Can be changed later 
        by adding more code.
        chrlen is a dictionary of chromosomes that is provided for variablestep wigs. 
        It specifies the length of each chrom.
    """
    resultlist={}
    f=open(wigfile)
    chrom=""
    steptype=-1         # Fixed is 1, Variable is 0, and not assigned is -1.
    for line in f:      # Should be okay to be checking past ends of lines. Sloppy but no problems yet.
        if line[0:9]=='fixedStep':  
            steptype=1
            if line[11:17]=="chrom=":
                chrom=parsenext(line[17:],'  ')
                chrom=cleanchrom(chrom)  # Remove any end of lines or quotations.
                resultlist[chrom]=[]
            else:
                print "Error - no chrom."
                return -1
            
        elif line[0:12]=='variableStep':
            assert chrlen
            steptype=0
            if line[14:20]=="chrom=":
                chrom=gentools.parsenext(line[20:],'  ')
                chrom=cleanchrom(chrom) # Remove any end of lines or quotations.
                resultlist[chrom]=[float(0) for x in range(chrlen)]
            else:
                print "Error - no chrom."
                return -1
            
        else:
            if steptype==1:
                resultlist[chrom].append(float(line))
            if steptype==0:
                
                resultlist[chrom][int(parsenext(line,' '))-1]=float(gentools.parselast(line,'  '))
            else:
                continue    # This happens on early comment lines.
                #print "Step type not specified. Error."
                #return -1

    f.close()
    return resultlist

def parsenext(string,delimiter):
    """ Returns all text from the start of the input string up to the 
    occurence of the delimiter."""
    delimlen=len(delimiter)
    i=0
    strlen=len(string)
    while(i<strlen):
        if string[i:i+delimlen]==delimiter:
            break
        i+=1
        
    return string[0:i]

def cleanchrom(chrom):
    """Support function for converting wig to countslist function above."""
    if chrom[-1]=="\n":
        chrom=chrom.strip()
    if chrom[0]=='\"':
        chrom=chrom[1:-1]
    return chrom


# Make a csv back into a dictionary. f is open file handle. Handles multiple hits for same gene by putting in a _# term for them.
def readindict(f):
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


# Function to write out csv files from makeposavg or avggene
def avgcsvout(inlist,outfilestring,seqwin):
    t=[]
    headers= ["position","avg"]
    t.append(headers)

    #for i in range (-seqwin[0]+1,seqwin[1]): # position starts at -seqwin[0]+1
        #newline= [i, inlist[i+seqwin[0]]]

    for i in range (seqwin[1]+seqwin[0]): # position starts at 0
        newline= [i, inlist[i]]
        t.append(newline)

    fa = open(outfilestring+".csv", "w")
    writer = csv.writer(fa)
    writer.writerows(t)
    fa.close()


