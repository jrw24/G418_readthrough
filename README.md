# G418_readthrough
All code used in analysis of G418 data is available here.

Exact code for generating figures is found in G418_readthrough/figures/figscripts, 
Example figures are found in G418_readthrough/figures/output_figures

Clonning this repository and running Wangen_G418_workflow.py will regenerate the figures from Wangen and Green. 2019. (https://www.biorxiv.org/content/10.1101/798579v1?rss=1)

**Change threadNumb to the desired number of threads, default is 40**


### *** WARNING ***
**Running this workflow requires at least 1.1 Terabytes of free storage space to generate all necessary files**

**The pipeline is designed to be run on a computational server with at least 40 threads, and may take an extremely long time to complete if run locally.**

### Overview of Wangen_G418_workflow.py
1. Generate Genomes
   - Download hg38 and Gencode annotations
   - Parse the GTF file and choose single isoform for each gene
   - Build all other necessary annotation files from parsed GTF
   - Parse the GTF file and create a second GTF file that contains all valid transcripts
   - Build non codon RNA depletion annotation files
   - Build indexes for STAR alignments
2. Raw Sequencing Data
   - Download FASTQ files from SRA ***Pending release of data***
   - create merged FASTQ files of replicates for select ribosome profiling experiments
3. Ribosome Profiling data processing
   - run main analysis pipeline on ribosome profiling data, building all files required for generation of figures
   - mapping to all possible transcripts is only performed for select datasets, as this takes a long time
4. RNAseq data processing
   - run main analysis pipeline on RNAseq data, building all files required for generation of figures
5. Plot Figures
   - Run individual scripts that generate all figures in the manuscript
   - Figures can be compared to examples in output_figures to validate successful completion of analysis



### Description of subfolders: 
- **Data** contains raw luciferase measurments. All raw and processed data files will be saved here
- **figures** contains scripts for generating all figures in *figures/figscripts* and example output figures in *figures/output_figures*
- **genomes** contains currated refseq annotations of rRNA in fasta format. All annotation files will be created in this directory
- **GFF** contains a python module that must be coppied to python2.7/site-packages for import into python. Can also be found here: https://github.com/chapmanb/bcbb/tree/master/gff/BCBio/GFF
- **riboseq** contains scripts for processing ribosome profiling data
- **RNAseq** contains scripts for processing RNAseq data
- **utils** contains scripts for generating annotation files


To run the analysis pipeline, clone the repository and run the **Wangen_G418_workflow.py script**

All command line utilities must be downloaded and added to $PATH


### REQUIREMENTS:
#### Command Line Utilities, added to $PATH:
- tally, 15-065: (https://www.ebi.ac.uk/research/enright/software/kraken) 
- seqtk, 1.0-r31: (https://github.com/lh3/seqtk)
- skewer, 0.2.2: (https://github.com/relipmoc/skewer)
- STAR, STAR_2.5.3a_modified: (https://github.com/alexdobin/STAR)
- pigz, 2.3.1: (https://zlib.net/pigz/)
- samtools, 0.1.19-96b5f2294a:(https://github.com/samtools/samtools)
- kpLogo (http://kplogo.wi.mit.edu/)
- kentUtils (https://github.com/ENCODE-DCC/kentUtils)


#### python 2.7, install using pip:

*built in modules do not have version number listed*
- os 
- sys 
- ftplib
- subprocess
- glob 
- struct 
- ast 
- time 
- datetime 
- collections 
- importlib
- math 
- argparse *1.1*
- pandas *0.22.0*
- numpy *1.14.0*
- pysam *0.13*
- scipy *1.0.1*
- statsmodels *0.9.0*
- twobitreader *3.1.5*
- Bio *1.58*
- multiprocessing *0.70a1*
- pathos *0.2.1*
- matplotlib *2.2.2*
- seaborn *0.9.0*
- csv *1.0*
- urllib *1.17*
- lifelines *0.19.5*
- logomaker *0.8*
- scikit-learn *0.20.3*
- GFF ***copy and paste GFF folder into ~/lib/python2.7/site-packages/***


#### R 3.4.3:
- ggplot2 *2.2.1*
- plyr *1.8.4*
- reshape2 *1.4.3*
- scales *0.5.0*
- xtail *1.1.5*
- glue *1.3.1*



