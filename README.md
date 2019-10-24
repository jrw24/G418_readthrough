# G418_readthrough
All code used in analysis of G418 data is available here.

Exact code for generating figures is found in G418_readthrough/figures/figscripts

Clonning this repository and running Wangen_G418_workflow.py will regenerate the figures from Wangen and Green. 2019. (https://www.biorxiv.org/content/10.1101/798579v1?rss=1)

### *** WARNING ***
**Running this workflow requires at least 1.1 Terabytes of free storage space to generate all necessary files**

**The pipeline is designed to be run on a computational server with at least 40 threads, and may take an extremely long time to complete if run locally.**

To run the analysis pipeline, clone the repository and run the **Wangen_G418_workflow.py script**

All command line utilities must be downloaded and added to $PATH


### REQUIREMENTS:
#### Command Line Utilities, added to $PATH:
- tally: 
- seqtk, 1.0-r31: 
- skewer, 0.2.2:
- STAR, STAR_2.5.3a_modified:
- pigz, 2.3.1:
- samtools, 0.1.19-96b5f2294a:
- kpLogo
- kentUtils


#### python 2.7, install using pip:

* *built in modules do not have version number listed* *
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



