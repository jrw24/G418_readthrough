# G418_readthrough
All code used in analysis of G418 data is available here.

Clonning this repository and running Wangen_G418_workflow.py will regenerate the figures from Wangen and Green 2019 (https://www.biorxiv.org/content/10.1101/798579v1?rss=1)

*** WARNING ***
Running this workflow requires at least 1.1 Terabytes of free storage space

The pipeline is designed to be run on a computational server with at least 40 threads, and may take an extremely long time to complete if run locally.

To run the analysis pipeline, clone the repository and run the Wangen_G418_workflow.py script
All command line utilities must be downloaded and added to $PATH:


REQUIREMENTS:
  Command Line Utilities, added to $PATH:
  tally: 
  seqtk, 1.0-r31: 
  skewer, 0.2.2:
  STAR, STAR_2.5.3a_modified:
  pigz, 2.3.1:
  samtools, 0.1.19-96b5f2294a:




python 2.7, install using pip:



R 3.4.3:


