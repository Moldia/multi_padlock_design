This is a padlock design pipeline for multiplexed assay with multiple probes per target in cDNA-based expression profiling
Mats Nilsson Lab, Stockholm University
Xiaoyan, 2018

How to install and run:
Install blast+, add to path
Install clustalW, add to path
Install Python 3, add to path
Add folder Scripts in Python installation folder to path
Open command window, type: pip install numpy
From ftp.ncbi.nlm.nih.gov download refseq mouse and human mRNA sequences
Download this repository and unzip
Modify config.py file to specify where you have saved mRNA sequence files
Modify config.py to specify maximum number of parallel threads for MSA (multiple sequence alignment) and blast search
Double click probedesign.py to run

The software allows the following input files:
1. A text file with a list of gene acronyms. Only works for humand and mouse transcriptomes. One row per gene.
2. A csv file with gene acronyms and corresponding linker and barcode sequences. One row per gene. No header in the file.
3. A text file with target sequences in fasta format.

Other parameters:
Species: the database against which the sequence specificity should be checked. This will also be the database to look for gene acronyms if input fasta sequences are not provided.
Padlock arm length: the length of each target arm of a padlock probe. The final target sequence will be twice this length.
Tm: lower and upper limit of melting temperatures, in a reaction containig 0.1 uM probe, 0.075 M monovalent salt, 0.01 M bivalent salt and 20% formamide.

Fixed parameters:
For specificity check, in order to say a target sequence is specific enough, it cannot have any sequence in the database with more than 50% of sequence coverage, covering ligation site +- 5 nucleotides, and 80% of homology.

