This is used to design multiple padlock probes per target for cDNA-based expression profiling.
Xiaoyan, 2017

How to run:
Install blast+, add to path
Install clustalW, add to path
Install Python, add to path
Add folder Scripts in Python installation folder to path
Open command window, type: pip install numpy
From ftp.ncbi.nlm.nih.gov download refseq mouse and human mRNA sequences and decompress to a desired location
Modify config.py file to specify where you have saved mRNA fasta files
Modify config.py to specify maximum number of parallel threads for MSA (multiple sequence alignment) and blast search
Double click probedesign.py to run
