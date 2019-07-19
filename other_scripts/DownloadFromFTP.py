from ftplib import FTP
import os

# connect to NCBI FTP server
ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()     # anonymous login

# move to refseq directory
ftp.cwd('refseq')

# list contents
ftp.retrlines('LIST')

# move to species specific subdirectory
ftp.cwd('H_sapiens/mRNA_Prot')
allfiles = ftp.nlst()

for file in allfiles:
    if 'rna.fna' in file:
        ftp.retrbinary('RETR ' + file,
                       open(os.path.join(r'F:\ProbeDesign\RefSeqDatabase', file), 'wb').write)
        print(file)

ftp.quit()