# for multi_padlock_design
# Xiaoyan, 2018

import os
from lib import readfastafile
from lib import retrieveseq
from lib import createoutput


def correctinput(string):
    """ Change from backslash to slash """
    string = string.encode('unicode-escape').decode()   # un-escape escape characters
    string = string.replace('\\', '/')  # change every unescaped backslash to slash
    string = string.replace('//', '/')  # change original double backslash (one for escaping) to one slash
    return string


def checkspecies(species):
    """ Only human and mouse are currently supported """
    success = False
    if not (species in ["mouse", "human"]):
        print("Could not identify species. Try again.")
    else:
        success = True
    return success


def readgenefile(genefile):
    """ Get the gene list """
    success = False
    try:
        with open(str(genefile), 'r') as f:
            lines = [line.rstrip('\n').split(',') for line in f]
            genes = [line[0] for line in lines]
            linkers = [line[1:] for line in lines]
            success = True
    except IOError:
        print("Gene list file not found. Try again.")
        genes = []
        linkers = []
    return success, genes, linkers


def readseqfile():
    """ Import sequences when gene list is not available """
    success_f = False
    while not success_f:  # get input sequences if it is given instead of gene acronyms
        seqfile = input("File containing target sequences in FASTA format: ")
        seqfile = correctinput(seqfile)
        success_f, headers, sequences = readfastafile.readfasta(seqfile)
        headers_wpos = []
        basepos = []
        for c, header in enumerate(headers):
            headers_wpos.append(header + ', 1 to ' + str(len(sequences[c])))
            basepos.append([0, len(sequences[c])])

    return success_f, seqfile, headers, sequences, headers_wpos, basepos


def makeoutputdir(outdir):
    """ Create output directory """
    success = False
    try:
        os.mkdir(outdir)
        success = True
    except WindowsError as e:
        if 'Error 183' in str(e):    # already existing
            print("Directory already existing.")
            success = True
        elif 'Error 3' in str(e):    # directory does not exist
            try:
                os.makedirs(outdir)     # create subdirectories
                success = True
            except:     # catch all exceptions
                print("Could not create output directory. Try again.")
    return success


def armlength(armlen):
    """ Length of a padlock probe arm """
    success = False
    if not armlen > 6:
        print("Padlock arm length too short. Should be at least 7. Try again")
    else:
        success = True
    return success


def spacing(interval):
    """ Minimum distance between two targets """
    success = False
    if not interval >= 0:
        print("Spacing must be at least 0. Try again.")
    else:
        success = True
    return success


def tmthreshold(t1, extreme):
    """ Lower threshold for Tm screening """
    success = False
    if t1 < extreme:
        print("Threshold must be higher than " + str(extreme) + ". Try again.")
    else:
        success = True
    return success


def nprobes(n):
    """ Fixed number of probes per gene """
    success = False
    if not n >= 1:
        print("Fixed number of probes per gene must be at least 1. Try again.")
    else:
        success = True
    return success


def getdesigninput():
    """ Check all keyboard inputs and format target sequences """
    success_s = False  # species
    success_g = False  # gene acronyms
    success_f = False  # fasta file
    success_d = False  # output directory
    success_a = False  # arm length
    success_i = False  # interval between targets
    success_t = False  # Tm threshold
    success_n = False  # fixed number of output sequences

    # loop until all the keyboard inputs are correct
    while not success_s:
        species = input("Specify the species (human or mouse): ").lower()
        success_s = checkspecies(species)

# if species in (["human", "mouse"]):
    # when human or mouse, possible to load only gene list
    while not success_g:
        genefile = input("File containing gene acronyms (text/csv file with one gene per row, or with linker sequences).\n"
                             "Just press Enter if input sequences can be provided:\n")
        if len(genefile):
            genefile = correctinput(genefile)
            success_g, allgenes, alllinkers = readgenefile(genefile)
            # remove gene name duplicates
            genes = []
            linkers = []
            for c, name in enumerate(allgenes):
                if name not in genes:
                    genes.append(name)
                    linkers.append(alllinkers[c])
            # genes = list(set(genes))
        else:
            success_g = True
            success_f, seqfile, headers, sequences, headers_wpos, basepos = readseqfile()
            toavoid = [':', '/', '\\', '[', ']', '?', '"', ' ', '<', '>']
            genes = []
            linkers = []
            for header in headers:
                for i in toavoid:
                    header = header.replace(i, '')
                genes.append(header)
                linkers.append([])
# else:
#     success_f, seqfile, headers, headers_wpos, sequences = readseqfile()
#     basepos = []
#     linkers = []

    while not success_d:
        outdir = input("Output directory: ")
        outdir = correctinput(outdir)
        if outdir.find(' ') > -1:
            print("No whitespace allowed in the directory path!")
        else:
            success_d = makeoutputdir(outdir)

    # create temporary output folder
    outdir_temp = createoutput.tempdir(outdir)
    t = outdir_temp.split('TempFolder')[1]
    os.mkdir(outdir_temp)

    while not success_a:
        armlen = input("Length of one padlock arm (nt): ")
        success_a = armlength(int(armlen))

    while not success_i:
        interval = input("The minimum number of nucleotides between targets: ")
        success_i = spacing(int(interval))

    while not success_t:
        t1 = input("The lower threshold for target Tm\n"
                       "(0.1 uM oligo conc., 0.075 M monovalent salt, 0.01 M bivalent salt, 20% formamide): ")
        temp = tmthreshold(int(t1), 0)
        if temp:
            while not success_t:
                t2 = input("The upper threshold for target Tm\n"
                               "(0.1 uM oligo conc., 0.075 M monovalent salt, 0.01 M bivalent salt, 20% formamide): ")
                success_t = tmthreshold(int(t2), int(t1))

    while not success_n:
        n = input("Number of probes per gene (skip by pressing Enter): ")
        if len(n):
            success_n = nprobes(int(n))
        else:
            success_n = True


    # find hits if no target sequence is given
    if not success_f:
        # find genes in the database
        hits = retrieveseq.querygenes(genes, species)

        # if any gene is not found in the RefSeq acronym list, write to file
        nohit = [i for i in range(len(genes)) if len(hits[i]) == 0]
        if len(nohit):
            with open(os.path.join(outdir, '0.AcronymNotFound_' + t + '.txt'), 'w') as f:
                for i in nohit[::-1]:
                    f.write("%s\n" % genes[i])
                    del genes[i]  # remove genes that are not found
                    del linkers[i]
                    del hits[i]

        # find sequences (MSA included if multiple variants)
        headers, basepos, sequences, msa, nocommon, variants = retrieveseq.findseq(genes, hits, outdir_temp)

        # genes that have no common sequence among all variants
        nocommon = [c for c, i in enumerate(genes) if i in nocommon]
        if len(nocommon):
            with open(os.path.join(outdir, '0.NoConsensusSequence_' + t + '.txt'), 'w') as f:
                for i in nocommon[::-1]:
                    f.write("%s\n" % genes[i])
                    del genes[i]  # remove genes that are not found
                    del linkers[i]
                    del hits[i]

        idxmsa = [c for c, i in enumerate(genes) if i in msa]
        geneorder = [c for c, i in enumerate(genes) if i not in msa]
        [geneorder.append(i) for i in idxmsa]
        linkers = [linkers[i] for i in geneorder]
        genes = [genes[i] for i in geneorder]

        # write found sequences to output file
        with open(os.path.join(outdir, '1.InputSeq_' + t + '.fasta'), 'w') as f:
            for i, base in enumerate(basepos):
                if isinstance(base[0], int):  # no variants found
                    f.write("%s, %d to %d\n" % (headers[i], base[0] + 1, base[1] + 1))
                    f.write("%s\n\n" % sequences[i])
                else:  # more than one variants found
                    for j, subbase in enumerate(base):
                        f.write("%s, %d to %d\n" % (headers[i], subbase[0] + 1, subbase[1] + 1))
                        f.write("%s\n\n" % sequences[i][j])
        temp = readfastafile.readfasta(os.path.join(outdir, '1.InputSeq_' + t + '.fasta'))
        headers_wpos = temp[1]
        sequences = temp[2]

    # write an overview log file
    with open(os.path.join(outdir, 'log_' + outdir_temp.split('TempFolder')[1] + '.txt'), 'w') as f:
        f.write("Multiple padlock design v1\nDate and time: %s\n" % outdir_temp[-14:])
        f.write("Species: %s\n" % species)
        f.write("Output directory: %s\n" % outdir)
        f.write("Padlock arm length: %s nt\n" % armlen)
        f.write("Space between targets: %s nt\n" % interval)
        f.write("Target Tm range: %d to %d after adjustment\n" % (int(t1), int(t2)))
        f.write("Number of probes per gene: %s\n" % n)
        if not success_f:
            f.write("Input file: %s\n" % genefile)
            f.write("Unique gene acronyms found: %d\n" % (len(genes) + len(nocommon)))
            f.write("Number of genes not found in the corresponding database: %d\n" % len(nohit))
            f.write("Number of genes with multiple sequence variants: %d\n"
                    % len(msa))
            f.write("Number of genes with no common sequence after multiple sequence alignment: %d\n"
                    % len(nocommon))
        else:
            f.write("Input file: %s\n" % seqfile)
            f.write("Number of sequences processed from the input file: %d\n" % len(headers))
    return (species, int(armlen), int(interval), int(t1), int(t2), n),\
           (outdir, outdir_temp), \
           (genes, linkers, headers, variants),\
           (basepos, headers_wpos, sequences)


def checkformat(headers):
    fmt = 'NCBI'
    c = 0
    while fmt == 'NCBI' and c < len(headers):
        if headers[c][:3] != '>gi' and headers[c][3] != '_':
            fmt = 'Unknown'
        else:
            c += 1
    return fmt
