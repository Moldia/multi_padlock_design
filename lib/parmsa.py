# run parallel multiple sequence alignment (ClustalW) to find common sequences of transcript variants
# for multi_padlock_design
# Xiaoyan, 2016-8-4


import subprocess
import config


Processes = []
NextProcess = 0


def runmsa(msafile):
    """ Start a new incident of ClustalW """
    global Processes
    global NextProcess

    msa_process = subprocess.Popen("clustalw2 -quicktree -quiet -infile=" + msafile)
    NextProcess += 1
    Processes.append(msa_process)


def readmsa(alnfile):
    """ Read alignment file """
    nvars = 0     # number of transcript variants
    alnresult = []
    # read file
    with open(alnfile, 'r') as f:
        for row, line in enumerate(f):
            if line[:7] == 'CLUSTAL':   # skip the first row
                pass
            elif line == '\n':
                if nvars == 0 and row > 2:
                    nvars = row
                    nvars -= 4  # first three rows are always dummy rows, and the last row before empty new line is the alignment line
            else:
                alnresult.append(line)

    # get the gi/MN # of the first sequence in alignment
    colexclude = alnresult[0].split(' ').count('') + len(alnresult[0].split(' ')[0]) + 1    # columns to exclude
    name = alnresult[0][:colexclude]
    if '|' in name:
        name = name.split('|')
        name = name[1]      # take only NM or gi number
    else:
        name = name.replace(' ', '')
        # name = name.split(' ')
        # name = name[1][1:]

    alnresult = [line[colexclude:].rstrip('\n') for line in alnresult]

    # combine multiple lines of alignment results from the first sequence
    alnline = str()
    sequence = str()
    for row, line in enumerate(alnresult[0:len(alnresult)-1:nvars+1]):
        sequence += line
        alnline += alnresult[(row+1)*(nvars+1)-1]

    # define starting and end positions of aligned regions, ignoring matches shorter than 20bp
    aligned = []
    alnindex_start = -1
    alnindex_end = -1
    for c, i in enumerate(alnline):
        if i == '*':
            if alnindex_start < 0:
                alnindex_start = c
            else:
                alnindex_end = c
        else:
            if alnindex_end > 0 and (alnindex_end-alnindex_start) >= 20:    # only no shorter than 20 bp
                aligned.append([alnindex_start, alnindex_end])
            alnindex_start = -1
            alnindex_end = -1

    # the last item
    if alnindex_end > 0:
        aligned.append([alnindex_start, alnindex_end])
    # remove the last item if the same as the second last one
    if len(aligned) > 1 and aligned[-1] == aligned[-2]:
        del aligned[-1]

    # convert to original base index (based on the first sequence in alignment) and get sequences
    baseindex = []
    alnseq = []
    for i in aligned:
        alnseq.append(sequence[i[0]:i[1]+1])
        base = sequence[:i[0]].count('A') + sequence[:i[0]].count('C') + \
               sequence[:i[0]].count('G') + sequence[:i[0]].count('T')      # skip deletions
        baseindex.append([base, base+i[1]-i[0]])

    return (name, baseindex, alnseq)


def runningmsa(dirname, msa):
    """ Check processes and distribute MSA jobs """
    global Processes
    global NextProcess

    for p in range(len(Processes)-1,-1,-1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
    while len(Processes) < config.max_MSA_processes and NextProcess < len(msa):  # more to do and some spare slots
        msafile = dirname + '/' + msa[NextProcess] + '_variants.fasta'
        runmsa(msafile)


def continuemsa(dirname, msa):
    """ Continue MSA until all done """
    global Processes
    global NextProcess

    runningmsa(dirname, msa)
    while len(Processes) > 0:
        runningmsa(dirname, msa)

    Names = []
    BasePos = []
    Seqs = []
    if NextProcess == len(msa) and len(Processes) == 0:
        for aln in msa:
            alnfile = dirname + '/' + aln + '_variants.aln'
            tempout = readmsa(alnfile)
            Names.append(tempout[0])
            BasePos.append(tempout[1])
            Seqs.append(tempout[2])
    return (Names, BasePos, Seqs)


