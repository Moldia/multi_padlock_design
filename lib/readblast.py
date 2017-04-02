# for multi_padlock_design
# Xiaoyan, 2017

import os
import createoutput
import numpy as np

notmapped = []
funmap = []


def readblastout(file, armlength):
    """ Read the results from blast """
    global funmap
    global notmapped
    specific = True
    mappable = False
    with open(file, 'r') as fh:
        for line in fh:
            if specific:
                scores = []
                columns = line.split('|')
                if len(columns) <= 3:
                    if columns[1][:2] == 'NM' or columns[1][:2] == 'NR':    # check if the hit is transcript (NM) or non-coding RNA (NR), remove all predicted (XM)
                        scores = columns[-1].split(',')
                else:
                    if columns[3][:2] == 'NM' or columns[3][:2] == 'NR':
                        scores = columns[-1].split(',')
                if len(scores):
                    if 2*armlength*.5 < int(scores[2]) < 2*armlength and float(scores[1]) > 80 and int(scores[5]) < armlength-4 and int(scores[6] > armlength+5):
                        # more than 50% coverage, 80% homology, and non-target sequence covers ligation site +- 5
                        specific = False

                    if not mappable:
                        if float(scores[1]) == 100 and int(scores[2]) == 2*armlength:
                            mappable = True
        if not mappable:
            with open(file[0:-10] + '.fasta', 'r') as f:
                seq = f.readlines()
                funmap.write('Could not map sequence in ' + file[:-10] + '!\n')
                funmap.write(seq[1])
                notmapped.append(int(file[:-10].split('_')[-1])-1)
    return specific


def getcandidates(listSiteChopped, headers, dirnames, armlength):
    """ Get specific fragments """
    global notmapped
    global funmap
    siteCandidates = []
    notMapped = []
    t = dirnames[1].split('TempFolder')[1]
    with open(os.path.join(dirnames[0], '2.Unmappable_' + t + '.txt'), 'w') as funmap:
        for i, sites in enumerate(listSiteChopped):
            funmap.write("%s\n" % headers[i])

            fname = createoutput.blastinfilename(dirnames[1], headers[i])

            notmapped = []
            blast_bw = []
            for j, target in enumerate(sites):
                fblast = fname + '_' + str(j + 1) + '_blast.txt'
                blast_bw.append(readblastout(fblast, armlength))

            # find sequences that are specific enough
            idxspecific = np.nonzero(blast_bw)[0]
            tempCandidates = np.array(sites)
            sitespecific = tempCandidates[idxspecific]

            # write unmappable sites
            notmapped = tempCandidates[notmapped]
            recorded = False
            for j, nomap in enumerate(notmapped):
                if j == 0:
                    funmap.write("\nUnmapped sequence starting position(s):\n%d" % (nomap + 1))
                    recorded = True
                else:
                    if nomap == temp + 1:
                        funmap.write("-")
                        recorded = False
                        if j == len(notmapped) - 1:
                            funmap.write("%d" % (nomap + 1))
                    else:
                        if recorded:
                            funmap.write(",%d" % (nomap + 1))
                        else:
                            funmap.write("%d,%d" % (temp + 1, nomap + 1))
                            recorded = True
                temp = nomap
            notMapped.append(notmapped)

            # continuous regions
            if len(sitespecific):
                idxPairStart = np.nonzero(sitespecific[1:] - sitespecific[0:-1] != 1)[0]
                if len(idxPairStart) == 0 and len(sitespecific):  # only one continuous region exists
                    idxPairStart = np.array([0])
                    idxPairEnd = np.array([len(sitespecific) - 1])
                else:
                    if idxspecific[0] == 0:
                        idxPairStart = np.append(idxspecific[0], np.add(idxPairStart, 1))
                    else:
                        idxPairStart = np.add(idxPairStart, 1)

                    idxPairEnd = np.zeros(len(idxPairStart), dtype=np.int)
                    for j in range(len(idxPairStart) - 1):
                        p = idxPairStart[j]
                        c = 0
                        while sitespecific[p + c + 1] == sitespecific[p + c] + 1:
                            c += 1
                        idxPairEnd[j] = p + c

                sitePairStart = sitespecific[idxPairStart]
                sitePairEnd = sitespecific[idxPairEnd]
                sitePairEnd[-1] = sitespecific[-1]
                siteCandidates.append(np.vstack((sitePairStart, sitePairEnd)))

            else:   # no usable fragment
                siteCandidates.append(np.zeros((2, 0)))
            funmap.write("\n\n")
    return siteCandidates, notMapped

