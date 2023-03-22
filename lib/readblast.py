# for multi_padlock_design
# Xiaoyan, 2017

import os
from lib import createoutput
import numpy as np

notmapped = []
funmap = []


def readblastout(file, armlength, variants):
    """Read the results from blast"""
    global funmap
    global notmapped
    specific = True
    mappable = False
    with open(file, "r") as fh:
        noblast = True
        for line in fh:
            noblast = False
            if specific:
                # scores = []
                if not ("XR" in line or "XM" in line):  # skip all predicted transcripts
                    if "|" in line:
                        columns = line.split("|")
                        if len(columns) <= 3:
                            hit = columns[1].split(".", 1)[0]
                            # if columns[1][:2] == 'NM' or columns[1][:2] == 'NR':    # check if the hit is transcript (NM) or non-coding RNA (NR), remove all predicted (XM)
                            scores = columns[-1].split(",")
                        else:
                            hit = columns[3].split(".", 1)[0]
                            # if columns[3][:2] == 'NM' or columns[3][:2] == 'NR':
                            scores = columns[-1].split(",")

                    else:
                        columns = line.split(",")
                        # if columns[1][:2] == 'NM' or columns[1][:2] == 'NR':
                        hit = columns[1].split(".", 1)[0]
                        scores = columns[2:]

                    if len(scores):
                        # remove the first empty column (somehow appears in some db and blast versions)
                        if "" in scores:
                            scores.remove("")

                        if (
                            2 * armlength * 0.5 < int(scores[1]) < 2 * armlength
                            and float(scores[0]) > 80
                            and int(scores[4]) < armlength - 4
                            and int(scores[5]) > armlength + 5
                        ):
                            # more than 50% coverage, 80% homology, and non-target sequence covers ligation site +- 5
                            specific = False

                        if float(scores[0]) == 100 and int(scores[1]) == 2 * armlength:
                            mappable = True
                            if len(variants) and hit not in variants:
                                if isinstance(variants, list):
                                    with open(
                                        os.path.join(
                                            os.path.dirname(file), "homology.txt"
                                        ),
                                        "a",
                                    ) as fsimilar:
                                        fsimilar.write("%s,%s\n" % (hit, variants[0]))
                                else:
                                    with open(
                                        os.path.join(
                                            os.path.dirname(file), "homology.txt"
                                        ),
                                        "a",
                                    ) as fsimilar:
                                        fsimilar.write("%s,%s\n" % (hit, variants))
                                specific = False
        if not mappable:
            with open(file[0:-10] + ".fasta", "r") as f:
                seq = f.readlines()
                funmap.write("Could not map sequence in " + file[:-10] + "!\n")
                funmap.write(seq[1])
                notmapped.append(int(file[:-10].split("_")[-1]) - 1)

    # if no blast results returned, ignore the sequence (can be due to error in blast+ or due to no similar sequence)
    if noblast:
        specific = False

    return specific


def getcandidates(listSiteChopped, headers, dirnames, armlength, accession):
    """Get specific fragments"""
    global notmapped
    global funmap
    siteCandidates = []
    notMapped = []
    t = dirnames[1].split("TempFolder")[1]
    with open(os.path.join(dirnames[0], "2.Unmappable_" + t + ".txt"), "w") as funmap:
        for i, sites in enumerate(listSiteChopped):
            funmap.write("%s\n" % headers[i])

            fname = createoutput.blastinfilename(dirnames[1], headers[i])

            notmapped = []
            blast_bw = []
            for j, target in enumerate(sites):
                fblast = fname + "_" + str(j + 1) + "_blast.txt"
                blast_bw.append(readblastout(fblast, armlength, accession[i]))

            # find sequences that are specific enough
            idxspecific = np.nonzero(blast_bw)[0]
            tempCandidates = np.array(sites)
            sitespecific = tempCandidates[idxspecific]

            # write unmappable sites
            notmapped = tempCandidates[notmapped]
            recorded = False
            for j, nomap in enumerate(notmapped):
                if j == 0:
                    funmap.write(
                        "\nUnmapped sequence starting position(s):\n%d" % (nomap + 1)
                    )
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
                if len(idxPairStart) == 0 and len(
                    sitespecific
                ):  # only one continuous region exists
                    idxPairStart = np.array([0])
                    idxPairEnd = np.array([len(sitespecific) - 1])
                else:
                    if idxspecific[0] == 0:
                        idxPairStart = np.append(
                            idxspecific[0], np.add(idxPairStart, 1)
                        )
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

            else:  # no usable fragment
                siteCandidates.append(np.zeros((2, 0)))
            funmap.write("\n\n")
    return siteCandidates, notMapped
