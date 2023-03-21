# for multi_padlock_design
# Xiaoyan, 2017


import numpy as np


def asmanyprobes(siteCandidates, originalmap, sequences, designpars):
    Probes = []
    ProbesPos = []
    armlength = designpars[1]
    interval = designpars[2]
    for i, sites in enumerate(siteCandidates):
        probes = []
        probespos = []
        try:
            current = siteCandidates[i][1][-1]
            seq = sequences[i]
            while current >= siteCandidates[i][0][0]:  # possible to squeeze more in
                probes.append(seq[current : current + armlength * 2])
                probespos.append(current)
                current -= armlength * 2 + interval
                idx_max_start = np.amax(
                    np.append(np.nonzero(siteCandidates[i][0] <= current), -1)
                )
                idx_min_end = np.amin(
                    np.append(np.nonzero(siteCandidates[i][1] >= current), len(seq))
                )
                while (
                    idx_max_start != idx_min_end and current >= siteCandidates[i][0][0]
                ):
                    current -= 1
                    idx_max_start = np.amax(np.nonzero(siteCandidates[i][0] <= current))
                    idx_min_end = np.amin(np.nonzero(siteCandidates[i][1] >= current))
        except IndexError:  # no usable fragment
            probes = []
            probespos = []

        Probes.append(probes)
        ProbesPos.append(probespos)

    # map to Tm list index
    mapTmlistnew = []
    for i, pos in enumerate(ProbesPos):
        tempmap = []
        for j in pos:
            tempmap.append(originalmap[i].index(j))
        mapTmlistnew.append(tempmap)

    return Probes, ProbesPos, mapTmlistnew
