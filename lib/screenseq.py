# for multi_padlock_design
# Xiaoyan, 2017


import numpy as np
import math
import multiprocessing
from lib import createoutput


def chopseq(seq, window, step):
    """Moving window to chop target sequence"""
    seqChopped = []
    while len(seq) >= window:
        seqChopped.append(seq[0:window])
        seq = seq[step:]
    return seqChopped


def calculatetm(seq):
    """Calculate Tm of a target candidate, nearest neighbor model"""
    NNlist = chopseq(seq, 2, 1)
    NNtable = [
        "AA",
        "AC",
        "AG",
        "AT",
        "CA",
        "CC",
        "CG",
        "CT",
        "GA",
        "GC",
        "GG",
        "GT",
        "TA",
        "TC",
        "TG",
        "TT",
    ]
    NNendtable = ["A", "C", "G", "T"]
    NNcount = np.zeros(16)
    NNend = np.zeros(4)
    for c, NN in enumerate(NNtable):
        NNcount[c] = NNlist.count(NN)
    for c, NN in enumerate(NNendtable):
        NNend[c] = seq[0].count(NN)
    # numbers below from Sugimoto et al. NAR (1996)
    NNEnthalpy = np.array(
        [
            -8.0,
            -9.4,
            -6.6,
            -5.6,
            -8.2,
            -10.9,
            -11.8,
            -6.6,
            -8.8,
            -10.5,
            -10.9,
            -9.4,
            -6.6,
            -8.8,
            -8.2,
            -8.0,
        ]
    )
    NNEntropy = np.array(
        [
            -21.9,
            -25.5,
            -16.4,
            -15.2,
            -21.0,
            -28.4,
            -29.0,
            -16.4,
            -23.5,
            -26.4,
            -28.4,
            -25.5,
            -18.4,
            -23.5,
            -21.0,
            -21.9,
        ]
    )
    NNendEnthalpy = np.array([0.6, 0.6, 0.6, 0.6])
    NNendEntropy = np.array([-9.0, -9.0, -9.0, -9.0])

    sumEnthalpy = np.sum(np.multiply(NNcount, NNEnthalpy)) + np.sum(
        np.multiply(NNend, NNendEnthalpy)
    )
    sumEntropy = np.sum(np.multiply(NNcount, NNEntropy)) + np.sum(
        np.multiply(NNend, NNendEntropy)
    )
    Tm = (sumEnthalpy * 1000) / (
        sumEntropy + (1.9872 * math.log(1e-7))
    ) - 273.15  # oligo concentration: 1e-7 M
    sumSalt = 0.075 + (3.795 * 0.01**0.5)  # monovalent: 0.075 M, bivalent: 0.01 M
    Tm += 16.6 * math.log10(sumSalt)  # salt correction
    Tm -= 0.72 * 20  # formamide correction
    return Tm


def runscreen(argin):
    """Screen Tm and write sequences that fulfill Tm requirement to file, to use as blast input"""
    dirname, header, seq, armlen, t1, t2 = argin
    fname = createoutput.blastinfilename(dirname, header)

    Tm = []
    siteChopped = []

    c = 0
    listSeqChopped = chopseq(seq, armlen * 2, 1)
    for j, seqChopped in enumerate(listSeqChopped):
        tm = calculatetm(seqChopped)
        Tm.append(tm)
        if t1 < tm < t2:  # Tm thresholding
            c += 1
            siteChopped.append(j)

            # write files that can be used as input in blastn
            with open(fname + "_" + str(c) + ".fasta", "w") as fblast:
                fblast.write(">target_%d\n%s\n" % (c, seqChopped))
    return Tm, siteChopped


def thresholdtm(headers, sequences, dirname, designpars):
    """Parallel process of Tm thresholding for input sequences"""

    # pack up inputs for multiprocessing
    inputs = []
    for c, i in enumerate(headers):
        inputs.append(
            (dirname, i, sequences[c], designpars[1], designpars[3], designpars[4])
        )
    pool = multiprocessing.Pool(6)

    # unpack output from pool
    Tm = []
    siteChopped = []
    for argout in pool.map(runscreen, inputs):
        Tm.append(argout[0])
        siteChopped.append(argout[1])
    return Tm, siteChopped
