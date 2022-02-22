# for multi_padlock_design
# Xiaoyan, 2018

import random


def correctpos(basepos, targets, targetpos, notMapped, mapTmlist, Tm, siteChopped):
    """ Correct fragment coordinates to full length mRNA coordinates """
    targetposnew = []
    Tmnew = []
    notMappednew = []
    targetsnew = []
    c = -1
    for base in basepos:
        if isinstance(base[0], int):    # only one variant
            c += 1
            targetsnew.append(targets[c])
            targetposnew.append(targetpos[c])
            notMappednew.append(notMapped[c])
            subTm = [tm for i, tm in enumerate(Tm[c]) if i in siteChopped[c]]
            Tmnew.append([subTm[i] for i in mapTmlist[c]])
        else:
            temptargets = []
            temppos = []
            tempnomap = []
            temptm = []

            for subbase in base:
                c += 1
                temptargets = temptargets + targets[c]
                for i, pos in enumerate(targetpos[c]):
                    temppos.append(targetpos[c][i] + subbase[0])
                for i, pos in enumerate(notMapped[c]):
                    tempnomap.append(notMapped[c][i] + subbase[0])
                subTm = [tm for i, tm in enumerate(Tm[c]) if i in siteChopped[c]]
                temptm = temptm + [subTm[i] for i in mapTmlist[c]]

            targetsnew.append(temptargets)
            targetposnew.append(temppos)
            notMappednew.append(tempnomap)
            Tmnew.append(temptm)

    return targetsnew, targetposnew, notMappednew, Tmnew


def assembleprobes(targets, genepars, armlength):
    """ Fill backbone sequences """
    linkers = genepars[1]
    Padlocks = []
    for c, probes in enumerate(targets):
        try:
            linker1 = linkers[c][0]
            linker1[0]
        except:
            linker1 = 'LINKERFIRST'

        try:
            barcode = linkers[c][1]
            barcode[0]
        except:
            barcode = 'XXXX'

        try:
            linker2 = linkers[c][2]
            linker2[0]
        except:
            linker2 = 'LINKERSECOND'

        padlocks = []
        for probe in probes:
            padlocks.append(probe[armlength:] + linker1 + barcode + linker2 + probe[0:armlength])

        Padlocks.append(padlocks)
    return Padlocks


def removeunmapped(notmapped, targetpos, headers, targets, Tm, probes):
    for i, header in enumerate(headers):
        if len(notmapped[i]):
            for j in list(reversed(range(len(targetpos[i])))):
                if targetpos[i][j] in notmapped[i]:
                    del targets[i][j]
                    del Tm[i][j]
                    del targetpos[i][j]
                    del probes[i][j]
    return (probes, Tm, targetpos, targets)


def selectprobes(n, finals, headers):
    """ Prioritize probes with no homopolymer sequences and choose randomly n candidates """
    probes = finals[0]
    Tm = finals[1]
    targetpos = finals[2]
    targets = finals[3]

    for i, header in enumerate(headers):
        if len(targets[i]) > n:
            # probes with homopolymers
            wAAAA = [c for c, j in enumerate(probes[i]) if "AAAA" in j]
            wCCCC = [c for c, j in enumerate(probes[i]) if "CCCC" in j]
            wGGGG = [c for c, j in enumerate(probes[i]) if "GGGG" in j]
            wTTTT = [c for c, j in enumerate(probes[i]) if "TTTT" in j]
            wHomo = set(wAAAA + wCCCC + wGGGG + wTTTT)

            # without homopolymers
            noHomo = list(set(range(0, len(targets[i]))) - wHomo)

            # prioritize sequence without homopolymers
            if len(noHomo) > n:
                deletei = random.sample(noHomo, len(noHomo)-n)
                deletei = deletei + list(wHomo)
            else:
                deletei = random.sample(wHomo, len(wHomo)-n+len(noHomo))

            deletei.sort(reverse=True)
            for j in deletei:
                del targets[i][j]
                del Tm[i][j]
                del targetpos[i][j]
                del probes[i][j]
    return (probes, Tm, targetpos, targets)
	
