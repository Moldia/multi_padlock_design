# for multi_padlock_design
# Xiaoyan, 2017

import datetime
import os
import numpy as np

def usetime():
    """ Use current time to create temporary folder"""
    t = str(datetime.datetime.now())[:19]
    t = t.translate(None, '-')
    t = t.translate(None, ':')
    t = t.translate(None, ' ')
    return t


def tempdir(outdir):
    t = usetime()
    if outdir[:-1] == '/':
        outdir += 'TempFolder' + t
    else:
        outdir += '/TempFolder' + t
    return outdir


def blastinfilename(dirname, filename):
    toavoid = [':', '/', '\\', '[', ']', '?', '"', ' ', '<', '>']
    for i in toavoid:
        filename = filename.translate(None, i)

    try:
        filenamefrag = filename.split('|')
        filename = os.path.join(dirname, filenamefrag[3] + filenamefrag[4] + '_target')  # access number and name
    except:
        filename = filename.translate(None, '|')
        filename = os.path.join(dirname, filename + '_target')  # first 20 characters
    return filename


def writetargetfile(designinput, sites, Tm, armlength, dirnames, fname):
    t = dirnames[1].split('TempFolder')[1]
    headers = designinput[1]
    sequences = designinput[2]
    with open(os.path.join(dirnames[0], fname + t + '.csv'), 'w') as f:
        f.write("target,Tm,startpos\n")
        for i, header in enumerate(headers):
            f.write("%s\n" % header)
            try:
                for j in range(len(sites[i][0])):
                    for k in range(sites[i][0][j], sites[i][1][j]+1):
                        f.write("%s,%f,%d\n"
                                % (sequences[i][k:k+armlength*2], Tm[i][k], k+1))
            except:
                pass
            f.write("\n")


def writeprobefile(headers, probes, Tm, targetpos, targets, dirnames, armlength, fname):
    """ Write file with original header, target sequence, Tm and final probe sequence """
    t = dirnames[1].split('TempFolder')[1]
    with open(os.path.join(dirnames[0], fname + t + '.csv'), 'w') as f:
        f.write("target,Tm,startpos,endpos,padlock\n")
        for i, header in enumerate(headers):
            f.write("%s\n" % header)
            idxsort = np.argsort(targetpos[i])
            for j in idxsort:
                f.write("%s,%f,%d,%d,%s\n"
                        % (targets[i][j], Tm[i][j],
                           targetpos[i][j] + 1, targetpos[i][j] + armlength*2,
                           probes[i][j]))
            f.write("\n")

