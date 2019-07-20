# padlock design pipeline for multiplexed assay with multiple probes per target in cDNA-based expression profiling
# Xiaoyan, 2017

from lib import checkinput
from lib import screenseq
from lib import formatrefseq
from lib import parblast
from lib import readblast
from lib import createoutput
from lib import distributeprobes
from lib import finalizeprobes


if __name__ == "__main__":
    try:
        # get keyboard inputs and prepare sequences
        designpars, outpars, genepars, designinput = \
            checkinput.getdesigninput()
        # fmt = checkinput.checkformat(genepars[2])

        # Tm screening
        Tm, siteChopped = screenseq.thresholdtm(designinput[1], designinput[2], outpars[1], designpars)

        # blast
        formatrefseq.blastdb(designpars[0])
        parblast.continueblast(siteChopped, designinput[1], outpars[1], designpars)

        # specific targets
        siteCandidates, notMapped = readblast.getcandidates(siteChopped, designinput[1], outpars, designpars[1], designinput[3])
        createoutput.writetargetfile(designinput, siteCandidates, Tm, designpars[1], outpars, '3.AllSpecificTargets_')

        # non-overlapping candidates
        targets, targetpos, mapTmlist = distributeprobes.asmanyprobes(siteCandidates, siteChopped, designinput[2], designpars)

        # correct positions
        targets, targetpos, notMapped, Tm = finalizeprobes.correctpos(
            designinput[0], targets, targetpos, notMapped, mapTmlist, Tm, siteChopped)

        # write genes with no candidates
        createoutput.emptyentries(targets, genepars[2], outpars)

        # fill up linker sequence and write
        probes = finalizeprobes.assembleprobes(targets, genepars, designpars[1])
        createoutput.writeprobefile(
            genepars[0], genepars[2], probes, Tm, targetpos, targets, outpars, designpars[1], '4.NonOverlappingProbes_')

        # remove targets cannot be found in database
        finallist = finalizeprobes.removeunmapped(notMapped, targetpos, genepars[2], targets, Tm, probes)
        createoutput.writeprobefile(
            genepars[0],genepars[2], finallist[0], finallist[1], finallist[2], finallist[3],
            outpars, designpars[1], '5.ProbesDBMappable_')

        # prioritize sequences without homopolymers and randomly select the fixed number of probes per gene (if applicable)
        if len(designpars[5]):
            sublist = finalizeprobes.selectprobes(int(designpars[5]), finallist, genepars[2])
            createoutput.writeprobefile(
                genepars[0], genepars[2], sublist[0], sublist[1], sublist[2], sublist[3],
                outpars, designpars[1], '6.ProbesRandomSubsetN=' + designpars[5] + '_')

        print("All finished!")

    except:
        import sys
        print (sys.exc_info()[0])
        import traceback
        print (traceback.format_exc())
    finally:
        print("Press Enter to continue ...")
        input()
