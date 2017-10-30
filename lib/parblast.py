# for multi_padlock_design
# Xiaoyan, 2017


import os
import subprocess
import config
from lib import createoutput


Processes = []
NextProcess = 0
sites = 0
fname = str()
species = []
notMapped = []


def newblast():
    """ Start a new blastn subprocess if there is work to do """
    global Processes
    global NextProcess
    global fname
    fastadir = (config.fastadir_mouse, config.fastadir_human)

    if NextProcess < len(sites):
        blastf = fname + '_' + str(NextProcess+1) + '.fasta'

        if not os.path.isfile(fname + '_' + str(NextProcess+1) + '_blast.txt'):     # skip already existing files
            txtcmd = ' '.join(('blastn', '-query', '"' + blastf + '"',
                               '-db', '"' + os.path.join(fastadir[species], ('mouse', 'human')[species] + '.transcriptome' + '"'),
                               '-outfmt 10',
                               '-out ', '"' + fname + '_' + str(NextProcess+1) + '_blast.txt' + '"',
                               '-word_size 7 -strand plus'))

            blastnprocess = subprocess.Popen(txtcmd, shell=True)
            Processes.append(blastnprocess)
        NextProcess += 1


def runningblast():
    """ Check any running processes and start new ones if there are spare slots """
    global Processes
    global NextProcess

    for p in range(len(Processes)-1, -1, -1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
    while (len(Processes) < config.max_MSA_processes) and (NextProcess < len(sites)):  # more to do and some spare slots
        newblast()


def continueblast(listSiteChopped, headers, dirname, designpars):
    global Processes
    global NextProcess
    global sites
    global fname
    global species

    if designpars[0] == "mouse":
        species = 0
    elif designpars[0] == "human":
        species = 1

    print ("Starting blast..")
    for i, sites in enumerate(listSiteChopped):
        print(headers[i])
        NextProcess = 0
        fname = createoutput.blastinfilename(dirname, headers[i])

        runningblast()  # start the max processes running
        while len(Processes) > 0:   # still going on
            runningblast()
    print ("Blast finished!")


