# for multi_padlock_design
# Xiaoyan, 2016-8-2


def readfasta(fastafile):
    success = False
    sequences = []
    headers = []

    try:
        with open(str(fastafile), 'r') as f:
            sequence = [line.rstrip('\n') for line in f]
            sequence.append([])

        # find the fasta header lines starting with '>'
        headerlines = [i for i in range(len(sequence)) if (sequence[i] and sequence[i][0] == '>')]

        # format sequences
        if len(headerlines) == 0:
            print('No sequence in fasta format found. Try again.')

        elif len(headerlines) >= 1:
            for i in range(len(headerlines) - 1):
                sequences.append(''.join(sequence[headerlines[i]+1:headerlines[i+1]-1]))
                headers.append(sequence[headerlines[i]])

            # the last entry (the only one in the case of one sequence in the whole file)
            sequences.append(''.join(sequence[headerlines[-1] + 1:len(sequence) - 1]))
            headers.append(sequence[headerlines[-1]])
            success = True

        del sequence

    except IOError:
        print('Fasta file not found. Try again.')

    return success, headers, sequences

