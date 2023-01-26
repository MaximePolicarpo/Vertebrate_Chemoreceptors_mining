import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'r')
FastaDroppedFile = open(sys.argv[2], 'w')
drop_cutoff = float(sys.argv[3])

if (drop_cutoff > 1) or (drop_cutoff < 0):
    print('\n Sequence drop cutoff must be in 0-1 range !\n')
    sys.exit(1)

for seqs in SeqIO.parse(FastaFile, 'fasta'):
    name = seqs.id
    seq = seqs.seq
    seqLen = len(seqs)
    gap_count = 0
    for z in range(seqLen):
        if seq[z]=='-':
            gap_count += 1
    if (gap_count/float(seqLen)) >= drop_cutoff:
        print(' %s was removed.' % name)
    else:
        SeqIO.write(seqs, FastaDroppedFile, 'fasta')

FastaFile.close()
FastaDroppedFile.close()
