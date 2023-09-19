import numpy
import random
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


_nuc = numpy.array(["A", "T", "C", "G"])


def randSeq(length):
    seqArray = _nuc[[random.randint(0, 3) for i in range(length)]]
    return (MutableSeq("".join(seqArray)))


# parameters
N = 500000
N_reads = 500
read_size = 10000

# reference sequence
ref = randSeq(N)
# write as "reference"
SeqIO.write([SeqRecord(ref, id='1', description='')], "ref.fa", "fasta")

# genome with a deletion
samp = ref[:int(N/4-1000)] + ref[int(N/4+1000):int(N/2)] + \
    ref[int(3*N/4):]
fastq_f = open('reads.fastq', 'wt')
readid = 0
for rr in range(N_reads):
    pos = random.randint(0, len(samp) - read_size)
    read = samp[pos:(pos+read_size)]
    fastq_f.write('@r{}\n{}\n+\n{}\n'.format(readid, read, '~'*len(read)))
    readid += 1
fastq_f.close()
