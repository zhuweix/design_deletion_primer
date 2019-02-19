import os
from biofile import simple_fasta_load
from sequence_lib import rc_seq
from sequence_lib import calc_tm

"""
Design deletion primer of given C.glabrata strain
For each gene, upstream primer started ~ 500 bp to the ATG;
downstream primer started ~ 500 bp to the stop codon.
The last 8 bp of the foward primer in both streams were unqiue
to avoid off-target alignments. In addition, one primer pair 
within the gene was also designed to verify gene deletion.
Tm ~ 55C
GC ratio 30%-60%
length 20-30 bp
"""

# global parameters
gc_low = 0.3
gc_high = 0.6

min_primerlen = 20
max_primerlen = 20

unique_size = 10

# Overlap with deletion marker, NAT gene
leftnat = ''
righnat = ''


def generate_unique_kmer(fasta: str, prefix: str, ksize=unique_size):
    '''
    Generate locations for unqiue kmer locations in the strain
    :param fasta: name of the fasta file
    :param prefix: prefix of the output file
    :param ksize: size for unique kmer
    :return: None
    '''
    # first filter for unqiue kmers
    _, seqs = simple_fasta_load(fasta)
    seqs = map(lambda x: x.upper(), seqs)
    kmer_count = {}
    # perform kmer counting from W and C strand to
    # filter for kmer unique in both strand
    # palindromic kmer saved as single copy
    # non-palindromic kmer in both directions were saved.
    for s in seqs:
        for i in range(len(s) - ksize + 1):
            kmer = s[i: i + ksize]
            kmer_count.setdefault(kmer, 0)
            kmer_count[kmer] += 1
        rs = rc_seq(s)
        for i in range(len(s) - ksize + 1):
            kmer = rs[i: i + ksize]
            kmer_count.setdefault(kmer, 0)
            kmer_count[kmer] += 1
    unique_kmer = []
    for kmer, count in kmer_count.items():
        if count == 2:
            if kmer == rc_seq(kmer):
                unique_kmer.append(kmer)
        elif count == 1:
            unique_kmer.append(kmer)
    print('%d unqiue %d mers identified!' % (len(unique_kmer), ksize))
    # save kmers
    with open(prefix + '.unique.%d.mer.tab' % ksize, 'w') as filep:
        filep.write('\n'.join(unique_kmer))





def test_function():
    # Test the function by BG2 data
    os.chdir('/home/zhuwei/data/development/collection/data/190214_deletion_primer/')
    generate_unique_kmer('./data/bg2.05.fa', './data/bg2.05')


if __name__ == '__main__':
    test_function()