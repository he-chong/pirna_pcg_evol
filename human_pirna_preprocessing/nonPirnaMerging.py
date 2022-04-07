from Bio import SeqIO
from Bio.Seq import Seq
import os, glob, gzip


def non_pirna_merging(non_pirna_dir, out_fa, taxon):
    out_dir = os.path.dirname(out_fa)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    # if not os.path.isfile(out_fa):
    with gzip.open(out_fa, 'wt') as out_handle:
        for non_pirna_fa in glob.glob(os.path.join(non_pirna_dir, '*fa.gz')):
            with gzip.open(non_pirna_fa, 'rt') as non_pirna_handle:
                for record in SeqIO.parse(non_pirna_handle, 'fasta'):
                    # print(record.description)
                    if taxon in record.description:
                        record.seq = Seq(str(record.seq).replace('U', 'T'))
                        SeqIO.write(record, out_handle, 'fasta')
        for non_pirna_fa in glob.glob(os.path.join(non_pirna_dir, '*txt.gz')):
            with gzip.open(non_pirna_fa, 'rt') as non_pirna_handle:
                for record in SeqIO.parse(non_pirna_handle, 'fasta'):
                    SeqIO.write(record, out_handle, 'fasta')
                        
if __name__ == '__main__':
    non_pirna_dir = '../../original_data/non_pirnas/human'
    out_fa = '../../generated_data/non_pirnas/human_non_pirnas.fa.gz'
    non_pirna_merging(non_pirna_dir, out_fa, 'Homo')
