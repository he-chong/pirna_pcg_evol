from Bio import SeqIO
import os, glob, gzip


def pirna_size_selecting(fq_dir, out_dir, phred_cutoff, len_lower, len_upper):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    for fq_file in glob.glob(os.path.join(fq_dir, '*.cutadapt.fastq.gz')):
        out_base = os.path.basename(fq_file).replace('.fastq.gz', '.refined.fastq.gz')
        out_fq = os.path.join(out_dir, out_base)
        if not os.path.isfile(out_fq):
            with gzip.open(fq_file, 'rt') as fq_handle, gzip.open(out_fq, 'wt') as out_handle:
                for record in SeqIO.parse(fq_handle, 'fastq'):
                    qualified = True
                    for each in record.letter_annotations['phred_quality']:
                        if each < phred_cutoff:
                            qualified = False
                    if len(record.seq) >= len_lower and len(record.seq) <= len_upper and qualified:
                        record.description = ''
                        record.name = ''
                        SeqIO.write(record, out_handle, 'fastq')


def refine_human():
    fq_dir = r'../../generated_data/adaptcut_pirnas/human'
    out_dir = r'../../generated_data/refined_pirnas/human'
    pirna_size_selecting(fq_dir, out_dir, 10, 25, 32)


if __name__ == '__main__':
    refine_human()
