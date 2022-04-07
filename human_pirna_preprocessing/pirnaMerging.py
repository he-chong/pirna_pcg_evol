from Bio import SeqIO
import os, glob
import gzip
import subprocess


def pirna_merging(pirna_dir, merged_fq, selected=None):
    merged_dir = os.path.dirname(merged_fq)
    if not os.path.isdir(merged_dir):
        os.makedirs(merged_dir)
    with gzip.open(merged_fq, 'wt') as merged_handle:
        if selected:
            pirna_fq_list = []
            for srr in selected:
                pirna_fq_list.append(os.path.join(pirna_dir, srr+'.cutadapt.refined.fastq.gz'))
        else:
            pirna_fq_list = glob.glob(os.path.join(pirna_dir, '*.cutadapt.refined.fastq.gz'))
        for pirna_fq in pirna_fq_list:
            with gzip.open(pirna_fq, 'rt') as pirna_handle:
                for line in pirna_handle:
                    merged_handle.write(line)


if __name__ == '__main__':
    pirna_dir = '../../generated_data/refined_pirnas/human'
    merged_fq = '../../generated_data/merged_pirnas/human_pirnas.healthy.fastq.gz'
    selected = [
        'SRR2156539',    # Williams et al., 2015
        'SRR2156540',
        'SRR835324',     # Ha et al., 2014
        'SRR950451',
                         # Ozata et al., 2020, Adult Group 1
        'SRR8575349',    # A2
        'SRR8575350',    # A1
        'SRR8575385',    # A5
        'SRR8575386',    # A6
        'SRR8575387',    # A3
        'SRR8575388',    # A4
                         # Ozata et al., 2020, Juvenile
        'SRR8575409',    # J1
        'SRR8575410',
        'SRR8575352',    # J2
        'SRR8575351',    # J3
    ]
    pirna_merging(pirna_dir, merged_fq, selected)