import os
import re
import subprocess
import multiprocessing
import gzip
from Bio import SeqIO


def autoRunBwa(pirna_seq, target_rna_fa, working_dir, edit_dist, thread=None):
    if not os.path.isfile(target_rna_fa):
        raise IOError('{} does not exist'.format(target_rna_fa))

    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    if not thread:
        thread = int(multiprocessing.cpu_count()/2)

    sai_file = os.path.join(working_dir, 'mm{}.sai'.format(edit_dist))
    sam_file = os.path.join(working_dir, 'mm{}.sam'.format(edit_dist))
    raw_bam_file = os.path.join(working_dir, 'mm{}.bam'.format(edit_dist))
    sorted_bam_file = os.path.join(working_dir, 'mm{}.sorted.bam'.format(edit_dist))

    if not (os.path.isfile(target_rna_fa+'.amb') and os.path.isfile(target_rna_fa+'.ann') and os.path.isfile(target_rna_fa+'.bwt')):
        p_bwa_index = subprocess.Popen(['bwa', 'index', '-a', 'bwtsw', target_rna_fa])
        p_bwa_index.communicate()

    if not os.path.isfile(sai_file):
        with open(sai_file, 'w') as sai_handle:
            p_bwa_aln = subprocess.Popen(['bwa', 'aln', '-t', str(thread), '-n', str(edit_dist), '-N', '-o', '0', '-R', '999999', target_rna_fa, pirna_seq], stdout=sai_handle)
            p_bwa_aln.communicate()

    if not os.path.isfile(sam_file):
        with open(sam_file, 'w') as sam_handle:
            p_bwa_samse = subprocess.Popen(['bwa', 'samse', '-n', '999999', target_rna_fa, sai_file, pirna_seq], stdout=sam_handle)
            p_bwa_samse.communicate()

    if not os.path.isfile(raw_bam_file):
        p_samtools_view = subprocess.Popen(['samtools', 'view', '-bS', sam_file, '-o', raw_bam_file])
        p_samtools_view.communicate()

    if not os.path.isfile(sorted_bam_file):
        p_samtools_sort = subprocess.Popen(['samtools', 'sort', raw_bam_file, '-o', sorted_bam_file])
        p_samtools_sort.communicate()

    if not os.path.isfile(sorted_bam_file+'.bai'):
        p_samtools_index = subprocess.Popen(['samtools', 'index', sorted_bam_file])
        p_samtools_index.communicate()


def run_for_human_simple():
    pirna_seq = os.path.abspath('../generated_data/cleaned_pirnas/human_pirnas.healthy.dedup.retained.fastq.gz')
    pcg_fa = os.path.abspath('../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa')
    working_dir = os.path.abspath('../generated_data/pirna_pss_bwa/human')
    if not os.path.isfile(pirna_seq):
        raise Exception('{} does not exist.'.format(pirna_seq))
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    autoRunBwa(pirna_seq, pcg_fa, working_dir, edit_dist=4, thread=12)


def run_for_human_complicated():
    pirna_seq = os.path.abspath('../generated_data/cleaned_pirnas/human_pirnas.healthy.dedup.retained.fastq.gz')
    pcg_fa = os.path.abspath('../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa')
    working_dir = os.path.abspath('../generated_data/pirna_pss_bwa/human')
    if not os.path.isfile(pirna_seq):
        raise Exception('{} does not exist.'.format(pirna_seq))
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    autoRunBwa(pirna_seq, pcg_fa, working_dir, edit_dist=5, thread=12)


def run_for_d_melanogaster():
    pirna_seq = os.path.abspath('../original_data/pirna_sequences/dme.v3.0.fa.gz')
    pcg_fa = os.path.abspath('../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna.info_seperated.fa')
    working_dir = os.path.abspath('../generated_data/pirna_pss_bwa/d_melanogaster')
    if not os.path.isfile(pirna_seq):
        raise Exception('{} does not exist.'.format(pirna_seq))
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    autoRunBwa(pirna_seq, pcg_fa, working_dir, edit_dist=4, thread=6)


if __name__ == '__main__':
    run_for_human_simple()
    run_for_human_complicated()
    run_for_d_melanogaster()