import os, glob, multiprocessing, subprocess
import pysam
import gzip
from Bio import SeqIO


def run_bwa(pirna, target_rna_fa, working_dir, name, thread=None):
    if not os.path.isfile(pirna):
        raise IOError('{} does not exist'.format(pirna))

    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    if not thread:
        thread = int(multiprocessing.cpu_count()/2)

    sai_file = os.path.join(working_dir, name+'.sai')
    sam_file = os.path.join(working_dir, name+'.sam')
    raw_bam_file = os.path.join(working_dir, name+'.bam')
    sorted_bam_file = os.path.join(working_dir, name+'.sorted.bam')

    if not (os.path.isfile(target_rna_fa+'.amb') and os.path.isfile(target_rna_fa+'.ann') and os.path.isfile(target_rna_fa+'.bwt')):
        p_bwa_index = subprocess.Popen(['bwa', 'index', '-a', 'bwtsw', target_rna_fa])
        p_bwa_index.communicate()

    if not os.path.isfile(sai_file):
        with open(sai_file, 'w') as sai_handle:
            p_bwa_aln = subprocess.Popen(['bwa', 'aln', '-t', str(thread), '-l', '4', '-n', '0', '-N', '-o', '0', target_rna_fa, pirna], stdout=sai_handle)
            p_bwa_aln.communicate()

    if not os.path.isfile(sam_file):
        with open(sam_file, 'w') as sam_handle:
            print(sai_file)
            print(target_rna_fa)
            p_bwa_samse = subprocess.Popen(['bwa', 'samse', target_rna_fa, sai_file, pirna], stdout=sam_handle)
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


def non_pirna_removing(raw_pirna, excluded_list, working_dir):
    print(excluded_list)
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    pirna_base = os.path.basename(raw_pirna)
    retained_pirna = os.path.join(working_dir, pirna_base.replace('fastq.gz','retained.fastq.gz'))
    if not os.path.isfile(retained_pirna):
        excluded = set()
        for excluded_fa in excluded_list:
            print(excluded_fa)
            excluded_base = os.path.basename(excluded_fa)
            run_name = 'exclude_{}'.format(excluded_base)
            exclude_bam = os.path.join(working_dir, run_name + '.sorted.bam')
            if not os.path.isfile(exclude_bam):
                run_bwa(raw_pirna, excluded_fa, working_dir, run_name, thread=3)      
            with pysam.AlignmentFile(exclude_bam, 'rb') as exclude_aln:     
                for segment in exclude_aln.fetch():
                    if segment.flag == 0 or segment.flag == 16:
                        excluded.add(segment.query_name)
            print('excluded:', len(excluded))
        with gzip.open(raw_pirna, 'rt') as in_handle, gzip.open(retained_pirna, 'wt') as out_handle:
            record = None
            is_seq = False
            output = False
            for line in in_handle:
                if line.startswith('@SRR'):
                    if record != None and output == True: #output the last record
                        out_handle.write(record)
                        output = False
                    record = line
                    if line.strip().strip('@') not in excluded:
                        output = True
                else:
                    record += line


if __name__ == '__main__':
    raw_pirna = r'../../generated_data/merged_pirnas/human_pirnas.healthy.dedup.fastq.gz'
    excluded_list = ['../../generated_data/non_pirnas/human_non_pirnas.fa.gz']
    working_dir = r'../../generated_data/cleaned_pirnas'
    non_pirna_removing(raw_pirna, excluded_list, working_dir)