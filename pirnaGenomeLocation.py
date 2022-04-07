import os, multiprocessing, subprocess
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


def out_put_in_genome(pirna_name, chromosome, start, end, pirna_to_genome_strand, pirna_record, out_handle):
    pirna_record_id = pirna_record.id
    if pirna_name != pirna_record_id:
        raise Exception('{} != {}'.format(pirna_name,pirna_record_id))
    pirna_genome_loc = '{}:{}-{}'.format(chromosome, start, end)
    out_handle.write('\t'.join([pirna_name, pirna_genome_loc, str(pirna_to_genome_strand)]))
    out_handle.write('\n')


def find_pirna_genome_location(pirna_seq, genome_fa, working_dir, reference_chromosomes, thread=None): # output an table recording putative genome locations for each pirna
    pirna_base = os.path.basename(pirna_seq)
    out_base = '.'.join(pirna_base.split('.')[:-2]) + '.pirna_pos.txt'
    genome_out_table = os.path.join(working_dir, out_base)
    genome_bwa_name = 'pirna_genome_mapping_{}'.format(pirna_base)
    genome_bam_file = os.path.join(working_dir, genome_bwa_name+'.bam')

    if not os.path.isfile(genome_bam_file):
        run_bwa(pirna_seq, genome_fa, working_dir, genome_bwa_name, thread)

    with pysam.AlignmentFile(genome_bam_file, 'rb') as genome_bam_aln, gzip.open(pirna_seq, 'rt') as pirna_handle, open(genome_out_table, 'w') as genome_out_handle:
        genome_out_handle.write('\t'.join(['piRNA', 'piRNA_genome_location', 'piRNA_to_genome_strand']))
        genome_out_handle.write('\n')
        records = SeqIO.parse(pirna_handle, 'fastq') if 'fastq' in pirna_base else SeqIO.parse(pirna_handle, 'fasta')
        for pirna_record, segment in zip(records, genome_bam_aln):
            pirna_name = segment.query_name
            if (segment.flag == 0 or segment.flag == 16) and pirna_name:
                chromosome = genome_bam_aln.get_reference_name(segment.reference_id)
                if chromosome in reference_chromosomes:
                    start = segment.reference_start + 1 #change to be 1-based
                    end = start + segment.query_length - 1
                    pirna_to_genome_strand = 1 if segment.flag == 0 else -1
                    out_put_in_genome(pirna_name, chromosome, start, end, pirna_to_genome_strand, pirna_record, genome_out_handle)

                if segment.has_tag('XA'):
                    xa_info = segment.get_tag('XA')
                    for xa_target in xa_info.split(';'):
                        if xa_target:
                            xa_target_split = xa_target.split(',')
                            chromosome = xa_target_split[0]
                            if chromosome in reference_chromosomes:
                                start = int(xa_target_split[1][1:])
                                end = start + segment.query_length - 1
                                pirna_to_genome_strand = 1 if xa_target_split[1][0] == '+' else -1
                                out_put_in_genome(pirna_name, chromosome, start, end, pirna_to_genome_strand, pirna_record, genome_out_handle)

def find_human():
    pirna_fq = r'../generated_data/cleaned_pirnas/human_pirnas.healthy.dedup.retained.fastq.gz'
    genome_fa = r'../original_data/genomes/hg38.fa' 
    position_working_dir = r'../generated_data/pirna_pos/human_healthy' 
    reference_chromosomes = set(['chr'+str(i+1) for i in range(22)] + ['X']) 
    find_pirna_genome_location(pirna_fq, genome_fa, position_working_dir, reference_chromosomes)


def find_d_melanogaster():
    pirna_fa = '../original_data/pirna_sequences/dme.v3.0.fa.gz'
    genome_fa = '../original_data/genomes/dm6.fa' 
    position_working_dir = '../generated_data/pirna_pos/d_melanogaster'  
    reference_chromosomes = set(['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY']) 
    find_pirna_genome_location(pirna_fa, genome_fa, position_working_dir, reference_chromosomes, thread=2)


if __name__ == '__main__':
    find_human()
    find_d_melanogaster()