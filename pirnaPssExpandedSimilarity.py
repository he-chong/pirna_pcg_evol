import os, glob, random, subprocess, multiprocessing
from Bio import SeqIO


PIRNA_LOC = 5
PSS_LOC = 6

def pirna_pss_pair_file_size(pirna_pss_pair_file):
    with open(pirna_pss_pair_file) as pirna_pss_pair_handle:
        head = True
        sampled = False
        num = 0
        for line in pirna_pss_pair_handle:
            if head:
                head = False
                continue
            if line.strip():
                num += 1
    print(pirna_pss_pair_file, num)


def aln_each_expanded_pirna_pss_pair(expanded_seq_file, expanded_aln_dir):
    expanded_aln_file = os.path.join(expanded_aln_dir, os.path.basename(expanded_seq_file))
    with open(expanded_aln_file, 'w') as expanded_aln_handle:
        p = subprocess.Popen(['mafft', expanded_seq_file], stdout=expanded_aln_handle)
        p.communicate()


def aln_expanded_pirna_pss_pairs(seq_dir, aln_dir, processes):
    pool = multiprocessing.Pool(processes)
    for seq_file in glob.glob(os.path.join(seq_dir, '*.fa')):
        pool.apply_async(aln_each_expanded_pirna_pss_pair, (seq_file, aln_dir))
    pool.close()
    pool.join()


def expanded_pirna_pss_pairs(pirna_pss_pair_file, up_expanded_seq_dir, up_expanded_aln_dir, down_expanded_seq_dir, down_expanded_aln_dir, overall_expanded_seq_dir, overall_expanded_aln_dir, genome_fa, file_size, repeat_num, sample_size, expanded_size, processes):
    pirna_pss_strand = -1 if 'antisense' in pirna_pss_pair_file else 1
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fa, 'fasta'))
    line_num_list = list(range(file_size))
    for i in range(repeat_num):
        up_repeat_seq_dir = os.path.join(up_expanded_seq_dir, 'repeat_{}'.format(i+1))
        up_repeat_aln_dir = os.path.join(up_expanded_aln_dir, 'repeat_{}'.format(i+1))
        down_repeat_seq_dir = os.path.join(down_expanded_seq_dir, 'repeat_{}'.format(i+1))
        down_repeat_aln_dir = os.path.join(down_expanded_aln_dir, 'repeat_{}'.format(i+1))
        overall_repeat_seq_dir = os.path.join(overall_expanded_seq_dir, 'repeat_{}'.format(i+1))
        overall_repeat_aln_dir = os.path.join(overall_expanded_aln_dir, 'repeat_{}'.format(i+1))
        if not os.path.isdir(up_repeat_seq_dir):
            os.makedirs(up_repeat_seq_dir)
        if not os.path.isdir(up_repeat_aln_dir):
            os.makedirs(up_repeat_aln_dir)
        if not os.path.isdir(down_repeat_seq_dir):
            os.makedirs(down_repeat_seq_dir)
        if not os.path.isdir(down_repeat_aln_dir):
            os.makedirs(down_repeat_aln_dir)
        if not os.path.isdir(overall_repeat_seq_dir):
            os.makedirs(overall_repeat_seq_dir)
        if not os.path.isdir(overall_repeat_aln_dir):
            os.makedirs(overall_repeat_aln_dir)
        random.shuffle(line_num_list)
        sampled_line_num_list = sorted(line_num_list[:sample_size])
        # print(sampled_line_num_list)
        with open(pirna_pss_pair_file) as pirna_pss_pair_handle:
            head = True
            sampled = False
            num = 0
            for line in pirna_pss_pair_handle:
                if head:
                    head = False
                    continue
                if num in sampled_line_num_list:
                    sampled = True
                line_split = line.strip().split()
                pirna_loc_list = line_split[PIRNA_LOC].split(';')
                random.shuffle(pirna_loc_list)
                pirna_loc = pirna_loc_list[0]
                pss_loc = line_split[PSS_LOC]
                
                num += 1
                if pirna_loc == 'null':
                    continue
                if sampled and '~' not in pirna_loc and '~' not in pss_loc:
                    print(num)
                    pirna_chrom = pirna_loc.split(':')[0]
                    if pirna_chrom not in genome_dict:
                        continue
                    pirna_start = int(pirna_loc.split(':')[1].split('-')[0])
                    pirna_end = int(pirna_loc.split(':')[1].split('-')[1].split('(')[0])
                    pirna_strand = int(pirna_loc.split('(')[-1].strip(')'))
                    pirna_expanded_start = pirna_start - expanded_size
                    pirna_expanded_end = pirna_end + expanded_size
                    pss_chrom = pss_loc.split(':')[0]
                    if pss_chrom not in genome_dict:
                        continue
                    pss_start = int(pss_loc.split(':')[1].split('-')[0])
                    pss_end = int(pss_loc.split(':')[1].split('-')[1].split('(')[0])
                    pss_strand = int(pss_loc.split('(')[-1].strip(')'))
                    pss_expanded_start = pss_start - expanded_size
                    pss_expanded_end = pss_end + expanded_size

                    if pirna_chrom == pss_chrom and ( \
                       (pirna_start >= pss_start and pirna_start <= pss_end) or (pirna_end >= pss_start and pirna_end <= pss_end) or \
                       (pss_start >= pirna_start and pss_start <= pirna_end) or (pss_end >= pirna_start and pss_end <= pirna_end)):
                        continue

                    if pirna_strand == 1:
                        up_pirna_expanded_record = genome_dict[pirna_chrom][pirna_expanded_start-1:pirna_end]
                        up_pirna_expanded_record.id = '{}:{}-{}({})'.format(pirna_chrom, pirna_expanded_start, pirna_end, pirna_strand)
                    else:
                        up_pirna_expanded_record = genome_dict[pirna_chrom][pirna_start-1:pirna_expanded_end]
                        up_pirna_expanded_record.id = '{}:{}-{}({})'.format(pirna_chrom, pirna_start, pirna_expanded_end, pirna_strand)
                        up_pirna_expanded_record.seq = up_pirna_expanded_record.seq.reverse_complement()
                    up_pirna_expanded_record.description = ''
                    up_pirna_expanded_record.name = ''                    

                    if pss_strand*pirna_pss_strand == 1:
                        up_pss_expanded_record = genome_dict[pss_chrom][pss_expanded_start-1:pss_end]
                        up_pss_expanded_record.id = '{}:{}-{}({})'.format(pss_chrom, pss_expanded_start, pss_end, pss_strand)
                    else:
                        up_pss_expanded_record = genome_dict[pss_chrom][pss_start-1:pss_expanded_end]
                        up_pss_expanded_record.id = '{}:{}-{}({})'.format(pss_chrom, pss_start, pss_expanded_end, pss_strand)
                        up_pss_expanded_record.seq = up_pss_expanded_record.seq.reverse_complement()
                    up_pss_expanded_record.description = ''
                    up_pss_expanded_record.name = ''

                    if pirna_strand == 1:
                        down_pirna_expanded_record = genome_dict[pirna_chrom][pirna_start-1:pirna_expanded_end]
                        down_pirna_expanded_record.id = '{}:{}-{}({})'.format(pirna_chrom, pirna_start, pirna_expanded_end, pirna_strand)
                    else:
                        down_pirna_expanded_record = genome_dict[pirna_chrom][pirna_expanded_start-1:pirna_end]
                        down_pirna_expanded_record.id = '{}:{}-{}({})'.format(pirna_chrom, pirna_expanded_start, pirna_end, pirna_strand)
                        down_pirna_expanded_record.seq = down_pirna_expanded_record.seq.reverse_complement()
                    down_pirna_expanded_record.description = ''
                    down_pirna_expanded_record.name = ''                    

                    if pss_strand*pirna_pss_strand == 1:
                        down_pss_expanded_record = genome_dict[pss_chrom][pss_start-1:pss_expanded_end]
                        down_pss_expanded_record.id = '{}:{}-{}({})'.format(pss_chrom, pss_start, pss_expanded_end, pss_strand)
                    else:
                        down_pss_expanded_record = genome_dict[pss_chrom][pss_expanded_start-1:pss_end]
                        down_pss_expanded_record.id = '{}:{}-{}({})'.format(pss_chrom, pss_expanded_start, pss_end, pss_strand)
                        down_pss_expanded_record.seq = down_pss_expanded_record.seq.reverse_complement()                    
                    down_pss_expanded_record.description = ''
                    down_pss_expanded_record.name = ''

                    overall_pirna_expanded_record = genome_dict[pirna_chrom][pirna_expanded_start-1:pirna_expanded_end]
                    overall_pirna_expanded_record.id = '{}:{}-{}({})'.format(pirna_chrom, pirna_expanded_start, pirna_expanded_end, pirna_strand)
                    if pirna_strand == -1:
                        overall_pirna_expanded_record.seq = overall_pirna_expanded_record.seq.reverse_complement()
                    overall_pirna_expanded_record.description = ''
                    overall_pirna_expanded_record.name = ''                    

                    overall_pss_expanded_record = genome_dict[pss_chrom][pss_expanded_start-1:pss_expanded_end]
                    overall_pss_expanded_record.id = '{}:{}-{}({})'.format(pss_chrom, pss_expanded_start, pss_expanded_end, pss_strand)
                    if pss_strand*pirna_pss_strand == -1:
                        overall_pss_expanded_record.seq = overall_pss_expanded_record.seq.reverse_complement()
                    overall_pss_expanded_record.description = ''
                    overall_pss_expanded_record.name = ''

                    up_expanded_seq_file = os.path.join(up_repeat_seq_dir, '{}_{}_{}.fa'.format(num, pirna_loc, pss_loc))
                    down_expanded_seq_file = os.path.join(down_repeat_seq_dir, '{}_{}_{}.fa'.format(num, pirna_loc, pss_loc))
                    overall_expanded_seq_file = os.path.join(overall_repeat_seq_dir, '{}_{}_{}.fa'.format(num, pirna_loc, pss_loc))                  
                    SeqIO.write([up_pirna_expanded_record, up_pss_expanded_record], up_expanded_seq_file, 'fasta')
                    SeqIO.write([down_pirna_expanded_record, down_pss_expanded_record], down_expanded_seq_file, 'fasta')
                    SeqIO.write([overall_pirna_expanded_record, overall_pss_expanded_record], overall_expanded_seq_file, 'fasta')

                    sampled = False    # reset sampled as False
                else:
                    # print(line)
                    pass
        aln_expanded_pirna_pss_pairs(up_repeat_seq_dir, up_repeat_aln_dir, processes)
        aln_expanded_pirna_pss_pairs(down_repeat_seq_dir, down_repeat_aln_dir, processes)
        aln_expanded_pirna_pss_pairs(overall_repeat_seq_dir, overall_repeat_aln_dir, processes)


def human_expanded():
    sense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/human/sense.txt'
    sense_up_expanded_seq_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/sense/up_expanded_seq'
    sense_up_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/sense/up_expanded_aln'
    sense_down_expanded_seq_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/sense/down_expanded_seq'
    sense_down_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/sense/down_expanded_aln'
    sense_expanded_seq_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/sense/expanded_seq'
    sense_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/sense/expanded_aln'

    antisense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/human/antisense.txt' 
    antisense_up_expanded_seq_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/antisense/up_expanded_seq'
    antisense_up_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/antisense/up_expanded_aln'
    antisense_down_expanded_seq_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/antisense/down_expanded_seq'
    antisense_down_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/antisense/down_expanded_aln'
    antisense_expanded_seq_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/antisense/expanded_seq'
    antisense_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_100/human/antisense/expanded_aln'

    sense_file_size = 160153074
    antisense_file_size = 138648496
    repeat_num = 30
    sample_size = 200
    expanded_size = 100
    processes = 4
    genome_fa = '../original_data/genomes/hg38.fa'
    pirna_pss_pair_file_size(sense_pirna_pss_pair_file)
    pirna_pss_pair_file_size(antisense_pirna_pss_pair_file)
    expanded_pirna_pss_pairs(sense_pirna_pss_pair_file, sense_up_expanded_seq_dir, sense_up_expanded_aln_dir, sense_down_expanded_seq_dir, sense_down_expanded_aln_dir, sense_expanded_seq_dir, sense_expanded_aln_dir, genome_fa, sense_file_size, repeat_num, sample_size, expanded_size, processes)
    expanded_pirna_pss_pairs(antisense_pirna_pss_pair_file, antisense_up_expanded_seq_dir, antisense_up_expanded_aln_dir, antisense_down_expanded_seq_dir, antisense_down_expanded_aln_dir, antisense_expanded_seq_dir, antisense_expanded_aln_dir, genome_fa, antisense_file_size, repeat_num, sample_size, expanded_size, processes)


if __name__ == '__main__':
    human_expanded()
