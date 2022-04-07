import os
from collections import defaultdict
from Bio import SeqIO


PIRNA_GENO_LOC = 5
PSS_GENO_LOC = 6

# GENO_NAME = 5
# GENO_START = 6
# GENO_END = 7
# STRAND = 9
# REP_NAME = 10
# REP_CLASS = 11
# REP_FAMILY = 12
# REP_START = 13
# REP_END = 14

QUERY_SEQUENCE = 4
POS_BEGIN = 5
POS_END = 6
RE_CLASS_FAMILY = 10

def split_sites(sites_str):
    return [int(i) for i in sites_str.split(';')]


def loc_TE_associated(loc_info, TE_dict):
    for each_loc in loc_info.split(';'):
        for loc_segment in each_loc.split('~'):
            chromosome = loc_segment.split(':')[0]
            segment_start = int(loc_segment.split(':')[1].split('(')[0].split('-')[0])
            segment_end = int(loc_segment.split(':')[1].split('(')[0].split('-')[1])
            if chromosome in TE_dict:
                TE_records = TE_dict[chromosome]
                i = int(len(TE_records)/2)
                j = -1
                while len(TE_records) > 0:
                    j = i
                    TE_geno_start, TE_geno_end = TE_records[i]
                    if (segment_start >= TE_geno_start and segment_start <= TE_geno_end) or \
                       (segment_end >= TE_geno_start and segment_end <= TE_geno_end) or \
                       (segment_start <= TE_geno_start and segment_end >= TE_geno_end):
                        return True
                    TE_records = TE_records[:i] if segment_start < TE_geno_start else TE_records[i+1:]
                    i = int(len(TE_records)/2)
    return False


def block_TE_overlap(block_geno_start, block_geno_end, TE_records_in_same_chrom):
    TE_records = TE_records_in_same_chrom
    TE_overlap_blocks = []
    found = False
    i = int(len(TE_records)/2)
    j = -1
    while len(TE_records) > 0:
        j = i
        TE_geno_start, TE_geno_end = TE_records[i]
        if (block_geno_start >= TE_geno_start and block_geno_start <= TE_geno_end) or \
           (block_geno_end >= TE_geno_start and block_geno_end <= TE_geno_end) or \
           (block_geno_start <= TE_geno_start and block_geno_end >= TE_geno_end):
            found = True
            TE_overlap_start = block_geno_start if block_geno_start >= TE_geno_start and block_geno_start <= TE_geno_end else TE_geno_start
            TE_overlap_end = block_geno_end if block_geno_end >= TE_geno_start and block_geno_end <= TE_geno_end else TE_geno_end
            TE_overlap_blocks.append((TE_overlap_start, TE_overlap_end))
            break
        TE_records = TE_records[:i] if block_geno_start < TE_geno_start else TE_records[i+1:]
        i = int(len(TE_records)/2)

    if found:
        for m in range(i+1, len(TE_records)):
            TE_geno_start, TE_geno_end = TE_records[m]
            # print(geno_end, re_geno_start, re_geno_end)
            if TE_geno_start > block_geno_end:
                break
            else:
                TE_overlap_start = block_geno_start if block_geno_start >= TE_geno_start and block_geno_start <= TE_geno_end else TE_geno_start
                TE_overlap_end = block_geno_end if block_geno_end >= TE_geno_start and block_geno_end <= TE_geno_end else TE_geno_end
                TE_overlap_blocks.append((TE_overlap_start, TE_overlap_end))

        for n in range(i-1, -1, -1):
            TE_geno_start, TE_geno_end = TE_records[n]
            if TE_geno_end < block_geno_start:
                break
            else:
                TE_overlap_start = block_geno_start if block_geno_start >= TE_geno_start and block_geno_start <= TE_geno_end else TE_geno_start
                TE_overlap_end = block_geno_end if block_geno_end >= TE_geno_start and block_geno_end <= TE_geno_end else TE_geno_end
                TE_overlap_blocks.append((TE_overlap_start, TE_overlap_end))
    return TE_overlap_blocks


def build_TE_dict(repeatmasker_file, target_family=None):
    TE_dict = defaultdict(list)
    with open(repeatmasker_file) as repeatmasker_handle:
        print('repeatmasker file dicting')
        head = True
        for line in repeatmasker_handle:
            if line.strip() == '':
                head = False
                continue
            if head:
                continue
            line_split = line.strip().split()
            try:
                re_class, re_family = line_split[RE_CLASS_FAMILY].split(r'/')
            except ValueError:
                re_class = line_split[RE_CLASS_FAMILY]
                re_family = re_class
            if  re_class != 'Simple_repeat' and re_class != 'Satellite' and re_class != 'Low_complexity' and re_class != 'Unknown' and 'RNA' not in re_class:
                if target_family == None:
                    chromosome = line_split[QUERY_SEQUENCE]
                    TE_geno_start = int(line_split[POS_BEGIN]) + 1    # change to 1-baseed
                    TE_geno_end = int(line_split[POS_END])
                    TE_dict[chromosome].append((TE_geno_start, TE_geno_end))
                else:
                    if re_family == target_family:
                        chromosome = line_split[QUERY_SEQUENCE]
                        TE_geno_start = int(line_split[POS_BEGIN]) + 1    # change to 1-baseed
                        TE_geno_end = int(line_split[POS_END])
                        TE_dict[chromosome].append((TE_geno_start, TE_geno_end))
        print('OK\n')
    return TE_dict


def merge_block_dict(block_dict):
    merged_dict = {}
    for block_host_id in block_dict.keys():
        merged_blocks = []
        sorted_blocks = sorted(block_dict[block_host_id], key=lambda i:i[0])
        for start, end in sorted_blocks:
            if not merged_blocks:
                merged_blocks.append([start, end])
            else:
                last_start, last_end = merged_blocks[-1]
                if start >= last_start and start <= last_end + 1:
                    if end >= last_end:
                        merged_blocks[-1][1] = end
                else:
                    merged_blocks.append([start, end])
        merged_dict[block_host_id] = merged_blocks
    return merged_dict


def build_pcg_TE_dict(TE_dict, gene_info):
    chrom_index = 5
    block_geno_start_index = (10,12,14)
    block_geno_end_index = (11,13,15)

    total_TE_len = 0
    total_len = 0
    pcg_TE_dict = defaultdict(list)
    with open(gene_info) as gene_handle:
        head = True
        for line in gene_handle:
            if head:
                head = False
                continue
            line_split = line.strip().split('\t')
            # print(line_split)
            chromosome = 'chr' + line_split[chrom_index]
            if chromosome not in TE_dict:
                continue
            block_bounds_in_genome = []
            for each_start_index, each_end_index in zip(block_geno_start_index, block_geno_end_index):
                if line_split[each_start_index] and line_split[each_end_index]:
                    block_geno_start_list = sorted(split_sites(line_split[each_start_index]))
                    block_geno_end_list = sorted(split_sites(line_split[each_end_index]))
                    block_bounds_in_genome += list(zip(block_geno_start_list, block_geno_end_list))   # 1-based
            transcript_strand = int(line_split[-1])
            TE_records = TE_dict[chromosome]
            for block_geno_start, block_geno_end in block_bounds_in_genome:
                for block_TE_overlap_start, block_TE_overlap_end in block_TE_overlap(block_geno_start, block_geno_end, TE_records):
                    total_TE_len += block_TE_overlap_end - block_TE_overlap_start + 1
                    pcg_TE_dict[chromosome].append((block_TE_overlap_start, block_TE_overlap_end))
                total_len += block_geno_end - block_geno_start + 1
    total_non_TE_len = total_len - total_TE_len
    return total_TE_len, total_non_TE_len, total_len, pcg_TE_dict


def sort_block_dict(block_dict):
    sorted_dict = {}
    for block_host_id in block_dict.keys():
        sorted_blocks = sorted(block_dict[block_host_id], key=lambda i:i[0])
        sorted_dict[block_host_id] = sorted_blocks
    return sorted_dict


def stat_pss_TE_preference(pirna_pss_pair_files, repeatmasker_file, gene_info, out_file):   
    print(pirna_pss_pair_files)
    print('piRNA-PSS pair file dicting')

    TE_dict = build_TE_dict(repeatmasker_file)
    Alu_dict = build_TE_dict(repeatmasker_file, target_family='Alu')
    merged_TE_dict = merge_block_dict(TE_dict)
    merged_Alu_dict = merge_block_dict(Alu_dict)
    total_TE_len, total_non_TE_len, total_len, pcg_TE_dict = build_pcg_TE_dict(merged_TE_dict, gene_info)
    total_Alu_len, total_non_Alu_len, total_len, pcg_Alu_dict = build_pcg_TE_dict(merged_Alu_dict, gene_info)
    sorted_pcg_TE_dict = sort_block_dict(pcg_TE_dict)
    sorted_pcg_Alu_dict = sort_block_dict(pcg_Alu_dict)

    pss_dict = defaultdict(list)
    TE_associated_count = 0
    Alu_associated_count = 0
    total = 0
    for pirna_pss_pair_file in pirna_pss_pair_files:
        head = True
        with open(pirna_pss_pair_file) as pirna_pss_pair_handle:        
            for line in pirna_pss_pair_handle:
                if head:
                    head = False
                    continue
                line_split = line.strip().split('\t')
                pirna_geno_loc = line_split[PIRNA_GENO_LOC]
                if pirna_geno_loc == 'null':
                    continue
                pss_geno_loc = line_split[PSS_GENO_LOC]
                if loc_TE_associated(pss_geno_loc, sorted_pcg_TE_dict):
                    TE_associated_count += 1
                if loc_TE_associated(pss_geno_loc, sorted_pcg_Alu_dict):
                    Alu_associated_count += 1
                for segment in pss_geno_loc.split('~'):
                    chromosome = segment.split(':')[0]
                    geno_start = int(segment.split(':')[1].split('(')[0].split('-')[0])
                    geno_end = int(segment.split(':')[1].split('(')[0].split('-')[1])
                    pss_dict[chromosome].append((geno_start, geno_end))
                if total % 100000 == 0:
                    print(total)
                total += 1
    non_TE_associated_count = total - TE_associated_count
    non_Alu_associated_count = total - Alu_associated_count
    merged_pss_dict = merge_block_dict(pss_dict)    # each block represent a continuous region similar to some piRNAs

    total_pss_TE_overlap_len = 0
    total_pss_Alu_overlap_len = 0
    total_pss_len = 0
    for chromosome in merged_pss_dict.keys():
        if chromosome not in sorted_pcg_TE_dict:
            continue
        TE_records = sorted_pcg_TE_dict[chromosome]
        try:
            Alu_records = sorted_pcg_Alu_dict[chromosome]
        except KeyError:
            Alu_records = []
        for pss_geno_start, pss_geno_end in merged_pss_dict[chromosome]:
            for pss_TE_overlap_start, pss_TE_overlap_end in block_TE_overlap(pss_geno_start, pss_geno_end, TE_records):
                total_pss_TE_overlap_len += pss_TE_overlap_end - pss_TE_overlap_start + 1
            if Alu_records:
                for pss_Alu_overlap_start, pss_Alu_overlap_end in block_TE_overlap(pss_geno_start, pss_geno_end, Alu_records):
                    total_pss_Alu_overlap_len += pss_Alu_overlap_end - pss_Alu_overlap_start + 1
            total_pss_len += pss_geno_end - pss_geno_start + 1
    total_pss_non_TE_overlap_len = total_pss_len - total_pss_TE_overlap_len
    total_pss_non_Alu_overlap_len = total_pss_len - total_pss_Alu_overlap_len

    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(out_file, 'w') as out_handle:
        out_handle.write('P(TE): {}\t{}\t{}\n'.format(total_TE_len, total_len, total_TE_len/total_len))
        # out_handle.write('P(PSS): {}\t{}\t{}\n'.format(total_pss_len, total_len, total_pss_len/total_len))
        # out_handle.write('P(TE|PSS): {}\t{}\t{}\n'.format(total_pss_TE_overlap_len, total_pss_len, total_pss_TE_overlap_len/total_pss_len))
        # out_handle.write('P(TE|non-PSS): {}\t{}\t{}\n'.format(total_pss_TE_overlap_len, total_non_TE_len, total_pss_TE_overlap_len/total_non_TE_len))
        out_handle.write('P(PSS|TE): {}\t{}\t{}\n'.format(total_pss_TE_overlap_len, total_TE_len, total_pss_TE_overlap_len/total_TE_len))
        out_handle.write('P(PSS|non-TE): {}\t{}\t{}\n'.format(total_pss_non_TE_overlap_len, total_non_TE_len, total_pss_non_TE_overlap_len/total_non_TE_len))
        out_handle.write('P(PSS|Alu): {}\t{}\t{}\n'.format(total_pss_Alu_overlap_len, total_Alu_len, total_pss_Alu_overlap_len/total_Alu_len))
        out_handle.write('P(PSS|non-Alu): {}\t{}\t{}\n'.format(total_pss_non_Alu_overlap_len, total_non_Alu_len, total_pss_non_Alu_overlap_len/total_non_Alu_len))
        out_handle.write('Average_density_of_PSS_for_TE: {}\t{}\t{}\n'.format(TE_associated_count, total_TE_len, TE_associated_count/total_TE_len))
        out_handle.write('Average_density_of_PSS_for_non-TE: {}\t{}\t{}\n'.format(non_TE_associated_count, total_non_TE_len, non_TE_associated_count/total_non_TE_len))
        out_handle.write('Average_density_of_PSS_for_Alu: {}\t{}\t{}\n'.format(Alu_associated_count, total_Alu_len, Alu_associated_count/total_Alu_len))
        out_handle.write('Average_density_of_PSS_for_non-Alu: {}\t{}\t{}\n'.format(non_Alu_associated_count, total_non_Alu_len, non_Alu_associated_count/total_non_Alu_len))


def stat_human():
    sense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/human/sense.txt'
    antisense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/human/antisense.txt'
    repeatmasker_file = '../original_data/repeatmasker_files/hg38.fa.out'
    gene_info = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    sense_out_file = '../generated_data/TE_PSS_preference_stat_repeatmasker/human/sense.txt'
    antisense_out_file = '../generated_data/TE_PSS_preference_stat_repeatmasker/human/antisense.txt'
    sense_antisense_out_file = '../generated_data/TE_PSS_preference_stat_repeatmasker/human/sense_antisense.txt'

    stat_pss_TE_preference([sense_pirna_pss_pair_file], repeatmasker_file, gene_info, sense_out_file)
    stat_pss_TE_preference([antisense_pirna_pss_pair_file], repeatmasker_file, gene_info, antisense_out_file)
    stat_pss_TE_preference([sense_pirna_pss_pair_file, antisense_pirna_pss_pair_file], repeatmasker_file, gene_info, sense_antisense_out_file)


# def stat_d_melanogaster():
#     sense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/d_melanogaster/sense.txt'
#     antisense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/d_melanogaster/antisense.txt'
#     rmsk_file = '../original_data/rmsk_files/dm6_rmsk.txt'
#     gene_info = '../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna_info.txt'
#     sense_out_file = '../generated_data/TE_PSS_preference_stat/d_melanogaster/sense.txt'
#     antisense_out_file = '../generated_data/TE_PSS_preference_stat/d_melanogaster/antisense.txt'
#     sense_antisense_out_file = '../generated_data/TE_PSS_preference_stat/d_melanogaster/sense_antisense.txt'

#     stat_pss_TE_preference([sense_pirna_pss_pair_file], rmsk_file, gene_info, sense_out_file)
#     stat_pss_TE_preference([antisense_pirna_pss_pair_file], rmsk_file, gene_info, antisense_out_file)
#     stat_pss_TE_preference([sense_pirna_pss_pair_file, antisense_pirna_pss_pair_file], rmsk_file, gene_info, sense_antisense_out_file)


if __name__ == '__main__':
    stat_human()
    # stat_d_melanogaster()
