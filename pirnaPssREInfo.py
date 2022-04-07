import os
import multiprocessing
from collections import defaultdict


PIRNA_LOC = 5
PSS_LOC = 6

QUERY_SEQUENCE = 4
POS_BEGIN = 5
POS_END = 6
RE_CLASS_FAMILY = 10

def split_sites(sites_str):
    return [int(i) for i in sites_str.split(';')]


def merge_block_dict(block_dict):
    merged_dict = {}
    print('Merging Dict')
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
    print('OK\n')
    return merged_dict


def sort_block_dict(block_dict):
    sorted_dict = {}
    for block_host_id in block_dict.keys():
        sorted_blocks = sorted(block_dict[block_host_id], key=lambda i:i[0])
        sorted_dict[block_host_id] = sorted_blocks
    return sorted_dict


def re_associated(loc_info, re_dict):
    if loc_info != 'null':
        for each_loc in loc_info.split(';'):
            for loc_segment in each_loc.split('~'):
                chromosome = loc_segment.split(':')[0]
                segment_start = int(loc_segment.split(':')[1].split('(')[0].split('-')[0])
                segment_end = int(loc_segment.split(':')[1].split('(')[0].split('-')[1])
                if chromosome in re_dict:
                    re_records = re_dict[chromosome]
                    i = int(len(re_records)/2)
                    j = -1
                    while len(re_records) > 0:
                        j = i
                        re_record = re_records[i]
                        re_geno_start = int(re_record[POS_BEGIN]) + 1    # change to 1-baseed
                        re_geno_end = int(re_record[POS_END])
                        if (segment_start >= re_geno_start and segment_start <= re_geno_end) or \
                           (segment_end >= re_geno_start and segment_end <= re_geno_end) or \
                           (segment_start <= re_geno_start and segment_end >= re_geno_end):
                            try:
                                re_class, re_family = re_record[RE_CLASS_FAMILY].split(r'/')
                                return (re_class, re_family, str(segment_start - re_geno_start))
                            except ValueError:
                                return (re_record[RE_CLASS_FAMILY], re_record[RE_CLASS_FAMILY], str(segment_start - re_geno_start))
                        re_records = re_records[:i] if segment_start < re_geno_start else re_records[i+1:]
                        i = int(len(re_records)/2)
        return ('None', 'None', '.')
    else:
        return ('Null', 'Null', '.')


def focused_overlapping(chromosome, geno_start, geno_end, focused_dict):    # geno_start and geno_end are 1-based, focused_dict should be sorted
    if chromosome in focused_dict:
        focused_blocks = focused_dict[chromosome]
        i = int(len(focused_blocks)/2)
        j = -1
        while len(focused_blocks) > 0:
            j = i
            focused_block = focused_blocks[i]
            block_geno_start, block_geno_end = focused_block 
            if (geno_start >= block_geno_start and geno_start <= block_geno_end) or \
               (geno_end >= block_geno_start and geno_end <= block_geno_end) or \
               (geno_start <= block_geno_start and geno_end >= block_geno_end):
                return True
            focused_blocks = focused_blocks[:i] if geno_start < block_geno_start else focused_blocks[i+1:]
            i = int(len(focused_blocks)/2)
        return False
    else:
        return None


def build_pcg_dict(gene_info_file):
    chrom_index = 5
    block_geno_start_index = (10,12,14)
    block_geno_end_index = (11,13,15)

    pcg_dict = defaultdict(list)
    with open(gene_info_file) as gene_info_handle:
        print('Building PCG Dict')
        head = True
        for line in gene_info_handle:
            if head:
                head = False
                continue
            line_split = line.strip().split('\t')
            # print(line_split)
            chromosome = 'chr' + line_split[chrom_index]
            block_bounds_in_genome = []
            for each_start_index, each_end_index in zip(block_geno_start_index, block_geno_end_index):
                if line_split[each_start_index] and line_split[each_end_index]:
                    block_geno_start_list = sorted(split_sites(line_split[each_start_index]))
                    block_geno_end_list = sorted(split_sites(line_split[each_end_index]))
                    block_bounds_in_genome += list(zip(block_geno_start_list, block_geno_end_list))   # 1-based
            for block_geno_start, block_geno_end in block_bounds_in_genome:
                pcg_dict[chromosome].append((block_geno_start, block_geno_end))
        print('OK\n')
    return pcg_dict


def build_pirna_dict(pirna_loc_file):
    pirna_loc_index = 1
    pirna_dict = defaultdict(list)
    with open(pirna_loc_file) as pirna_loc_handle:
        print('Building piRNA Dict')
        head = True
        for line in pirna_loc_handle:
            if head:
                head = False
                continue
            line_split = line.strip().split('\t')
            pirna_loc = line_split[pirna_loc_index]
            chromosome = pirna_loc.split(':')[0]
            pirna_start = int(pirna_loc.split(':')[1].split('-')[0])
            pirna_end = int(pirna_loc.split(':')[1].split('-')[1].split('(')[0])
            pirna_dict[chromosome].append((pirna_start, pirna_end))
        print('OK\n')
    return pirna_dict


def build_focused_re_dict(repeatmasker_file, foucused_dict):
    focused_re_dict = defaultdict(list)
    with open(repeatmasker_file) as repeatmasker_handle:
        print('Building focused-RE Dict')
        head = True
        for line in repeatmasker_handle:
            if line.strip() == '':
                head = False
                continue
            if head:
                continue
            line_split = line.strip().split()
            chromosome = line_split[QUERY_SEQUENCE]
            # print(line_split)
            re_geno_start = int(line_split[POS_BEGIN]) + 1  # change to 1-baseed
            re_geno_end = int(line_split[POS_END])
            if focused_overlapping(chromosome, re_geno_start, re_geno_end, foucused_dict):
                focused_re_dict[chromosome].append(line_split)
        print('OK\n')
    return focused_re_dict


def pirna_pss_overlapping(pirna_loc, pss_loc):
    overlapping = False
    if pirna_loc != 'null':
        pirna_chrom = pirna_loc.split(':')[0]
        pirna_start = int(pirna_loc.split(':')[1].split('-')[0])
        pirna_end = int(pirna_loc.split(':')[1].split('-')[1].split('(')[0])
    for pss_loc_segment in pss_loc.split('~'):
        pss_chrom = pss_loc_segment.split(':')[0]
        pss_start = int(pss_loc_segment.split(':')[1].split('-')[0])
        pss_end = int(pss_loc_segment.split(':')[1].split('-')[1].split('(')[0])
        if pirna_loc != 'null' and pirna_chrom == pss_chrom and ( \
           (pirna_start >= pss_start and pirna_start <= pss_end) or (pirna_end >= pss_start and pirna_end <= pss_end) or \
           (pss_start >= pirna_start and pss_start <= pirna_end) or (pss_end >= pirna_start and pss_end <= pirna_end)):
            overlapping = True
            break
    return overlapping


def compile_pirna_pss_re_info(pirna_pss_pair_file, rmsk_file, gene_info_file, pirna_loc_file, out_file):
    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    pirna_dict = build_pirna_dict(pirna_loc_file)
    merged_pirna_dict = merge_block_dict(pirna_dict)
    pirna_re_dict = build_focused_re_dict(rmsk_file, merged_pirna_dict)
    pcg_dict = build_pcg_dict(gene_info_file)
    sorted_pcg_dict = sort_block_dict(pcg_dict)
    pcg_re_dict = build_focused_re_dict(rmsk_file, sorted_pcg_dict)
    with open(pirna_pss_pair_file) as pirna_pss_pair_handle, open(out_file, 'w') as out_handle:
        head = True
        for pirna_pss_record in pirna_pss_pair_handle:
            if head:
                head = False
                continue
            pirna_pss_record_split = pirna_pss_record.strip().split('\t')
            pirna_loc_info = pirna_pss_record_split[PIRNA_LOC]
            pss_loc_info = pirna_pss_record_split[PSS_LOC]

            pirna_loc_list = pirna_loc_info.split(';')

            overlapping = False
            for each_pirna_loc in pirna_loc_list:
                if pirna_pss_overlapping(each_pirna_loc, pss_loc_info):
                    overlapping = True
                    break    
            if overlapping:
                continue    # if the PSS is overlapped with the piRNA, the piRNA-PSS pair will be excluded.

            pirna_re_info = re_associated(pirna_loc_info, pirna_re_dict)
            pss_re_info = re_associated(pss_loc_info, pcg_re_dict)
            out_info = (pirna_pss_record_split[0],) + pirna_re_info + pss_re_info
            out_handle.write('\t'.join(out_info))
            out_handle.write('\n')
        print('OK\n')


def compile_human():
    repeatmasker_file = '../original_data/repeatmasker_files/hg38.fa.out'
    gene_info_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    pirna_loc_file = '../generated_data/pirna_pos/human_healthy/human_pirnas.healthy.dedup.retained.pirna_pos.txt'
    sense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/human/sense.txt'
    antisense_pirna_pss_pair_file = '../generated_data/pirna_pss_pairs/human/antisense.txt'
    sense_out_file = '../generated_data/piRNA_PSS_RE_info_repeatmasker/human/sense.txt'
    antisense_out_file = '../generated_data/piRNA_PSS_RE_info_repeatmasker/human/antisense.txt'
    compile_pirna_pss_re_info(sense_pirna_pss_pair_file, repeatmasker_file, gene_info_file, pirna_loc_file, sense_out_file)
    compile_pirna_pss_re_info(antisense_pirna_pss_pair_file, repeatmasker_file, gene_info_file, pirna_loc_file, antisense_out_file)


if __name__ == '__main__':
    compile_human()

    

