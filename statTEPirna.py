import os
from collections import defaultdict
from statPssTEPreference import build_TE_dict, loc_TE_associated


ref_chromosomes = ['chr{}'.format(i) for i in list(range(22))+['X', 'Y']]

def stat_TE_pirna(pirna_pos_file, repeatmasker_file, out_file):
    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    TE_dict = build_TE_dict(repeatmasker_file)

    TE_associated_pirnas = set()
    all_pirnas = set()
    counted = 0
    with open(pirna_pos_file) as pirna_pos_handle:
        head = True
        for line in pirna_pos_handle:
            if head:
                head = False
                continue
            line_split = line.strip().split('\t')
            pirna = line_split[0]
            loc_info = '{}({})'.format(*line_split[1:3])
            chromosome = loc_info.split(':')[0]
            if chromosome in ref_chromosomes:
                all_pirnas.add(pirna)
                if loc_TE_associated(loc_info, TE_dict):
                    TE_associated_pirnas.add(pirna)
            counted += 1
            if (counted-1) % 100000 == 0:
                print(counted) 
    with open(out_file, 'w') as out_handle:
        TE_associated_count = len(TE_associated_pirnas)
        total = len(all_pirnas)
        out_handle.write('{}\t{}\t{}\n'.format(TE_associated_count, total, TE_associated_count/total))
        print(TE_associated_count, total, TE_associated_count/total)


def stat_human():
    pirna_pos_file = '../generated_data/pirna_pos/human_healthy/human_pirnas.healthy.dedup.retained.pirna_pos.txt'
    repeatmasker_file = '../original_data/repeatmasker_files/hg38.fa.out'
    out_file = '../generated_data/TE_associated_stat/human/proportion_TE_pirnas.txt'
    stat_TE_pirna(pirna_pos_file, repeatmasker_file, out_file)


if __name__ == '__main__':
    stat_human()
