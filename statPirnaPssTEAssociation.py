# This script is written for calculating the value shown in Fig. 1A at the right side.

import os
from collections import defaultdict


PIRNA_RE_CLASS = 1
PIRNA_RE_FAMILY = 2
PSS_RE_CLASS = 4
PSS_RE_FAMILY = 5

def stat_TE_association(re_associated_info_file_list, out_file):
    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    pirna_re_class_count_dict = defaultdict(int)
    pirna_re_family_count_dict = defaultdict(int)
    pss_re_class_count_dict = defaultdict(int)
    pss_re_family_count_dict = defaultdict(int)
    total = 0
    count = 0
    for re_associated_info_file in re_associated_info_file_list:
        with open(re_associated_info_file) as re_associated_info_handle:
            head = True
            for line in re_associated_info_handle:
                if head:
                    head = False
                    continue
                line_split = line.strip().split('\t')
                # print(line_split)
                try:
                    pirna_re_class = line_split[PIRNA_RE_CLASS].strip('?')
                    pirna_re_family = line_split[PIRNA_RE_FAMILY].strip('?') + '\t' + pirna_re_class
                    pss_re_class = line_split[PSS_RE_CLASS].strip('?')
                    pss_re_family = line_split[PSS_RE_FAMILY].strip('?') + '\t' + pss_re_class 

                    pirna_re_class_count_dict[pirna_re_class] += 1
                    pirna_re_family_count_dict[pirna_re_family] += 1
                    pss_re_class_count_dict[pss_re_class] += 1
                    pss_re_family_count_dict[pss_re_family] += 1
                    total += 1
                except IndexError:
                    print(line_split)
                    # raise IndexError

    with open(out_file, 'w') as out_handle:
        out_handle.write('#piRNA RE family stat\n')
        for key, value in sorted(pirna_re_family_count_dict.items(), key=lambda i:i[1], reverse=True):
            out_handle.write('{}\t{}\n'.format(key, value))
        out_handle.write('\n')

        out_handle.write('#piRNA RE class stat\n')
        for key, value in sorted(pirna_re_class_count_dict.items(), key=lambda i:i[1], reverse=True):
            out_handle.write('{}\t{}\n'.format(key, value))
        out_handle.write('\n')

        out_handle.write('#PSS RE family stat\n')
        for key, value in sorted(pss_re_family_count_dict.items(), key=lambda i:i[1], reverse=True):
            out_handle.write('{}\t{}\n'.format(key, value))

        out_handle.write('#PSS RE class stat\n')
        for key, value in sorted(pss_re_class_count_dict.items(), key=lambda i:i[1], reverse=True):
            out_handle.write('{}\t{}\n'.format(key, value))
        out_handle.write('\n')


def stat_human():
    sense_re_info_file = '../generated_data/piRNA_PSS_RE_info_repeatmasker/human/sense.txt'
    antisense_re_info_file = '../generated_data/piRNA_PSS_RE_info_repeatmasker/human/antisense.txt'
    sense_out_file = '../generated_data/piRNA_PSS_RE_stat_repeatmasker/human/sense.txt'
    antisense_out_file = '../generated_data/piRNA_PSS_RE_stat_repeatmaskerhuman/antisense.txt'
    sense_antisense_out_file = '../generated_data/piRNA_PSS_RE_stat_repeatmasker/human/sense_antisense.txt'
    stat_TE_association([sense_re_info_file], sense_out_file)
    stat_TE_association([antisense_re_info_file], antisense_out_file)
    stat_TE_association([sense_re_info_file, antisense_re_info_file], sense_antisense_out_file)


if __name__ == '__main__':
    stat_human()