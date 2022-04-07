import os, glob
import gzip

def pirna_redundant_removing(in_fq, out_fq, info_file):
    out_dir = os.path.dirname(out_fq)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isfile(out_fq):
        seq_dict = {}    # record unique sequences
        with gzip.open(in_fq, 'rt') as in_handle, gzip.open(out_fq, 'wt') as out_handle:
            record = None
            is_seq = False
            for line in in_handle:
                if line.startswith('@SRR'):
                    record = line               # start a new record
                    head = line.strip()
                    is_seq = True               # next line is the sequence line
                else:
                    record += line
                    if is_seq:
                        seq = line.strip()
                        is_seq = False
                    else:
                        if line.strip() != '+': # the only possibility for this is that this line is the quality line (last line of each record)
                            if seq not in seq_dict:
                                seq_dict[seq] = ['Not_to_output', 1]
                            else:
                                seq_dict[seq][1] += 1
                                if seq_dict[seq][1] == 2:    # Require appearing more than one time
                                    seq_dict[seq][0] = head 
                                    out_handle.write(record)
    if not os.path.isfile(info_file):
        with open(info_file, 'w') as info_handle:
            info_handle.write('\t'.join(['piRNA_sequence', 'representative_read_id', 'read_count']))
            info_handle.write('\n')
            for seq, (head, count) in seq_dict.items():
                if count > 1:
                    info_handle.write('\t'.join([seq, head, str(count)]))
                    info_handle.write('\n')
    print('ok')


if __name__ == '__main__':
    in_fq = '../../generated_data/merged_pirnas/human_pirnas.healthy.fastq.gz'
    out_fq = '../../generated_data/merged_pirnas/human_pirnas.healthy.dedup.fastq.gz'
    info_file = '../../generated_data/merged_pirnas/human_pirnas.healthy.info.txt'
    pirna_redundant_removing(in_fq, out_fq, info_file)
