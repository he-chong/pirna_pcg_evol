import os, glob, subprocess


def pirna_adaptors_removing(reads_root, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    raw_fq_list = glob.glob(os.path.join(reads_root, '*.fastq.gz'))
    for raw_fq in raw_fq_list:
        out_fq_base = os.path.basename(raw_fq).replace('.fastq', '.cutadapt.fastq')
        out_fq = os.path.join(out_dir, out_fq_base)
        if not os.path.isfile(out_fq):
            print('Raw fastq:', raw_fq)
            p_minion = subprocess.Popen(['minion', 'search-adapter', '-i', raw_fq], stdout=subprocess.PIPE)
            adaptor_seq = None
            for line_bytes in p_minion.stdout.readlines():
                line = line_bytes.decode('utf-8')
                if line.startswith('sequence='):
                    adaptor_seq = line.strip().split('=')[-1]
                    break
            p_minion.communicate()
            print('Adaptor:', adaptor_seq)        
            p_cut = subprocess.Popen(['cutadapt', '-j', '3', '-a', adaptor_seq, '-o', out_fq, raw_fq])
            p_cut.communicate() 


if __name__ == '__main__':
    reads_root = '../../original_data/pirna_reads/human'
    out_dir = '../../generated_data/adaptcut_pirnas/human'
    pirna_adaptors_removing(reads_root, out_dir)