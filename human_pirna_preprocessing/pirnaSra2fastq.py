import os
import subprocess


def pirna_sra2fastq(sra_root):
	abs_sra_root = os.path.abspath(sra_root)
	os.chdir(abs_sra_root)
	for each_sra_base in os.listdir(abs_sra_root):
		each_sra_dir = os.path.join(abs_sra_root, each_sra_base)
		sra_file = os.path.join(each_sra_dir, each_sra_base+'.sra')
		fq_file = each_sra_base+'.fastq.gz'
		if not os.path.isfile(fq_file):
			p = subprocess.Popen(['fastq-dump', '--gzip', sra_file])
			p.communicate()


if __name__ == '__main__':
	sra_root = r'../../original_data/pirna_reads/human'
	pirna_sra2fastq(sra_root)