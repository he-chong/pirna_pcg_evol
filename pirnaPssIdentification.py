from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam
import os, glob, sys
import re
import subprocess
import multiprocessing
from collections import defaultdict
import gzip
import inspect


class BeyondTheBoundsException(Exception):
    pass

class BeyondTheStartException(Exception):
    pass

class PirnaPssPair:

    def __init__(self, pirna, pss, strand, start_site=1):
        self.pirna = pirna
        self.pss = pss
        self.strand = strand
        self.all_mismatches = []
        self.mismatches_GU_1 = []
        self.mismatches_GU_2to7 = []
        self.mismatches_GU_8to21 = []
        self.mismatches_GU_outside = []
        self.mismatches_nonGU_1 = []
        self.mismatches_nonGU_2to7 = []
        self.mismatches_nonGU_8to21 = []       
        self.mismatches_nonGU_outside = []
        self.continuous_num = 0
        pirna_comparing_seq = str(pirna.seq).upper()
        pss_comparing_seq = str(pss.seq).upper() if strand == '+' else str(pss.seq.reverse_complement()).upper()
        for enum_index, (pi_nt, pss_nt) in enumerate(zip(pirna_comparing_seq, pss_comparing_seq)):
            site = enum_index + 1                     # all mismatch sites are 1-based
            real_site = enum_index + start_site       # if start_site == 2, the 1st nt corresponds to the 2nd nt in the piRNA        
            if pss_nt != pi_nt:
                self.all_mismatches.append(site)
                if (pi_nt == 'G' and pss_nt == 'A') or (pi_nt == 'T' and pss_nt == 'C'):
                    if real_site == 1:
                        self.mismatches_GU_1.append(site)
                    elif real_site >= 2 and real_site <= 7:
                        self.mismatches_GU_2to7.append(site)
                    else:
                        self.mismatches_GU_outside.append(site)
                else:
                    if real_site == 1:
                        self.mismatches_nonGU_1.append(site)
                    elif real_site >= 2 and real_site <= 7:
                        self.mismatches_nonGU_2to7.append(site)
                    else:
                        self.mismatches_nonGU_outside.append(site)

        self.introducing_sites = []

class RuleRepertoire:

    @staticmethod
    def mismatch_rule(pair):
        mismatch_num_GU = len(pair.mismatches_GU_2to7) + len(pair.mismatches_GU_outside)
        mismatch_num_nonGU = len(pair.mismatches_nonGU_2to7) + len(pair.mismatches_nonGU_outside)
        return (mismatch_num_GU*0.9 + mismatch_num_nonGU) < 3

    @staticmethod
    def seeding_rule(pair):
        return (len(pair.mismatches_GU_2to7)*0.9 + len(pair.mismatches_nonGU_2to7)) < 1

    @staticmethod
    def nearly_satisfying_mismatch_rule(pair):
        mismatch_num_GU = len(pair.mismatches_GU_2to7) + len(pair.mismatches_GU_outside)
        mismatch_num_nonGU = len(pair.mismatches_nonGU_2to7) + len(pair.mismatches_nonGU_outside)
        if int(mismatch_num_GU*0.9 + mismatch_num_nonGU) == 3:
            if mismatch_num_GU == 1:
                pair.introducing_sites += pair.mismatches_nonGU_2to7 + pair.mismatches_nonGU_outside
            elif mismatch_num_GU > 1:
                pair.introducing_sites += pair.mismatches_nonGU_2to7 + pair.mismatches_nonGU_outside + pair.mismatches_GU_2to7 + pair.mismatches_GU_outside
        return len(pair.introducing_sites) > 0

    @staticmethod
    def nearly_satisfying_seeding_rule(pair):
        all_2to7 = pair.mismatches_GU_2to7 + pair.mismatches_nonGU_2to7
        if len(pair.mismatches_GU_2to7) == 2 and len(pair.mismatches_nonGU_2to7) == 0:
            pair.introducing_sites += pair.mismatches_GU_2to7
        elif len(pair.mismatches_GU_2to7) <= 1 and len(pair.mismatches_nonGU_2to7) == 1:
            pair.introducing_sites += pair.mismatches_nonGU_2to7
        return len(pair.introducing_sites) > 0


class PirnaPssParser:

    def __init__(self, bam_file, pirna_file, target_rna_fa):
        print('Parser loading')
        self.bam_file = bam_file
        self.pirna_file = pirna_file
        self.target_rna_dict = SeqIO.to_dict(SeqIO.parse(target_rna_fa, 'fasta'))
        print('OK')

    def __make_pss_record(self, target_host_name, start, end, strand):
        host_record = self.target_rna_dict[target_host_name]
        if end <= len(host_record):
            pss_record = SeqRecord(host_record[start:end].seq)
            pss_record.id = "{}#:#{}-{}".format(host_record.id, start, end)
            pss_record.name = ''
            pss_record.description = ''
            return pss_record
        else:
            raise BeyondTheBoundsException()

    def parse(self, rule):
        print(inspect.getsource(rule))
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam_aln, gzip.open(self.pirna_file, 'rt') as pirna_handle:
            pirna_records = SeqIO.parse(pirna_handle, 'fastq') if 'fastq' in self.pirna_file else SeqIO.parse(pirna_handle, 'fasta')
            for pirna_record, segment in zip(pirna_records, bam_aln):
                if pirna_record.id != segment.query_name:
                    raise Exception('piRNA names differ: {}/{}'.format(pirna_record.id, segment.query_name))
                if segment.flag not in [0, 4, 16, 20]:
                    raise Exception("Unexpected bam flag: {}.".format(segment.flag))
                if segment.flag == 0 or segment.flag == 16:
                    target_host_name = bam_aln.get_reference_name(segment.reference_id)
                    start = segment.reference_start                                         # The ordinate is transformed to be 0-based by pysam
                    end = start + segment.query_length
                    strand = '-' if segment.flag == 16 else '+'                   
                    mismatch_num = int(segment.get_tag('NM'))
                    try:
                        pss_record = self.__make_pss_record(target_host_name, start, end, strand)
                        pair = PirnaPssPair(pirna_record, pss_record, strand)
                        if(mismatch_num != len(pair.all_mismatches)):
                            if 'N' in str(pair.pss.seq).upper():
                                pass
                            else:
                                print(pair.pirna)
                                print(pair.pss)
                                print(pair.strand)
                                raise Exception("Main diff: {} vs mismatch_num: {}".format(len(pair.all_mismatches), mismatch_num))
                        else:
                            if rule(pair):
                                yield pair
                    except BeyondTheBoundsException:
                        pass

                if segment.has_tag('XA'):
                    xa_info = segment.get_tag('XA')
                    for xa_target in xa_info.split(';'):
                        if xa_target:
                            xa_target_split = xa_target.split(',')
                            target_host_name = xa_target_split[0]
                            start = int(xa_target_split[1][1:]) - 1                    # The ordinate needs to be transformed manually
                            end = start + segment.query_length
                            strand = xa_target_split[1][0]
                            mismatch_num = int(xa_target_split[3])
                            try:
                                pss_record = self.__make_pss_record(target_host_name, start, end, strand)
                                pair = PirnaPssPair(pirna_record, pss_record, strand)                            
                                if(mismatch_num != len(pair.all_mismatches)):
                                    if 'N' in str(pair.pss.seq).upper():
                                        pass
                                    else:
                                        print(pair.pirna)
                                        print(pair.pss)
                                        raise Exception("XA diff: {} vs mismatch_num: {}".format(len(pair.all_mismatches), mismatch_num))
                                else:
                                    if rule(pair):
                                        yield pair
                            except BeyondTheBoundsException:
                                pass
            print('OK')


class PirnaOutputer:

    def __init__(self, pirna_pos_in_genome, target_rna_info):
        print('Outputer loading')
        self.target_info_dict = {}
        with open(target_rna_info) as target_info_handle:
            for line in target_info_handle:
                line_split = line.strip().split('\t')
                self.target_info_dict[line_split[0].split('|')[0]] = line_split

        self.pirna_info_dict = defaultdict(list)
        with open(pirna_pos_in_genome) as in_genome_handle:
            head = True
            for line in in_genome_handle:
                if head:
                    head = False
                    continue
                line_split = line.strip().split('\t')
                self.pirna_info_dict[line_split[0]].append(line_split)
        print('OK')

    def __pss_genome_loc(self, target_rna_id, start_in_spliced_rna, end_in_spliced_rna):
        chrom_index, block_geno_start_index, block_geno_end_index = 5, [10,12,14],[11,13,15]
        if start_in_spliced_rna > end_in_spliced_rna:
            return 'null'
        gene_info = self.target_info_dict[target_rna_id]
        chromosome = gene_info[chrom_index] if 'chr' in gene_info[chrom_index] else 'chr'+ gene_info[chrom_index]

        split_sites = lambda s: [int(i) for i in s.split(';')]
        gene_to_genome_strand = int(gene_info[-1])

        if isinstance(block_geno_start_index, int) and isinstance(block_geno_end_index, int):
            block_geno_start_list = sorted(split_sites(gene_info[block_geno_start_index]))
            block_geno_end_list = sorted(split_sites(gene_info[block_geno_end_index]))
            genome_bounds_list = list(zip(block_geno_start_list, block_geno_end_list))    # 1-based
        elif isinstance(block_geno_start_index, list) and isinstance(block_geno_end_index, list):
            genome_bounds_list = []
            for each_start_index, each_end_index in zip(block_geno_start_index, block_geno_end_index):
                if gene_info[each_start_index] and gene_info[each_end_index]:
                    block_geno_start_list = sorted(split_sites(gene_info[each_start_index]))
                    block_geno_end_list = sorted(split_sites(gene_info[each_end_index]))
                    genome_bounds_list += list(zip(block_geno_start_list, block_geno_end_list))    # 1-based    

        if gene_to_genome_strand == 1:
            genome_bounds_list.sort(key=lambda i:i[0])
        else:
            genome_bounds_list.sort(key=lambda i:i[1], reverse=True)

        last_block_end = None
        block_bounds_list = []
        for block_geno_start, block_geno_end in genome_bounds_list:
            block_start = 1 if last_block_end == None else last_block_end + 1
            block_end = block_start + block_geno_end - block_geno_start
            block_bounds_list.append([block_start, block_end])
            last_block_end = block_end

        if end_in_spliced_rna > block_bounds_list[-1][1]:
            return 'null'

        if start_in_spliced_rna < block_bounds_list[0][0]:
            return 'null'

        for start_rank, (block_start, block_end) in enumerate(block_bounds_list):
            if start_in_spliced_rna >= block_start and start_in_spliced_rna <= block_end:
                break

        start_to_block_start = start_in_spliced_rna - block_start
        block_genome_start, block_genome_end = genome_bounds_list[start_rank]
        if gene_to_genome_strand == 1:
            start_in_genome = block_genome_start + start_to_block_start #1-based
        else:
            start_in_genome = block_genome_end - start_to_block_start

        for end_rank, (block_start, block_end) in enumerate(block_bounds_list):
            if end_in_spliced_rna >= block_start and end_in_spliced_rna <= block_end:
                break

        end_to_block_start = end_in_spliced_rna - block_start
        block_genome_start, block_genome_end = genome_bounds_list[end_rank]
        if gene_to_genome_strand == 1:
            end_in_genome = block_genome_start + end_to_block_start #1-based
        else:
            end_in_genome = block_genome_end - end_to_block_start

        if start_rank == end_rank:
            if gene_to_genome_strand == 1:
                segment_genome_loc = '{}:{}-{}'.format(chromosome, start_in_genome, end_in_genome)
            else:
                segment_genome_loc = '{}:{}-{}'.format(chromosome, end_in_genome, start_in_genome)    # end_in_genome < start_in_genome
        elif start_rank < end_rank:
            block_genome_start, block_genome_end = genome_bounds_list[start_rank]
            if gene_to_genome_strand == 1:
                segment_genome_loc_parts = ['{}:{}-{}'.format(chromosome, start_in_genome, block_genome_end)]
            else:
                segment_genome_loc_parts = ['{}:{}-{}'.format(chromosome, block_genome_start, start_in_genome)]
            for i in range(start_rank+1, end_rank+1):
                block_genome_start, block_genome_end = genome_bounds_list[i]                                
                if end_in_genome >= block_genome_start and end_in_genome <= block_genome_end:
                    if gene_to_genome_strand == 1:
                        segment_genome_loc_parts.append('{}:{}-{}'.format(chromosome, block_genome_start, end_in_genome))
                    else:
                        segment_genome_loc_parts.append('{}:{}-{}'.format(chromosome, end_in_genome, block_genome_end))
                else:
                    segment_genome_loc_parts.append('{}:{}-{}'.format(chromosome, block_genome_start, block_genome_end))

            if gene_to_genome_strand == -1:
                segment_genome_loc_parts.reverse()
            segment_genome_loc = '~'.join(segment_genome_loc_parts)
        else:
            print(gene_info)
            print(start_rank, end_rank)
            print(start_in_spliced_rna, end_in_spliced_rna)
            raise Exception('start_rank > end_rank')
        return segment_genome_loc+'({})'.format(gene_to_genome_strand)

    def __pirna_genome_loc(self, pirna_id):
        if pirna_id in self.pirna_info_dict:
            pirna_genome_loc_set = set()
            for pirna_info in self.pirna_info_dict[pirna_id]:
                pirna_candidate_loc = pirna_info[1]
                pirna_to_genome_strand = pirna_info[-1]
                pirna_genome_loc_set.add(pirna_candidate_loc+'({})'.format(pirna_to_genome_strand))
            pirna_genome_loc = ';'.join(pirna_genome_loc_set)
        else:
            pirna_genome_loc = 'null'
        return pirna_genome_loc


    def output(self, pairs, table_file):
        out_dir = os.path.dirname(table_file)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        with open(table_file, 'w') as table_handle:
            table_handle.write('\t'.join(['row_names', 'piRNA', 'gene', 'start', 'end', 'piRNA_genome_loc', 'PSS_genome_loc', 'strand', 'num_of_mismatches', 'mismatch_sites', 'introducing_sites', 'piRNA_sequence', 'PSS_sequence']))
            table_handle.write('\n')
            for pair in pairs:
                pirna_id = pair.pirna.id
                target_info = pair.pss.id.split('#:#')
                target_gene = target_info[0]
                target_gene_id = target_gene.split('|')[0]
                start = int(target_info[1].split('-')[0])  # 0-based
                end = int(target_info[1].split('-')[1])
                strand = pair.strand
                row_names = '.'.join([target_gene, str(start+1), str(end), strand, pirna_id]) # output 1-based site
                mismatch_sites = ','.join([str(i) for i in pair.all_mismatches])
                num_of_mismatches = str(len(pair.all_mismatches))
                introducing_sites = ','.join([str(i) for i in sorted(pair.introducing_sites)]) if pair.introducing_sites else '.'
                pirna_sequence = str(pair.pirna.seq)
                pss_sequence = str(pair.pss.seq)

                try:
                    pss_genome_loc = self.__pss_genome_loc(target_gene_id, start+1, end)
                except BeyondTheBoundsException:
                    print(pirna_id, target_gene, start, end)
                    print(pirna_sequence)
                    print(pss_sequence)
                    continue

                pirna_genome_loc = self.__pirna_genome_loc(pirna_id)

                table_handle.write('\t'.join([row_names, pirna_id, target_gene, str(start+1), str(end), pirna_genome_loc, pss_genome_loc, strand, num_of_mismatches, mismatch_sites, introducing_sites, pirna_sequence, pss_sequence]))
                table_handle.write('\n')


    def output_sam(self, pairs, sam_file, bam_template):
        with pysam.AlignmentFile(bam_template, 'rb') as bam_aln:
            with pysam.AlignmentFile(sam_file, mode='w', template=bam_aln) as sam_handle:
                for pair in pairs:
                    segment = pysam.AlignedSegment(bam_aln.header)
                    segment.query_name = pair.pirna.id
                    segment.query_sequence = str(pair.pirna.seq.reverse_complement())
                    segment.flag = 16
                    target_split = pair.target.id.split('#:#')
                    segment.reference_name = target_split[0]
                    segment.reference_start = int(target_split[1].split('-')[0])
                    segment.mapping_quality = 0
                    segment.cigar = ((0,len(pair.pirna)),)
                    segment.tags = (('NM', len(pair.all_mismatches)),)
                    segment.next_reference_id = -1
                    segment.next_reference_start = -1
                    sam_handle.write(segment)

        raw_bam_file = sam_file.replace('.sam', '.bam')
        sorted_bam_file = sam_file.replace('.sam', '.sorted.bam')
        # if not os.path.isfile(raw_bam_file):
        p_samtools_view = subprocess.Popen(['samtools', 'view', '-bS', sam_file, '-o', raw_bam_file])
        p_samtools_view.communicate()

        # if not os.path.isfile(sorted_bam_file):
        p_samtools_sort = subprocess.Popen(['samtools', 'sort', raw_bam_file, '-o', sorted_bam_file])
        p_samtools_sort.communicate()

        # if not os.path.isfile(sorted_bam_file+'.bai'):
        p_samtools_index = subprocess.Popen(['samtools', 'index', sorted_bam_file])
        p_samtools_index.communicate()


def pirna_pss_obtaining(bam_file, pirna_file, pirna_pos_in_genome, target_rna_fa, target_rna_info, out_table, out_sam, rule):
    if os.path.isfile(out_table):
        print('Output exists:', out_table)
        return
    parser = PirnaPssParser(bam_file, pirna_file, target_rna_fa)
    outputer = PirnaOutputer(pirna_pos_in_genome, target_rna_info)

    parsed_pairs = parser.parse(rule)
    outputer.output(parsed_pairs, out_table)

    # parsed_pairs = parser.parse(rule)
    # parser.output_sam(parsed_pairs, out_sam, bam_file)


def simple_obtaining(bam_file, pirna_file, pirna_pos_in_genome, pcg_fa, pcg_info, out_dir):
    sense_rule = lambda p: p.strand == '+'
    antisense_rule = lambda p: p.strand == '-'

    sense_out_table = os.path.join(out_dir, 'sense.txt')
    antisense_out_table = os.path.join(out_dir, 'antisense.txt')
    sense_out_sam = os.path.join(out_dir, 'sense.sam')    
    antisense_out_sam = os.path.join(out_dir, 'antisense.sam')

    pirna_pss_obtaining(bam_file, pirna_file, pirna_pos_in_genome, pcg_fa, pcg_info, sense_out_table, sense_out_sam, sense_rule)
    pirna_pss_obtaining(bam_file, pirna_file, pirna_pos_in_genome, pcg_fa, pcg_info, antisense_out_table, antisense_out_sam, antisense_rule)


def complicated_obtaining(bam_file, pirna_file, pirna_pos_in_genome, pcg_fa, pcg_info, out_dir):
    sense_rule = lambda p: p.strand == '+'
    antisense_rule = lambda p: p.strand == '-'
    complicated_rule = lambda p: (RuleRepertoire.mismatch_rule(p) and RuleRepertoire.seeding_rule(p)) or \
                                (RuleRepertoire.nearly_satisfying_mismatch_rule(p) and RuleRepertoire.seeding_rule(p)) or \
                                (RuleRepertoire.mismatch_rule(p) and RuleRepertoire.nearly_satisfying_seeding_rule(p))
    complicated_sense_rule = lambda p: sense_rule(p) and complicated_rule(p)
    complicated_antisense_rule = lambda p: antisense_rule(p) and complicated_rule(p)

    compliated_sense_out_table = os.path.join(out_dir, 'sense.txt')
    compliated_antisense_out_table = os.path.join(out_dir, 'antisense.txt')
    compliated_sense_out_sam = os.path.join(out_dir, 'sense.sam')    
    compliated_antisense_out_sam = os.path.join(out_dir, 'antisense.sam')

    pirna_pss_info_obtaining(bam_file, pirna_file, pirna_pos_in_genome, pcg_fa, pcg_info, compliated_sense_out_table, compliated_sense_out_sam, complicated_sense_rule)
    pirna_pss_info_obtaining(bam_file, pirna_file, pirna_pos_in_genome, pcg_fa, pcg_info, compliated_antisense_out_table, compliated_antisense_out_sam, complicated_antisense_rule)


def compile_human_simple():
    bam_file = '../generated_data/pirna_pss_bwa/human/mm4.bam'
    pirna_fq = '../generated_data/cleaned_pirnas/human_pirnas.healthy.dedup.retained.fastq.gz'
    pirna_pos_in_genome = '../generated_data/pirna_pos/human_healthy/human_pirnas.healthy.dedup.retained.pirna_pos.txt'
    pcg_fa = '../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'
    pcg_info = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    out_dir = '../generated_data/pirna_pss_pairs/human_simple'
    simple_obtaining(bam_file, pirna_fq, pirna_pos_in_genome, pcg_fa, pcg_info, out_dir)


def compile_human_complicated():
    bam_file = '../generated_data/pirna_pss_bwa/human/mm5.bam'
    pirna_fq = '../generated_data/cleaned_pirnas/human_pirnas.healthy.dedup.retained.fastq.gz'
    pirna_pos_in_genome = '../generated_data/pirna_pos/human_healthy/human_pirnas.healthy.dedup.retained.pirna_pos.txt'
    pcg_fa = '../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'
    pcg_info = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    out_dir = '../generated_data/pirna_pss_pairs/human_complicated'
    complicated_obtaining(bam_file, pirna_fq, pirna_pos_in_genome, pcg_fa, pcg_info, out_dir)


def compile_d_melanogaster():
    bam_file = '../generated_data/pirna_pss_bwa/d_melanogaster/mm4.bam'
    pirna_fa = '../original_data//pirna_sequences/dme.v3.0.fa.gz'
    pirna_pos_in_genome = '../generated_data/pirna_pos/d_melanogaster/dme.v3.0.pirna_pos.txt'
    pcg_fa = '../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna.info_seperated.fa'
    pcg_info = '../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna_info.txt'
    out_dir = '../generated_data/pirna_pss_pairs/d_melanogaster'
    simple_obtaining(bam_file, pirna_fa, pirna_pos_in_genome, pcg_fa, pcg_info, out_dir)


if __name__ == '__main__':
    compile_human_simple()
    compile_human_complicated()
    compile_d_melanogaster()
