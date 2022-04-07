from cyvcf2 import VCF
import os, glob
import random
import numpy as np
from collections import OrderedDict, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing 


human_reference_chromosomes = set([str(i+1) for i in range(22)] + ['X'])
d_melanogaster_reference_chromosomes = set(['2L', '2R', '3L', '3R', 'chr4', 'chrX', 'chrY']) 

rc_dict = {
    'a':'t',
    'A':'T',
    't':'a',
    'T':'A',
    'c':'g',
    'C':'G',
    'g':'c',
    'G':'C',
    }


def reverse_complement(nt):
    return rc_dict[nt]


def build_gene_info_dict(gene_info_file):
    gene_info_dict = {}
    with open(gene_info_file) as gene_info_handle:
        for line in gene_info_handle:
            line_split = line.strip().split('\t')
            # self.target_info_dict[line_split[0]] = line_split
            gene_info_dict[line_split[0].split('|')[0]] = line_split
    # print(gene_info_dict.keys())
    return gene_info_dict


def split_sites(sites_str):
    return [int(i) for i in sites_str.split(';')]

def categorize(alt, site_in_cdna, gene_info, cdna_seq, coding_start_index, coding_end_index):    # site_in_cdna is 1-based
    coding_bounds_in_cdna = list(zip(sorted(split_sites(gene_info[coding_start_index])), sorted(split_sites(gene_info[coding_end_index]))))    # 1-based
    coding_start, coding_end = coding_bounds_in_cdna[0][0], coding_bounds_in_cdna[-1][-1]
    try:
        transcript_strand = int(gene_info[STRAND])
    except:
        transcript_strand = 1 if gene_info[STRAND] == '+' else -1
    if site_in_cdna < coding_start:
        return '5\'UTR'
    elif site_in_cdna > coding_end:
        return '3\'UTR'
    else: 
        site_in_coding_region = site_in_cdna - coding_start    # 0-based
        codon_start = int(site_in_coding_region/3) * 3 + coding_start - 1
        codon_pos = site_in_coding_region % 3
        ref_codon = cdna_seq[codon_start:codon_start+3]
        alt_codon = cdna_seq[codon_start:codon_start+codon_pos] + Seq(alt) + cdna_seq[codon_start+codon_pos+1:codon_start+3] if transcript_strand == 1 else cdna_seq[codon_start:codon_start+codon_pos] + Seq(reverse_complement(alt)) + cdna_seq[codon_start+codon_pos+1:codon_start+3] 
        # if 'ENSG00000185960' in gene_info[0]:
        #     print(gene_info) 
        #     print(site_in_cdna)
        #     print(codon_start)
        #     print(cdna_seq[1192:1195])
        #     print(ref_codon, alt_codon)
        if ref_codon.translate() == alt_codon.translate():
            return 'Synonymous'
        else:
            return 'Non-synonymous'

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


class Stat:
    def __init__(self):
        self.magnetizing_var_count = 0
        self.non_magnetizing_var_count = 0
        self.magnetizing_var_dict = defaultdict(set)
        self.all_var_dict = defaultdict(set)
        self.magnetizing_var_list = []
        self.all_var_list = []
        self.var_info_dict = {}

    def update(self, var, alt, var_info, magnetizing):
        var_str = '{}:{}-{}'.format(var.CHROM, var.start+1, var.start+1)
        self.all_var_dict[var.CHROM].add((var_str,alt))
        if (var_str,alt) not in self.var_info_dict:
            self.var_info_dict[(var_str,alt)] = var_info

        if magnetizing:
            self.magnetizing_var_count += 1        
            self.magnetizing_var_dict[var.CHROM].add((var_str,alt))
        else:
            self.non_magnetizing_var_count += 1

    def show(self, magnetizing_var_info_file, non_magnetizing_var_info_file, stat_file):
        out_dir = os.path.dirname(magnetizing_var_info_file)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        var_sum = self.magnetizing_var_count + self.non_magnetizing_var_count
        with open(stat_file, 'w') as stat_handle:
            stat_handle.write('{}\t{}\t{}\n'.format(self.magnetizing_var_count, self.non_magnetizing_var_count, var_sum))
            stat_handle.write('{}\t{}\n'.format(self.magnetizing_var_count/var_sum, self.non_magnetizing_var_count/var_sum))

        non_magnetizing_var_dict = defaultdict(set)
        for key in self.all_var_dict.keys():
            non_magnetizing_var_dict[key] = self.all_var_dict[key].difference(self.magnetizing_var_dict[key])

        magnetizing_allele_list = []
        magnetizing_var_info_list = []       
        for chrom, magnetizing_var_set in self.magnetizing_var_dict.items():
            for item in magnetizing_var_set:
                var_info = self.var_info_dict[item]
                magnetizing_allele_list.append(item)
                magnetizing_var_info_list.append(var_info)

        non_magnetizing_allele_list = []
        non_magnetizing_var_info_list = []
        for chrom, non_magnetizing_var_set in non_magnetizing_var_dict.items():
            for item in non_magnetizing_var_set:
                var_info = self.var_info_dict[item]
                non_magnetizing_allele_list.append(item)
                non_magnetizing_var_info_list.append(var_info)

        head_line = '\t'.join(['LOC(1-based)', 'ALT']+list(var_info.keys()))
        if magnetizing_var_info_file:
            with open(magnetizing_var_info_file, 'w') as magnetizing_var_info_handle:
                magnetizing_var_info_handle.write(head_line)
                magnetizing_var_info_handle.write('\n')
                for (var_str, alt), var_info in zip(magnetizing_allele_list, magnetizing_var_info_list):
                    out_info = [var_str, alt] + [str(i) for i in var_info.values()]
                    magnetizing_var_info_handle.write('\t'.join(out_info))
                    magnetizing_var_info_handle.write('\n')

        if non_magnetizing_var_info_file:
            with open(non_magnetizing_var_info_file, 'w') as non_magnetizing_var_info_handle:
                non_magnetizing_var_info_handle.write(head_line)
                non_magnetizing_var_info_handle.write('\n')
                for (var_str, alt), var_info in zip(non_magnetizing_allele_list, non_magnetizing_var_info_list):
                    out_info = [var_str, alt] + [str(i) for i in var_info.values()]
                    non_magnetizing_var_info_handle.write('\t'.join(out_info))
                    non_magnetizing_var_info_handle.write('\n')


def get_var_corresponding_to_pirna_sites(each_genome_loc, pirna_sites, vcf_data, pirna_to_loc_transcript_strand, reference_chromosomes): 
    # e.g. chr1:123456-123478~chr1:123567-123576(-1) or chr1:123456-123488(1)
    pre_var_list = []
    var_list = []
    # print(each_genome_loc)
    loc_transcript_strand = int(each_genome_loc.split('(')[1].strip(')'))              # relative strand relationship between the genomic sequence and the mRNA or piRNA
    pirna_to_genome_strand = loc_transcript_strand * pirna_to_loc_transcript_strand    # relative strand relationship between the genomic sequence and the piRNA

    buoy_in_pirna = 1
    loc_segments = each_genome_loc.split('(')[0].split('~') if pirna_to_genome_strand == 1 else reversed(each_genome_loc.split('(')[0].split('~')) 
    for loc_segment in loc_segments:
        chromosome = loc_segment.split(':')[0]

        if chromosome.strip('chr') not in reference_chromosomes:
            # print(genome_loc)
            continue

        if isinstance(vcf_data, dict):
            if chromosome in vcf_data:
                vcf_obj = vcf_data[chromosome]
            else:
                vcf_obj = None
        else:
            vcf_obj = vcf_data
        if vcf_obj == None:
            continue

        if pirna_to_genome_strand == 1:            
            segment_start = int(loc_segment.split(':')[1].split('-')[0])               # 1-based; the buoy corresponds to segment_start
            segment_end = int(loc_segment.split(':')[1].split('-')[1])                 # 1-based
            for variant in vcf_obj(loc_segment):
                if ('VT' in variant.INFO and variant.INFO['VT'] == 'SNP') or (len(variant.REF) == 1 and len(variant.ALT) == 1 and len(variant.ALT[0]) == 1 and variant.REF.upper() in ['A', 'T', 'G', 'C'] and variant.ALT[0].upper() in ['A', 'T', 'G', 'C']):
                    dist_to_segment_start = (variant.start + 1) - segment_start        # variant.start is 0-based; counting from the forward direction
                    var_site_in_pirna = buoy_in_pirna + dist_to_segment_start          # 1-based
                    if var_site_in_pirna in pirna_sites:
                        pre_var_list.append((var_site_in_pirna, variant, loc_transcript_strand))
            buoy_in_pirna += segment_end - segment_start + 1                           # the buoy moves to the piRNA site that corresponds to next segment
        else:
            segment_end = int(loc_segment.split(':')[1].split('-')[0])                 # 1-based
            segment_start = int(loc_segment.split(':')[1].split('-')[1])               # 1-based; the buoy corresponds to segment_start
            for variant in reversed(list(vcf_obj(loc_segment))):
                if ('VT' in variant.INFO and variant.INFO['VT'] == 'SNP') or (len(variant.REF) == 1 and len(variant.ALT) == 1 and len(variant.ALT[0]) == 1 and variant.REF.upper() in ['A', 'T', 'G', 'C'] and variant.ALT[0].upper() in ['A', 'T', 'G', 'C']):
                    dist_to_segment_start = segment_start - (variant.start + 1)        # variant.start is 0-based; counting from the backward direction 
                    var_site_in_pirna = buoy_in_pirna + dist_to_segment_start          # 1-based
                    if var_site_in_pirna in pirna_sites:
                        pre_var_list.append((var_site_in_pirna, variant, loc_transcript_strand))
            buoy_in_pirna += segment_start - segment_end + 1                           # the buoy moves to the piRNA site that corresponds to next segment

    # The following lines are set to build a list in which None means corresponds to no SNP
    pirna_len = buoy_in_pirna
    if not pre_var_list:
        var_list += [None] * pirna_len
    else:
        var_list = []
        last_var_site = 0
        for var_site_in_pirna, variant, loc_transcript_strand in pre_var_list:
            # if last_var_site = x and var_site_in_pirna = y, then not include last_var_site but include var_site_in_pirna, there are y-x nucleotides
            # last_var_site and var_site_in_pirna should not be "None", thus add y-x-1 "None" to the list
            var_list += [None] * (var_site_in_pirna - last_var_site - 1)
            var_list.append((var_site_in_pirna, variant, loc_transcript_strand))
            last_var_site = var_site_in_pirna
        var_list += [None] * (pirna_len - var_site_in_pirna)
    return var_list 


def build_vcf_dict(vcf_dir):
    vcf_dict = {}
    for vcf_file in glob.glob(os.path.join(vcf_dir, '*.vcf.gz')):
        # key = os.path.basename(vcf_file).split('.')[1]
        key = os.path.basename(vcf_file).split('.')[0].split('_')[-1]
        vcf_dict[key] = VCF(vcf_file)
    return vcf_dict


PIRNA_ID = 1
CDNA_ID = 2
START_IN_CDNA = 3
END_IN_CDNA = 4
LOC_OF_PIRNA = 5
LOC_OF_PSS = 6
MISMATCH_SITES = 9
HYBRID_TAR_SEQ = -2
HYBRID_PIRNA_SEQ = -3

GENO_NAME = 5
GENO_START = 6
GENO_END = 7
STRAND = 9
REP_NAME = 10
REP_CLASS = 11
REP_FAMILY = 12
REP_START = 13
REP_END = 14

PIRNA_SEQ = -2
PSS_SEQ = -1

def magnetizing_alleles_in_pirnas(pirna_pss_pair_file, vcf_input, reference_chromosomes, pirna_to_pss_strand, magnetizing_var_info_file, non_magnetizing_var_info_file, stat_file):
    print('magnetizing_alleles_in_pirnas: '+pirna_pss_pair_file)
    vcf_data = build_vcf_dict(vcf_input) if os.path.isdir(vcf_input) else VCF(vcf_input)
    stat = Stat()
    len_total = 0
    pirna_pss_pair_count = 0
    count = 0

    with open(pirna_pss_pair_file) as pirna_pss_pair_handle:
        head = True
        for line in pirna_pss_pair_handle:
            if head:
                head = False
                continue
            if count % 1000000 == 0:
                print('Analyzed: {}'.format(count), stat_file)
            count += 1
            pirna_pss_pair_info = line.strip().split('\t')
            pirna_seq = pirna_pss_pair_info[PIRNA_SEQ]
            pss_seq = pirna_pss_pair_info[PSS_SEQ]
            
            pss_loc = pirna_pss_pair_info[LOC_OF_PSS]
            pirna_loc = pirna_pss_pair_info[LOC_OF_PIRNA]
            if pirna_loc == 'null':
                continue

            pirna_loc_list = pirna_loc.split(';')

            overlapping = False
            for each_pirna_loc in pirna_loc_list:
                if pirna_pss_overlapping(each_pirna_loc, pss_loc):
                    overlapping = True
                    break    
            if overlapping:
                continue    # if the PSS is overlapped with any possible position of the piRNA, the piRNA-PSS pair will be excluded.

            if pirna_pss_pair_info[MISMATCH_SITES]: 
                mismatch_sites = [int(site) for site in pirna_pss_pair_info[MISMATCH_SITES].split(',')]
            else:
                continue

            pirna_seq_len = len(pirna_seq)
            len_total += pirna_seq_len
            pirna_pss_pair_count += 1 
            for each_pirna_loc in pirna_loc_list:
                var_list = get_var_corresponding_to_pirna_sites(each_pirna_loc, mismatch_sites, vcf_data, 1, reference_chromosomes)   # pirna_to_loc_transcript_strand = 1
                for item in var_list:
                    if item != None:
                        var_site_in_pirna, variant, loc_transcript_strand = item    #var_site_in_pirna is 1-based
                        alt = variant.ALT[0]
                        var_info = OrderedDict()
                        var_info['piRNA_transcript_strand'] = loc_transcript_strand
                        try:
                            var_info['Allele_frequency'] = variant.INFO['AF']
                        except KeyError:
                            pass
                        try:
                            var_info['ALTCOUNT'] = variant.INFO['ALTCOUNT']
                        except KeyError:
                            pass
                        try:
                            var_info['REFCOUNT'] = variant.INFO['REFCOUNT']
                        except KeyError:
                            pass                        

                        if pirna_to_pss_strand == 1:
                            pss_nt = pss_seq[var_site_in_pirna-1]
                            if loc_transcript_strand == 1:
                                stat.update(variant, alt, var_info, pss_nt.upper() == alt.upper())
                                # print(var_site_in_pirna)
                                # print(variant.CHROM, variant.start, variant.REF)
                                # print(each_pirna_loc)
                            else:
                                stat.update(variant, alt, var_info, pss_nt.upper() == reverse_complement(alt.upper()))
                        else:
                            pss_nt = pss_seq[-var_site_in_pirna]
                            if loc_transcript_strand == 1:
                                stat.update(variant, alt, var_info, pss_nt.upper() == reverse_complement(alt.upper()))
                            else:
                                stat.update(variant, alt, var_info, pss_nt.upper() == alt.upper())

        avg_len_total = len_total/pirna_pss_pair_count 
        print(len_total, pirna_pss_pair_count, avg_len_total)
        stat.show(magnetizing_var_info_file, non_magnetizing_var_info_file, stat_file)
    print('OK')


def magnetizing_alleles_in_pss(pirna_pss_pair_file, vcf_input, reference_chromosomes, pirna_to_pss_strand, cdna_file, cdna_info_file, magnetizing_var_info_file, non_magnetizing_var_info_file, stat_file):  
    print('magnetizing_alleles_in_pss: '+pirna_pss_pair_file)
    coding_start_index, coding_end_index = 8, 9
    vcf_data = build_vcf_dict(vcf_input) if os.path.isdir(vcf_input) else VCF(vcf_input)
    cdna_dict = SeqIO.to_dict(SeqIO.parse(cdna_file, 'fasta'), key_function=lambda i:i.id.split('|')[0])
    cdna_info_dict = build_gene_info_dict(cdna_info_file)
    stat = Stat()
    len_total = 0
    pirna_pss_pair_count = 0
    count = 0

    with open(pirna_pss_pair_file) as pirna_pss_pair_handle:
        head=True
        for line in pirna_pss_pair_handle:
            if head:
                head = False
                continue
            if count % 1000000 == 0:
                print('Analyzed: {}'.format(count), stat_file)
            count += 1
            pirna_pss_pair_info = line.strip().split('\t')
            pirna_seq = pirna_pss_pair_info[PIRNA_SEQ]

            gene_id = pirna_pss_pair_info[CDNA_ID].split('|')[0]
            cdna_seq = cdna_dict[gene_id].seq
            cdna_info = cdna_info_dict[gene_id]

            pss_loc = pirna_pss_pair_info[LOC_OF_PSS]
            pirna_loc = pirna_pss_pair_info[LOC_OF_PIRNA]
            if pss_loc == 'null':
                raise Exception('PSS location cannot be null.')
            
            pirna_loc_list = pirna_loc.split(';')
            
            overlapping = False
            for each_pirna_loc in pirna_loc_list:
                if pirna_pss_overlapping(each_pirna_loc, pss_loc):
                    overlapping = True
                    break    
            if overlapping:
                continue    # if the PSS is overlapped with any possible position of the piRNA, the piRNA-PSS pair will be excluded.
            
            if pirna_pss_pair_info[MISMATCH_SITES]:                
                mismatch_sites = [int(site) for site in pirna_pss_pair_info[MISMATCH_SITES].split(',')]
            else:
                continue
            
            pirna_seq_len = len(pirna_seq)
            len_total += pirna_seq_len
            pirna_pss_pair_count += 1
            var_list = get_var_corresponding_to_pirna_sites(pss_loc, mismatch_sites, vcf_data, pirna_to_pss_strand, reference_chromosomes)
            for item in var_list:
                if item != None:
                    var_site_in_pirna, variant, loc_transcript_strand = item    #var_site_in_pirna is 1-based
                    tar_pirna_nt = pirna_seq[var_site_in_pirna-1]
                    alt = variant.ALT[0]
                    var_site_in_cdna = int(pirna_pss_pair_info[START_IN_CDNA]) + var_site_in_pirna - 1 if pirna_to_pss_strand == 1 else int(pirna_pss_pair_info[END_IN_CDNA]) - (var_site_in_pirna - 1)
            
                    var_info = OrderedDict()
                    var_info['Category'] = categorize(alt, var_site_in_cdna, cdna_info, cdna_seq, coding_start_index, coding_end_index)
                    try:
                        var_info['Allele_frequency'] = variant.INFO['AF']
                    except KeyError:
                        pass
                    try:
                        var_info['ALTCOUNT'] = variant.INFO['ALTCOUNT']
                    except KeyError:
                        pass
                    try:
                        var_info['REFCOUNT'] = variant.INFO['REFCOUNT']
                    except KeyError:
                        pass  
                        
                    if pirna_to_pss_strand == 1:
                        if loc_transcript_strand == 1:
                            stat.update(variant, alt, var_info, tar_pirna_nt.upper() == alt.upper())
                        else:
                            stat.update(variant, alt, var_info, tar_pirna_nt.upper() == reverse_complement(alt.upper()))
                    else:
                        if loc_transcript_strand == 1:
                            stat.update(variant, alt, var_info, tar_pirna_nt.upper() == reverse_complement(alt.upper()))
                        else:
                            stat.update(variant, alt, var_info, tar_pirna_nt.upper() == alt.upper())

        avg_len_total = len_total/pirna_pss_pair_count 
        print(len_total, pirna_pss_pair_count, avg_len_total)
        stat.show(magnetizing_var_info_file, non_magnetizing_var_info_file, stat_file)
    print('OK')


def shared_polymorphism(pirna_pss_pair_file, vcf_input, reference_chromosomes, pirna_to_pss_strand, out_file):
    print('shared_polymorphism: '+pirna_pss_pair_file)
    coding_start_index, coding_end_index = 8, 9
    rand_num = 30
    vcf_data = build_vcf_dict(vcf_input) if os.path.isdir(vcf_input) else VCF(vcf_input)

    len_total = 0
    pirna_pss_pair_count = 0

    shared_count = 0
    non_shared_count = 0
    query_snp_count = 0
    hit_snp_count = 0
    total = 0
    shuffled_shared_count_list = [0]*rand_num
    shuffled_non_shared_count_list = [0]*rand_num
    shuffled_query_snp_count_list = [0]*rand_num
    shuffled_hit_snp_count_list = [0]*rand_num    
    shuffled_total_list = [0]*rand_num

    with open(pirna_pss_pair_file) as pirna_pss_pair_handle:
        head=True
        for line in pirna_pss_pair_handle:
            if head:
                head = False
                continue
            pirna_pss_pair_info = line.strip().split('\t')
            pirna_seq = pirna_pss_pair_info[PIRNA_SEQ]
            pirna_seq_len = len(pirna_seq)
            len_total += pirna_seq_len
            pirna_pss_pair_count += 1
            
            gene_id = pirna_pss_pair_info[0].split('|')[0]

            pss_loc = pirna_pss_pair_info[LOC_OF_PSS]
            pirna_loc = pirna_pss_pair_info[LOC_OF_PIRNA]
            if pirna_loc == 'null':
                continue

            pirna_loc_list = pirna_loc.split(';')

            overlapping = False
            for each_pirna_loc in pirna_loc_list:
                if pirna_pss_overlapping(each_pirna_loc, pss_loc):
                    overlapping = True
                    break    
            if overlapping:
                continue    # if the PSS is overlapped with the piRNA, the piRNA-PSS pair will be excluded.

            for each_pirna_loc in pirna_loc_list:                
                if pirna_pss_pair_info[MISMATCH_SITES]:
                    mismatch_sites = [int(i) for i in pirna_pss_pair_info[MISMATCH_SITES].split(',')]
                    match_sites = []
                    for i in range(pirna_seq_len):
                        site = i + 1
                        if site not in mismatch_sites:
                            match_sites.append(site)
                else:
                    match_sites = list(range(pirna_seq_len))

                pirna_var_list = get_var_corresponding_to_pirna_sites(each_pirna_loc, match_sites, vcf_data, 1, reference_chromosomes)   # pirna_to_loc_transcript_strand = 1
                pss_var_list = get_var_corresponding_to_pirna_sites(pss_loc, match_sites, vcf_data, pirna_to_pss_strand, reference_chromosomes)

                for pirna_item, pss_item in zip(pirna_var_list, pss_var_list):    #var_site_in_pirna is 1-based
                    total += 1
                    if pirna_item != None or pss_item != None:
                        shared = False
                        if pirna_item != None and pss_item != None:
                            pirna_var_site_in_pirna, pirna_var, pirna_loc_transcript_strand = pirna_item
                            pss_var_site_in_pirna, pss_var, pss_loc_transcript_strand = pss_item
                            pirna_alt = pirna_var.ALT[0].upper()
                            pss_alt = pss_var.ALT[0].upper()
                            if pirna_to_pss_strand * pirna_loc_transcript_strand * pss_loc_transcript_strand == 1:
                                if pirna_alt == pss_alt:
                                    shared = True
                            else:
                                if reverse_complement(pirna_alt) == pss_alt:
                                    shared = True
                            # print(match_sites)
                            # print(pirna_var_list, pss_var_list)
                            # print(pirna_var_site_in_pirna, pss_var_site_in_pirna)
                            # print(pirna_alt, pss_alt)
                        if shared:
                            shared_count += 1
                        else:
                            non_shared_count += 1

                for k in range(rand_num):
                    random.shuffle(pirna_var_list)
                    random.shuffle(pss_var_list)
                    for pirna_item, pss_item in zip(pirna_var_list, pss_var_list):
                        shuffled_total_list[k] += 1
                        if pirna_item != None or pss_item != None:
                            shared = False
                            if pirna_item != None and pss_item != None:
                                pirna_var_site_in_pirna, pirna_var, pirna_loc_transcript_strand = pirna_item
                                pss_var_site_in_pirna, pss_var, pss_loc_transcript_strand = pss_item                        
                                pirna_alt = pirna_var.ALT[0].upper()
                                pss_alt = pss_var.ALT[0].upper()
                                if pirna_to_pss_strand * pirna_loc_transcript_strand * pss_loc_transcript_strand == 1:
                                    if pirna_alt == pss_alt:
                                        shared = True
                                else:
                                    if reverse_complement(pirna_alt) == pss_alt:
                                        shared = True
                            if shared:
                                shuffled_shared_count_list[k] += 1
                            else:
                                shuffled_non_shared_count_list[k] += 1
                # print(pirna_pss_pair_count)
                if pirna_pss_pair_count % 1000000 == 0:
                    print(out_file)
                    print(shared_count, non_shared_count, total)
                    print()
                    for k in range(rand_num):
                        print(shuffled_shared_count_list[k], shuffled_non_shared_count_list[k], shuffled_total_list[k])
                    print()

    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(out_file, 'w') as out_handle:
        out_handle.write('\t'.join(str(each) for each in (shared_count, non_shared_count, total)))
        out_handle.write('\n\n')
        for k in range(rand_num):
            out_handle.write('\t'.join(str(each) for each in (shuffled_shared_count_list[k], shuffled_non_shared_count_list[k], shuffled_total_list[k])))
            out_handle.write('\n')
    print('OK')


def analyze(vcf_data, cdna_info_file, cdna_file, pirna_pss_pairs_dir, var_dir, shared_polymorphism_dir, reference_chromosomes, processes=6):
    sense_pirna_pss_pair_file = os.path.join(pirna_pss_pairs_dir, 'sense.txt')
    antisense_pirna_pss_pair_file = os.path.join(pirna_pss_pairs_dir, 'antisense.txt')

    pirna_sense_magnetizing_var_info_file = os.path.join(var_dir, 'pirna.sense.magnetizing_var_info.txt')
    pirna_sense_non_magnetizing_var_info_file = os.path.join(var_dir, 'pirna.sense.non_magnetizing_var_info.txt')
    pirna_sense_stat_file = os.path.join(var_dir, 'pirna.sense.stat.txt')
    pirna_antisense_magnetizing_var_info_file = os.path.join(var_dir, 'pirna.antisense.magnetizing_var_info.txt')
    pirna_antisense_non_magnetizing_var_info_file = os.path.join(var_dir, 'pirna.antisense.non_magnetizing_var_info.txt')
    pirna_antisense_stat_file = os.path.join(var_dir, 'pirna.antisense.stat.txt')

    pss_sense_magnetizing_var_info_file = os.path.join(var_dir, 'pss.sense.magnetizing_var_info.txt')
    pss_sense_non_magnetizing_var_info_file = os.path.join(var_dir, 'pss.sense.non_magnetizing_var_info.txt')
    pss_sense_stat_file = os.path.join(var_dir, 'pss.sense.stat.txt')
    pss_antisense_magnetizing_var_info_file = os.path.join(var_dir, 'pss.antisense.magnetizing_var_info.txt')
    pss_antisense_non_magnetizing_var_info_file = os.path.join(var_dir, 'pss.antisense.non_magnetizing_var_info.txt')
    pss_antisense_stat_file = os.path.join(var_dir, 'pss.antisense.stat.txt')

    sense_shared_polymorphism_file = os.path.join(shared_polymorphism_dir, 'sense_shared_polymorphism.txt')
    antisense_shared_polymorphism_file = os.path.join(shared_polymorphism_dir, 'antisense_shared_polymorphism.txt')

    magnetizing_alleles_in_pirnas(sense_pirna_pss_pair_file, vcf_data, reference_chromosomes, 1, pirna_sense_magnetizing_var_info_file, pirna_sense_non_magnetizing_var_info_file, pirna_sense_stat_file)
    magnetizing_alleles_in_pss(sense_pirna_pss_pair_file, vcf_data, reference_chromosomes, 1, cdna_file, cdna_info_file, pss_sense_magnetizing_var_info_file, pss_sense_non_magnetizing_var_info_file, pss_sense_stat_file)
    magnetizing_alleles_in_pirnas(antisense_pirna_pss_pair_file, vcf_data, reference_chromosomes, -1, pirna_antisense_magnetizing_var_info_file, pirna_antisense_non_magnetizing_var_info_file, pirna_antisense_stat_file)
    magnetizing_alleles_in_pss(antisense_pirna_pss_pair_file, vcf_data, reference_chromosomes, -1, cdna_file, cdna_info_file, pss_antisense_magnetizing_var_info_file, pss_antisense_non_magnetizing_var_info_file, pss_antisense_stat_file)
    shared_polymorphism(sense_pirna_pss_pair_file, vcf_data, reference_chromosomes, 1, sense_shared_polymorphism_file)
    shared_polymorphism(antisense_pirna_pss_pair_file, vcf_data, reference_chromosomes, -1, antisense_shared_polymorphism_file)
    
    pool = multiprocessing.Pool(processes=processes)

    pool.apply_async(magnetizing_alleles_in_pirnas, (sense_pirna_pss_pair_file, vcf_data, reference_chromosomes, 1, pirna_sense_magnetizing_var_info_file, pirna_sense_non_magnetizing_var_info_file, pirna_sense_stat_file))
    pool.apply_async(magnetizing_alleles_in_pss, (sense_pirna_pss_pair_file, vcf_data, reference_chromosomes, 1, cdna_file, cdna_info_file, pss_sense_magnetizing_var_info_file, pss_sense_non_magnetizing_var_info_file, pss_sense_stat_file))
    pool.apply_async(magnetizing_alleles_in_pirnas, (antisense_pirna_pss_pair_file, vcf_data, reference_chromosomes, -1, pirna_antisense_magnetizing_var_info_file, pirna_antisense_non_magnetizing_var_info_file, pirna_antisense_stat_file))
    pool.apply_async(magnetizing_alleles_in_pss, (antisense_pirna_pss_pair_file, vcf_data, reference_chromosomes, -1, cdna_file, cdna_info_file, pss_antisense_magnetizing_var_info_file, pss_antisense_non_magnetizing_var_info_file, pss_antisense_stat_file))
    pool.apply_async(shared_polymorphism, (sense_pirna_pss_pair_file, vcf_data, reference_chromosomes, 1, sense_shared_polymorphism_file))
    pool.apply_async(shared_polymorphism, (antisense_pirna_pss_pair_file, vcf_data, reference_chromosomes, -1, antisense_shared_polymorphism_file))

    pool.close()
    pool.join()
    print('var_dir:', var_dir)


def analyze_1000_genomes_simple(): 
    vcf_data = '../original_data/polymorphism_data/1000_Genomes_high_cov'
    cdna_info_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    cdna_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'
    pirna_pss_pairs_dir = '../generated_data/pirna_pss_pairs/human_simple'
    var_dir = '../generated_data/var_in_pirna_pss_pairs/1000_Genomes_simple'
    shared_polymorphism_dir = '../generated_data/pirna_pss_shared_polymorphism/1000_Genomes_simple'
    analyze(vcf_data, cdna_info_file, cdna_file, pirna_pss_pairs_dir, var_dir, shared_polymorphism_dir, human_reference_chromosomes)


def analyze_uk10k_simple(): 
    vcf_data = '../original_data/polymorphism_data/UK10K/UK10K_COHORT.20160215.sites.chr_added.lifted.vcf.gz'
    cdna_info_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    cdna_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'
    pirna_pss_pairs_dir = '../generated_data/pirna_pss_pairs/human_simple'
    var_dir = '../generated_data/var_in_pirna_pss_pairs/UK10K_simple'
    shared_polymorphism_dir = '../generated_data/pirna_pss_shared_polymorphism/UK10K_simple'
    analyze(vcf_data, cdna_info_file, cdna_file, pirna_pss_pairs_dir, var_dir, shared_polymorphism_dir, human_reference_chromosomes)


def analyze_1000_genomes_complicated(): 
    vcf_data = '../original_data/polymorphism_data/1000_Genomes_high_cov'
    cdna_info_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    cdna_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'
    pirna_pss_pairs_dir = '../generated_data/pirna_pss_pairs/human_complicated'
    var_dir = '../generated_data/var_in_pirna_pss_pairs/1000_Genomes_complicated'
    shared_polymorphism_dir = '../generated_data/pirna_pss_shared_polymorphism/1000_Genomes_complicated'
    analyze(vcf_data, cdna_info_file, cdna_file, pirna_pss_pairs_dir, var_dir, shared_polymorphism_dir, human_reference_chromosomes)


def analyze_uk10k_complicated(): 
    vcf_data = '../original_data/polymorphism_data/UK10K/UK10K_COHORT.20160215.sites.chr_added.lifted.vcf.gz'
    cdna_info_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    cdna_file = '../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'
    pirna_pss_pairs_dir = '../generated_data/pirna_pss_pairs/human_complicated'
    var_dir = '../generated_data/var_in_pirna_pss_pairs/UK10K_complicated'
    shared_polymorphism_dir = '../generated_data/pirna_pss_shared_polymorphism/UK10K_complicated'
    analyze(vcf_data, cdna_info_file, cdna_file, pirna_pss_pairs_dir, var_dir, shared_polymorphism_dir, human_reference_chromosomes)


def analyze_dgrp2(): 
    vcf_data = '../original_data/polymorphism_data/D_melanogaster_Genetic_Reference_Panel_2/dgrp2_dm6.vcf.gz'
    cdna_info_file = r'../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna_info.txt'
    cdna_file = '../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna.info_seperated.fa'    
    pirna_pss_pairs_dir = '../generated_data/pirna_pss_pairs/d_melanogaster'
    var_dir = '../generated_data/var_in_pirna_pss_pairs/DGRP2'
    shared_polymorphism_dir = '../generated_data/pirna_pss_shared_polymorphism/DGRP2'
    analyze(vcf_data, cdna_info_file, cdna_file, pirna_pss_pairs_dir, var_dir, shared_polymorphism_dir, d_melanogaster_reference_chromosomes)



if __name__ == '__main__':
    analyze_1000_genomes_simple()
    analyze_uk10k_simple()    
    analyze_1000_genomes_complicated()
    analyze_uk10k_complicated()
    analyze_dgrp2()

