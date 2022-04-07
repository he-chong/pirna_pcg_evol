import os
from Bio import SeqIO
from collections import OrderedDict


def cdna_longest_seq(ensembl_fa, gene_info_file, seq_fa):
    out_dir = os.path.dirname(gene_info_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    out_gene_dict = OrderedDict()
    gene_info_dict = OrderedDict()
    for record in SeqIO.parse(ensembl_fa, 'fasta'):
        gene_info = record.id.split('|')
        gene_id = gene_info[0]
        if gene_id in out_gene_dict and len(record.seq) < len(out_gene_dict[gene_id].seq):  # Select the longest sequence
            continue
        gene_name = gene_info[2]
        record.id = gene_id + '|' + gene_name
        record.description = ''
        record.name = ''
        out_gene_dict[gene_id] = record
        gene_info_dict[gene_id] = gene_info

    with open(gene_info_file, 'w') as gene_info_handle, open(seq_fa, 'w') as seq_handle:
        for gene_id in out_gene_dict:
            gene_info = gene_info_dict[gene_id]
            record = out_gene_dict[gene_id]
            gene_info_handle.write('\t'.join(gene_info))
            gene_info_handle.write('\n')
            SeqIO.write(record, seq_handle, 'fasta')


if __name__ == '__main__':
    human_ensembl_fa = r'../original_data/protein_coding_genes/human_protein_coding_cdna.fa'
    human_gene_info_file = r'../generated_data/protein_coding_genes/human_protein_coding_cdna_info.txt'
    human_seq_fa = r'../generated_data/protein_coding_genes/human_protein_coding_cdna.info_seperated.fa'

    fly_ensembl_fa = r'../original_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna.fa'
    fly_gene_info_file = r'../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna_info.txt'
    fly_seq_fa = r'../generated_data/protein_coding_genes/drosophila_melanogaster_protein_coding_cdna.info_seperated.fa'

    cdna_longest_seq(human_ensembl_fa, human_gene_info_file, human_seq_fa)
    cdna_longest_seq(fly_ensembl_fa, fly_gene_info_file, fly_seq_fa)
