import os, glob
from Bio import AlignIO
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt


def most_similar_block(aln, block_size):
    max_identity = 0
    for i in range(aln.get_alignment_length()-block_size):
        sub_aln = aln[:,i:i+block_size]
        ident_num = 0
        total = 0
        for nt_1, nt_2 in zip(sub_aln[0].seq, sub_aln[1].seq):
            total += 1
            # print(nt_1.upper(), nt_2.upper())
            if nt_1.upper() == nt_2.upper():
                ident_num += 1
        identity = ident_num/total
        if identity > max_identity:
            max_identity = identity
    return i, round(max_identity, 4)


def stat_expanded_similarity(up_expanded_aln_dir, down_expanded_aln_dir, overall_expanded_aln_dir, up_block_size, down_block_size, overall_block_size, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    for up_repeat in sorted(os.listdir(up_expanded_aln_dir), key=lambda i:int(i.split('_')[1])):
        out_file = os.path.join(out_dir, up_repeat+'.txt')
        repeat_dir = os.path.join(up_expanded_aln_dir, up_repeat)
        with open(out_file, 'w') as out_handle:
            for expanded_aln_base in os.listdir(repeat_dir):
                up_expanded_aln_file = os.path.join(up_expanded_aln_dir, up_repeat, expanded_aln_base)
                down_expanded_aln_file = os.path.join(down_expanded_aln_dir, up_repeat, expanded_aln_base)
                overall_expanded_aln_file = os.path.join(down_expanded_aln_dir, up_repeat, expanded_aln_base)
                up_aln = AlignIO.read(up_expanded_aln_file, 'fasta')
                down_aln = AlignIO.read(down_expanded_aln_file, 'fasta')
                overall_aln = AlignIO.read(overall_expanded_aln_file, 'fasta')
                up_most_similar_block_start, up_ident = most_similar_block(up_aln, up_block_size)
                down_most_similar_block_start, down_ident = most_similar_block(down_aln, down_block_size)
                overall_most_similar_block_start, overall_ident = most_similar_block(overall_aln, overall_block_size)
                out_handle.write('\t'.join(map(str, [expanded_aln_base.strip('.fa'), up_most_similar_block_start, up_ident, down_most_similar_block_start, down_ident, overall_most_similar_block_start, overall_ident])))
                out_handle.write('\n')


def draw_expanded_similarity(similarity_dir, out_pdf_file):
    up_most_similar_block_start_index, up_ident_index, down_most_similar_block_start_index, down_ident_index, overall_most_similar_block_start_index, overall_ident_index = range(1,7)
    category_order = ['Up', 'Down', 'Overall']
    out_dir = os.path.dirname(out_pdf_file)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with PdfPages(out_pdf_file) as out_pdf:
        row_num = 5
        col_num = 6
        fig = plt.figure(figsize=(12, 9))
        axes = fig.subplots(row_num, col_num, sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0.4)
        i = 0
        j = 0
        n = 1
        for similarity_base in os.listdir(similarity_dir):
            similarity_file = os.path.join(similarity_dir, similarity_base)
            rows_for_drawing = []
            total_up_down_prop = 0
            up_down_90_prop = 0 
            with open(similarity_file) as similarity_handle:
                for line in similarity_handle:
                    line_split = line.strip().split('\t')
                    up_most_similar_block_start = int(line_split[up_most_similar_block_start_index])
                    up_ident = round(float(line_split[up_ident_index]), 4)
                    down_most_similar_block_start = int(line_split[down_most_similar_block_start_index])
                    down_ident = round(float(line_split[down_ident_index]), 4)
                    overall_most_similar_block_start = int(line_split[overall_most_similar_block_start_index])
                    overall_ident = round(float(line_split[overall_ident_index]), 4)
                    rows_for_drawing.append(pd.Series({'Block_start':up_most_similar_block_start, 'Identity':up_ident, 'Category':'Up'}))
                    rows_for_drawing.append(pd.Series({'Block_start':down_most_similar_block_start, 'Identity':down_ident, 'Category':'Down'}))
                    rows_for_drawing.append(pd.Series({'Block_start':overall_most_similar_block_start, 'Identity':overall_ident, 'Category':'Overall'}))
                    total_up_down_prop += 1
                    if float(line_split[up_ident_index]) > 0.9 and float(line_split[down_ident_index]) > 0.9:
                        up_down_90_prop += 1
            df_for_drawing = pd.concat(rows_for_drawing, axis=1).T.reset_index(drop=True)


            ax = axes[i][j]
            ax.set_title('Repeat_{}'.format(n))
            # mean = round(np.mean(similarity_list)*100, 2)
            up_mean = round(np.median(df_for_drawing.loc[df_for_drawing['Category']=='Up']['Identity'])*100, 2)
            down_mean = round(np.median(df_for_drawing.loc[df_for_drawing['Category']=='Down']['Identity'])*100, 2)
            overall_mean = round(np.median(df_for_drawing.loc[df_for_drawing['Category']=='Overall']['Identity'])*100, 2)
            sns.boxplot(data=df_for_drawing, x='Category', y='Identity', ax=ax, order=category_order, palette='Set2', flierprops={'marker':'x', 'markersize': 4, 'markerfacecolor':'gray', 'markeredgecolor':'lightgray'})
            ax.text(0.02, 0.314, 'Up: {}%\nDown: {}%\nOverall: {}%\nFlanking_prop: {}'.format(up_mean, down_mean, overall_mean, up_down_90_prop/total_up_down_prop), fontsize=6)
            # ax.set_xticklabels(['Repeat_{}'.format(n)])
            ax.set_xlabel('')
            if j > 0:
                ax.set_ylabel('')
            j += 1
            if j >= col_num:
                j = j % col_num
                i += 1
            n += 1
        out_pdf.savefig()
    

def stat_human():
    sense_up_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_150/human/sense/up_expanded_aln'
    sense_down_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_150/human/sense/down_expanded_aln'
    sense_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_150/human/sense/expanded_aln'
    antisense_up_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_150/human/antisense/up_expanded_aln'
    antisense_down_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_150/human/antisense/down_expanded_aln'
    antisense_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs_150/human/antisense/expanded_aln'
    sense_out_dir = '../generated_data/expanded_similarities_stat/human/sense'
    antisense_out_dir = '../generated_data/expanded_similarities_stat/human/antisense'
    sense_out_pdf = '../generated_data/expanded_similarities_pdf/human/sense.pdf'
    antisense_out_pdf = '../generated_data/expanded_similarities_pdf/human/antisense.pdf'
    up_block_size = 40
    down_block_size = 40
    overall_block_size = 80
    stat_expanded_similarity(sense_up_expanded_aln_dir, sense_down_expanded_aln_dir, sense_expanded_aln_dir, up_block_size, down_block_size, overall_block_size, sense_out_dir)
    stat_expanded_similarity(antisense_up_expanded_aln_dir, antisense_down_expanded_aln_dir, antisense_expanded_aln_dir, up_block_size, down_block_size, overall_block_size, antisense_out_dir)
    draw_expanded_similarity(sense_out_dir, sense_out_pdf)
    draw_expanded_similarity(antisense_out_dir, antisense_out_pdf)


def stat_d_melanogaster():
    sense_up_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs/d_melanogaster/sense/up_expanded_aln'
    sense_down_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs/d_melanogaster/sense/down_expanded_aln'
    sense_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs/d_melanogaster/sense/expanded_aln'
    antisense_up_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs/d_melanogaster/antisense/up_expanded_aln'
    antisense_down_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs/d_melanogaster/antisense/down_expanded_aln'
    antisense_expanded_aln_dir = '../generated_data/expanded_pirna_pss_pairs/d_melanogaster/antisense/expanded_aln'
    sense_out_dir = '../generated_data/expanded_similarities_stat/d_melanogaster/sense'
    antisense_out_dir = '../generated_data/expanded_similarities_stat/d_melanogaster/antisense'
    sense_out_pdf = '../generated_data/expanded_similarities_pdf/d_melanogaster/sense.pdf'
    antisense_out_pdf = '../generated_data/expanded_similarities_pdf/d_melanogaster/antisense.pdf'
    up_block_size = 35
    down_block_size = 35
    all_block_size = 100
    stat_expanded_similarity(sense_up_expanded_aln_dir, sense_down_expanded_aln_dir, sense_expanded_aln_dir, up_block_size, down_block_size, all_block_size, sense_out_dir)
    stat_expanded_similarity(antisense_up_expanded_aln_dir, antisense_down_expanded_aln_dir, antisense_expanded_aln_dir, up_block_size, down_block_size, all_block_size, antisense_out_dir)
    draw_expanded_similarity(sense_out_dir, sense_out_pdf)
    draw_expanded_similarity(antisense_out_dir, antisense_out_pdf)


if __name__ == '__main__':
    stat_human()
    # stat_d_melanogaster()

            