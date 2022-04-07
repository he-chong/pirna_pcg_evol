import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import seaborn as sns
import numpy as np
from collections import defaultdict


def af_compare(sense_af_list, antisense_af_list):
    sense_mean = np.mean(sense_af_list)
    antisense_mean = np.mean(antisense_af_list)
    u1, pvalue = stats.mannwhitneyu(sense_af_list, antisense_af_list, alternative='greater')
    pvalue_str_sep = str(pvalue).split('e')
    if len(pvalue_str_sep) > 1:
        pvalue_str = '{}e{}'.format(round(float(pvalue_str_sep[0]), 1), pvalue_str_sep[1])
    else:
        pvalue_str = str(round(pvalue, 4))
    u2 = len(sense_af_list)*len(antisense_af_list)-u1
    print(u1, u2, pvalue)
    return str(round(sense_mean, 4)), str(round(antisense_mean, 4)), str(u1), str(u2), pvalue_str


def build_df_for_analysis(sense_var_info, antisense_var_info, af_cutoff=0):
    sense_raw_df = pd.read_csv(sense_var_info, sep='\t')
    antisense_raw_df = pd.read_csv(antisense_var_info, sep='\t')
    sense_af_dict = defaultdict(list)
    antisense_af_dict = defaultdict(list)
    antisense_alleles = set()

    antisense_rows = []
    antisense_total = 0
    for i, row in antisense_raw_df.iterrows():
        site = row['LOC(1-based)']
        alt = row['ALT']
        try:
            af = row['Allele_frequency']
        except KeyError:
            af = row['ALTCOUNT']/(row['REFCOUNT']+row['ALTCOUNT'])
            row['Allele_frequency'] =  af
        if af <= 0:
            continue
        row['log10(Allele_frequency)'] = np.log10(af)
        try:
            category = row['Category']
        except KeyError:
            category = None
        row['Strand'] = 'Antisense'
        antisense_alleles.add((site, alt))
        antisense_total += 1
        if af > af_cutoff:
            antisense_rows.append(row)
    print(len(antisense_rows), antisense_total, len(antisense_rows)/antisense_total)
    
    sense_rows = []
    sense_total = 0
    for i, row in sense_raw_df.iterrows():
        site = row['LOC(1-based)']
        alt = row['ALT']
        try:
            af = row['Allele_frequency']
        except KeyError:
            af = row['ALTCOUNT']/(row['REFCOUNT']+row['ALTCOUNT'])
            row['Allele_frequency'] = af
        if af <= 0:
            continue
        row['log10(Allele_frequency)'] = np.log10(af)
        try:
            category = row['Category']
        except KeyError:
            category = None
        row['Strand'] = 'Sense'
        if (site, alt) not in antisense_alleles:
            sense_total += 1
            if af > af_cutoff:
                sense_rows.append(row)
    print(len(sense_rows), sense_total, len(sense_rows)/sense_total)

    sense_df = pd.concat(sense_rows, axis=1).T.reset_index(drop=True)
    antisense_df = pd.concat(antisense_rows, axis=1).T.reset_index(drop=True)
    total_df = pd.concat([sense_df, antisense_df]).reset_index(drop=True)
    return sense_df, antisense_df, total_df


def compare_sense_and_antisense(sense_df, antisense_df, total_df, hist_ax, box_ax, out_table, category_order=None, category_labels=None):
    sns.histplot(data=total_df, x='log10(Allele_frequency)', hue='Strand', ax=hist_ax, stat='density', common_bins=True, common_norm=False, multiple='layer', bins=20, log_scale=False, palette='Set2')
    # hist_ax.legend(loc='upper right')
    if category_order:
        sns.boxplot(data=total_df, x='Category', y='log10(Allele_frequency)', hue='Strand', ax=box_ax, order=category_order, palette='Set2')
    
        box_ax.legend(loc='lower right')
        box_ax.set_xticklabels(category_labels)
        box_ax.set_xlabel('')
    with open(out_table, 'w') as out_handle:
        out_handle.write('\t'.join(['Category', 'Sense mean', 'Antisense mean', 'U1', 'U2', 'Mann-Whitney U-test\'s p-value']))
        out_handle.write('\n')
        result = af_compare(sense_df['log10(Allele_frequency)'], antisense_df['log10(Allele_frequency)'])
        out_handle.write('All')
        out_handle.write('\t')
        out_handle.write('\t'.join(result))
        out_handle.write('\n')
        if category_order:
            for category in category_order:
                category_sense_df = sense_df.loc[sense_df['Category']==category]
                category_antisense_df = antisense_df.loc[antisense_df['Category']==category]
                category_result = af_compare(category_sense_df['log10(Allele_frequency)'], category_antisense_df['log10(Allele_frequency)'])
                out_handle.write(category)
                out_handle.write('\t')
                out_handle.write('\t'.join(category_result))
                out_handle.write('\n')
    

def compare_MAs_and_non_MAs(sense_ma_df, antisense_ma_df, sense_nonma_df, antisense_nonma_df, out_table, category_list=None):
    with open(out_table, 'w') as out_handle:
        out_handle.write('Sense MAs vs. sense non-MAs\n')
        out_handle.write('\t'.join(['Category', 'MA mean', 'Non-MA mean', 'U1', 'U2', 'Mann-Whitney U-test\'s p-value']))
        out_handle.write('\n')
        sense_result = af_compare(sense_ma_df['log10(Allele_frequency)'], sense_nonma_df['log10(Allele_frequency)'])
        out_handle.write('All')
        out_handle.write('\t')
        out_handle.write('\t'.join(sense_result))
        out_handle.write('\n')
        if category_list:
            for category in category_list:
                category_sense_ma_df = sense_ma_df.loc[sense_ma_df['Category']==category]
                category_sense_nonma_df = sense_nonma_df.loc[sense_nonma_df['Category']==category]
                sense_category_result = af_compare(category_sense_ma_df['log10(Allele_frequency)'], category_sense_nonma_df['log10(Allele_frequency)'])
                out_handle.write(category)
                out_handle.write('\t')
                out_handle.write('\t'.join(sense_category_result))
                out_handle.write('\n')
        out_handle.write('\n\n')

        out_handle.write('Antisense MAs vs. antisense non-MAs\n')
        out_handle.write('\t'.join(['Category', 'MA mean', 'Non-MA mean', 'U1', 'U2', 'Mann-Whitney U-test\'s p-value']))
        out_handle.write('\n')
        antisense_result = af_compare(antisense_ma_df['log10(Allele_frequency)'], antisense_nonma_df['log10(Allele_frequency)'])
        out_handle.write('All')
        out_handle.write('\t')
        out_handle.write('\t'.join(antisense_result))
        out_handle.write('\n')
        if category_list:
            for category in category_list:
                category_antisense_ma_df = antisense_ma_df.loc[antisense_ma_df['Category']==category]
                category_antisense_nonma_df = antisense_nonma_df.loc[antisense_nonma_df['Category']==category]
                antisense_category_result = af_compare(category_antisense_ma_df['log10(Allele_frequency)'], category_antisense_nonma_df['log10(Allele_frequency)'])
                out_handle.write(category)
                out_handle.write('\t')
                out_handle.write('\t'.join(antisense_category_result))
                out_handle.write('\n')


def stat_af(var_info_dir, out_dir, var_type):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    var_info_base = os.path.basename(var_info_dir)
    out_pdf_file = os.path.join(out_dir, '{}.{}.pdf'.format(var_info_base, var_type))
    ma_out_table = os.path.join(out_dir, '{}_magnetizing_stat.{}.txt'.format(var_info_base, var_type))
    nonma_out_table = os.path.join(out_dir, '{}_non_magnetizing_stat.{}.txt'.format(var_info_base, var_type))
    ma_nonma_out_table = os.path.join(out_dir, '{}_MA_non_MA_stat.{}.txt'.format(var_info_base, var_type))
    with PdfPages(out_pdf_file) as out_pdf:
        sense_magnetizing_var_info = os.path.join(var_info_dir, var_type+'.sense.magnetizing_var_info.txt')
        antisense_magnetizing_var_info = os.path.join(var_info_dir, var_type+'.antisense.magnetizing_var_info.txt')
        sense_non_magnetizing_var_info = os.path.join(var_info_dir, var_type+'.sense.non_magnetizing_var_info.txt')
        antisense_non_magnetizing_var_info = os.path.join(var_info_dir, var_type+'.antisense.non_magnetizing_var_info.txt')

        sense_ma_df, antisense_ma_df, total_ma_df = build_df_for_analysis(sense_magnetizing_var_info, antisense_magnetizing_var_info)
        sense_nonma_df, antisense_nonma_df, total_nonma_df = build_df_for_analysis(sense_non_magnetizing_var_info, antisense_non_magnetizing_var_info)

        if var_type == 'pss':
            fig = plt.figure(figsize=(9, 5))
            fig.subplots_adjust(hspace=0.3)
            axes = fig.subplots(nrows=2, ncols=2)   
            category_order = ['5\'UTR', 'Synonymous', 'Non-synonymous','3\'UTR']
            category_labels = ['5\'UTR', 'Syn.', 'Non-syn.','3\'UTR']
            compare_sense_and_antisense(sense_ma_df, antisense_ma_df, total_ma_df, axes[0][0], axes[0][1], ma_out_table, category_order, category_labels)
            compare_sense_and_antisense(sense_nonma_df, antisense_nonma_df, total_nonma_df,  axes[1][0], axes[1][1], nonma_out_table, category_order, category_labels)
        else:
            fig = plt.figure(figsize=(4.5, 5))
            fig.subplots_adjust(hspace=0.3)
            axes = fig.subplots(nrows=2, ncols=1)  
            compare_sense_and_antisense(sense_ma_df, antisense_ma_df, total_ma_df, axes[0], None, ma_out_table)
            compare_sense_and_antisense(sense_nonma_df, antisense_nonma_df, total_nonma_df, axes[1], None, nonma_out_table)
        
        out_pdf.savefig()
    plt.close()
    if var_type == 'pss':
        compare_MAs_and_non_MAs(sense_ma_df, antisense_ma_df, sense_nonma_df, antisense_nonma_df, ma_nonma_out_table, category_list=['5\'UTR', 'Synonymous', 'Non-synonymous','3\'UTR'])
    else:
        compare_MAs_and_non_MAs(sense_ma_df, antisense_ma_df, sense_nonma_df, antisense_nonma_df, ma_nonma_out_table)



def stat_human_1000_genomes():
    var_info_dir = r'..\raw_data\var_in_pirna_pss_pairs\1000_Genomes_all'
    out_dir = r'..\raw_data\var_stat'
    stat_af(var_info_dir, out_dir, 'pss')
    stat_af(var_info_dir, out_dir, 'pirna')


def stat_human_1000_genomes_complicated():
    var_info_dir = r'..\raw_data\var_in_pirna_pss_pairs\1000_Genomes_complicated'
    out_dir = r'..\raw_data\var_stat'
    stat_af(var_info_dir, out_dir, 'pss')


def stat_human_uk10k():
    var_info_dir = r'..\raw_data\var_in_pirna_pss_pairs\UK10K_all'
    out_dir = r'..\raw_data\var_stat'
    stat_af(var_info_dir, out_dir, 'pss')
    stat_af(var_info_dir, out_dir, 'pirna')


def stat_human_uk10k_complicated():
    var_info_dir = r'..\raw_data\var_in_pirna_pss_pairs\UK10K_complicated_mm5'
    out_dir = r'..\raw_data\var_stat'
    stat_af(var_info_dir, out_dir, 'pss')


def stat_d_melanogaster():
    var_info_dir = r'..\raw_data\var_in_pirna_pss_pairs\DGRP2_simple'
    out_dir = r'..\raw_data\var_stat'
    stat_af(var_info_dir, out_dir, 'pss')


if __name__ == '__main__':
    stat_human_1000_genomes()
    stat_human_1000_genomes_complicated()
    stat_human_uk10k()
    stat_human_uk10k_complicated()
    stat_d_melanogaster()