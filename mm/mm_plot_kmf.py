import pickle
import pandas as pd
from survival.surve_kaplanmeier import compare_groups


def plot_kmf_groups_s100():

    cut_percent_s = [0.25, 0.33, 0.5]
    pateint_group_dict = pickle.load(open('results/mm/subgroup/pateint_group_log_dict_2', 'rb'))
    g_list = list(pateint_group_dict.keys())
    cens_list = ['pfs', 'os']
    genes = ['S100A1',
            'S100A2',
            'S100A3',
            'S100A4',
            'S100A5',
            'S100A6',
            'S100A7',
            'S100A7A',
            'S100A8',
            'S100A9',
            'S100A10',
            'S100A11',
            'S100A12',
            'S100A13',
            'S100A14',
            'S100A16',
            'TCHHL1',
            'S100B',
            'S100G',
            'S100P',
            'S100Z',
            'MCL1',
            'IL6R',
            'CKS1B',
            'S100PBP']

    gene_data_log = pd.read_csv('data/mm/gene_processed.csv', index_col=0)
    genes_g = list(gene_data_log.columns)
    anno_gene_dict = pickle.load(open('results/mm/ensem_symbol_dict_mm', 'rb'))
    col_genes = [k for k, v in anno_gene_dict.items() if v in genes and k in genes_g]
    path_saved = 'results/mm/mm_group_kmf_S100_2/'
    data_soure = 'mm'
    for gene in col_genes:
        for g in ['all_patients']:
            for cens_term in cens_list:
                cens = 'cens{}'.format(cens_term)
                if cens_term == 'pfs':
                    tt = 'ttc{}'.format(cens_term)
                elif cens_term == 'os':
                    tt = 'tt{}'.format(cens_term)
                for cut_percent in cut_percent_s:
                    compare_groups(gene, tt, cens, g, pateint_group_dict, cut_percent, anno_gene_dict, data_soure, path_saved)


def main():
    plot_kmf_groups_s100()