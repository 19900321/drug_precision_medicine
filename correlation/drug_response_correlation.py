import pandas as pd
import os
import pickle
from src.scripts import get_correlations_2, prepare_drug_gene
from ml_model.ml import ensembol_gene_symbol_pd

def prepare_mm_ens_symbol():
    gene_data = pd.read_csv('data/MM_log2RPKM_expression_123_20200918.txt', sep='\t', )
    gene_ens_sym_dict= dict(zip(gene_data['EnsemblID'], gene_data['Gene']))
    with(open('results/mm/ensem_symbol_dict_mm', 'wb')) as handle:
        pickle.dump(gene_ens_sym_dict, handle)


def ensembol_gene_symbol_pd_folder_mm(path_selected, col_gene):
    symol_dict = pickle.load(open('results/mm/ensem_symbol_dict_mm', 'rb'))

    for file in os.listdir(path_selected):
        data = pd.read_csv(path_selected + file)
        data = ensembol_gene_symbol_pd(data, col_gene, symol_dict=symol_dict)

        data.to_csv(path_selected + file, index=None)


def get_drug_cor_mm(path_saved):
    gene_data = pd.read_csv('data/MM_log2RPKM_expression_123_20200918.txt', sep='\t')
    patient_drug_data = pd.read_csv('data/DrugSensitivity.tsv', sep='\t')
    expr_df, drug_df = prepare_drug_gene(gene_data, patient_drug_data)
    # drugs = ['Midostaurin',
    #          'Tamatinib',
    #          'Fostamatinib'
    #          ]
    # drug_df = drug_df[drugs]
    cor_pd, cor_pd_stack, pvalue_pd, fdr_pd = get_correlations_2(expr_df, drug_df)
    cor_pd.name = 'correlation'
    cor_pd_stack.name = 'extend_all_cor_pvalue_fdr'
    pvalue_pd.name = 'pvalue'
    fdr_pd.name = 'fdr'

    if not os.path.exists(path_saved):
        os.mkdir(path_saved)

    for df in [cor_pd, cor_pd_stack, pvalue_pd, fdr_pd]:
        if df.name == 'extend_all_cor_pvalue_fdr':
            df.to_csv('{}{}.csv'.format(path_saved, df.name), index=None)
        else:
            df.to_csv('{}{}.csv'.format(path_saved, df.name))

    return cor_pd, cor_pd_stack, pvalue_pd, fdr_pd


def main():
    prepare_mm_ens_symbol()
    path_saved = 'results/mm/correlation/'
    cor_pd, cor_pd_stack, pvalue_pd, fdr_pd = get_drug_cor_mm(path_saved)
    ensembol_gene_symbol_pd_folder_mm(path_saved, 'gene')