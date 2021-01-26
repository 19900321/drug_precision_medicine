import os
import pandas as pd
import pickle
from src.scripts import prepare_drug_gene
from scipy.stats import ttest_rel

def prepare_file_dict_by_fimm(drugs, dataset_drug, dataset_gene):
    path_saved = 'results/mm/t_test/'
    if not os.path.exists(path_saved):
        os.mkdir(path_saved)
    anno_gene_dict = pickle.load(open('results/commpass/ensembol_symbol_dict_commpass', 'rb'))
    col_genes = list(set(['ENSG00000116497',
                 'ENSG00000143546',
                 'ENSG00000160678',
                 'ENSG00000163221',
                 'ENSG00000197747',
                 'ENSG00000163191',
                 'ENSG00000197956',
                 'ENSG00000163220',
                 'ENSG00000196154',
                 'ENSG00000189171']) & set(dataset_gene.columns))

    patient_label_dict = {}
    p_value_pd= []
    for d in drugs:
        patients_n = list(dataset_drug.sort_values(by=[d]).index[0:10])
        patients_p = list(dataset_drug.sort_values(by=[d]).index[-10:])

        # genrate gene expression data
        dataset_gene_pos = dataset_gene.loc[patients_p, col_genes]
        dataset_gene_n = dataset_gene.loc[patients_n, col_genes]

        patient_label_dict[d] = {'negative': patients_n,
                                 'positive': patients_p,
                                 'dataset_gene_pos': dataset_gene_pos,
                                 'dataset_gene_n': dataset_gene_n}
        for g in col_genes:
            p_value = ttest_rel(dataset_gene_pos[g], dataset_gene_n[g]).pvalue
            p_value_pd.append([d, g, p_value])

    p_value_pd = pd.DataFrame(p_value_pd, columns=['drug', 'gene','p_value'])
    p_value_pd['symbol'] = p_value_pd['gene'].apply(lambda x:anno_gene_dict.get(x))
    p_value_pd.to_csv('{}mm_pateint_test.csv'.format(path_saved), index=None)
    with(open('{}mm_pateint_group_dict'.format(path_saved), 'wb')) as handle:
        pickle.dump(patient_label_dict, handle)

    return patient_label_dict


def main():
    drugs = ['Bortezomib',
                     'Carfilzomib',
                     'Ixazomib',
                     'Oprozomib']
    gene_data = pd.read_csv('data/MM_log2RPKM_expression_123_20200918.txt', sep='\t')
    patient_drug_data = pd.read_csv('data/DrugSensitivity.tsv', sep='\t')
    expr_df, drug_df = prepare_drug_gene(gene_data, patient_drug_data)
    patient_label_dict = prepare_file_dict_by_fimm(drugs, drug_df, expr_df)
