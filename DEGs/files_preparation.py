import pandas as pd
import os
import pickle
from ml_model.ml import ensembol_gene_symbol,get_final_columns_dict
from commonpass.commpass_data_process import ensembol_gene_symbol_pd_folder

def generate_annotation( patients_n, patients_p):
    lable_pd = []
    condition_record_n = list(zip(patients_n, ['negative'] * len(patients_n), ['single-read'] * len(patients_n)))
    lable_pd += condition_record_n
    condition_record_p = list(
        zip(patients_p, ['positive'] * len(patients_p), ['single-read'] * len(patients_p)))
    lable_pd += condition_record_p

    lable_pd = pd.DataFrame(lable_pd, columns=['file', 'condition', 'type'])
    return lable_pd


def prepare_file_dict_by_fimm(drugs, dataset_drug, dataset_gene):
    path_saved = 'DEGs'
    if not os.path.exists(path_saved):
        os.mkdir(path_saved)
    patient_label_dict = {}
    for d in drugs:
        patients_n = list(dataset_drug.sort_values(by=[d]).index[0:10])
        patients_p = list(dataset_drug.sort_values(by=[d]).index[-10:])
        patient_label_dict[d] = {'negative': patients_n, 'positive': patients_p}
        lable_pd = generate_annotation(patients_n, patients_p)
        lable_pd.to_csv('{}/{}_annotation.txt'.format(path_saved, d), sep='\t', index=None)
        # genrate gene expression data
        patients_selected = patients_n + patients_p
        dataset_gene_selected = dataset_gene.loc[patients_selected,:].T

        dataset_gene_selected.to_csv('{}/{}_gene.txt'.format(path_saved, d),  sep='\t')
    return patient_label_dict



def prepare_by_fimm():
    dataset_gene = pd.read_csv('data/gene_processed.csv', index_col=0)
    dataset_drug = pd.read_csv('data/drug_processed.csv', index_col=0)
    drugs = ['Bortezomib',
                  'Carfilzomib',
                  'Ixazomib',
                  'Oprozomib']

    patient_label_dict = prepare_file_dict_by_fimm(drugs, dataset_drug, dataset_gene)
    return patient_label_dict


def select_partients(sheet_name, p_term, n_term):
    pateints_pd = pd.read_excel('data/sample selection_compass.xlsx',
                  sheet_name=sheet_name)
    patients_n = list(pateints_pd.loc[pateints_pd['bestrespsh'].isin(n_term), 'public_id'])
    patients_p = list(pateints_pd.loc[pateints_pd['bestrespsh'].isin(p_term), 'public_id'])
    return patients_p, patients_n


def prepare_groups_by_commonpass(sheet_name, p_term, n_term,drug_type_name, dataset_gene, type):
    '''
    Very Good Partial Response	VGPR
    Stable Disease	SD
    Stringent Complete Response	sCR
    Partial Response	PR
    Progressive Disease	PD
    Complete Response	CR
    type = 'count', '
    '''

    patients_p, patients_n = select_partients(sheet_name, p_term, n_term)

    patients_n_selected = [c for c in dataset_gene.columns if c.startswith(tuple(patients_n))]
    patients_p_selected = [c for c in dataset_gene.columns if c.startswith(tuple(patients_p))]
    label_pd = generate_annotation(patients_n_selected, patients_p_selected)

    path_saved = 'DEGs'
    if not os.path.exists(path_saved):
        os.mkdir(path_saved)

    dataset_gene_count_selected = dataset_gene[patients_n_selected + patients_p_selected]
    if type == 'DEG':
        label_pd.to_csv('{}/{}_annotation.txt'.format(path_saved, drug_type_name), sep='\t', index=None)
        dataset_gene_count_selected.to_csv('{}/{}_gene.txt'.format(path_saved, drug_type_name), sep='\t')
    elif type == 'cluster':
        label_pd.to_csv('{}/{}_annotation_cluster.txt'.format(path_saved, drug_type_name), sep='\t', index=None)
        dataset_gene_count_selected.to_csv('{}/{}_gene_cluster.txt'.format(path_saved, drug_type_name), sep='\t')
        return label_pd,dataset_gene_count_selected


def prepare_groups_by_commonpass_processed(sheet_name, dataset_gene, type):
    '''
    Very Good Partial Response	VGPR
    Stable Disease	SD
    Stringent Complete Response	sCR
    Partial Response	PR
    Progressive Disease	PD
    Complete Response	CR
    type = 'count', '
    '''

    pateints_pd = pd.read_excel('data/sample selection_compass.xlsx',
                                sheet_name=sheet_name)
    pateints_pd = pateints_pd.drop_duplicates(subset=['group', 'gene_patient_id'], keep='first')

    patients_n_selected = list(pateints_pd.loc[pateints_pd['group'] == 'NR', 'gene_patient_id'])
    patients_p_selected = list(pateints_pd.loc[pateints_pd['group'] == 'R', 'gene_patient_id'])
    label_pd = generate_annotation(patients_n_selected, patients_p_selected)

    path_saved = 'DEGs'
    if not os.path.exists(path_saved):
        os.mkdir(path_saved)

    dataset_gene_count_selected = dataset_gene[patients_n_selected + patients_p_selected]
    if type == 'DEG':
        label_pd.to_csv('{}/{}_annotation.txt'.format(path_saved, sheet_name), sep='\t', index=None)
        dataset_gene_count_selected.to_csv('{}/{}_gene.txt'.format(path_saved, sheet_name), sep='\t')
    elif type == 'cluster':
        label_pd.to_csv('{}/{}_annotation_cluster.txt'.format(path_saved, sheet_name), sep='\t', index=None)
        dataset_gene_count_selected.to_csv('{}/{}_gene_cluster.txt'.format(path_saved, sheet_name), sep='\t')
        return label_pd, dataset_gene_count_selected


def get_feature_genes(path):
    drugs = ['Bortezomib',
             'Carfilzomib',
             'Ixazomib',
             'Oprozomib']

    data_saved = pickle.load(open(path + 'result.out', 'rb'))
    final_columns_dict = get_final_columns_dict(drugs, data_saved)

    return final_columns_dict



def map_DEG(r_method,drug_type_name):
    group_gene = pd.read_csv('DEGs/{}/{}_0.05.csv'.format(r_method,drug_type_name), index_col=0)
    sembol_genes = list(group_gene.index)
    symol_dict = ensembol_gene_symbol(sembol_genes)
    group_gene['ensmbol'] = list(group_gene.index)
    group_gene['symbol'] = group_gene['ensmbol'].apply(lambda x:symol_dict[x] if x in symol_dict else None)
    group_gene.to_csv('DEGs/{}/{}_0.05.csv'.format(r_method, drug_type_name))


def map_DEG_2():
    ensembol_gene_symbol_pd_folder('results/deg/edgeR/','Unnamed: 0')