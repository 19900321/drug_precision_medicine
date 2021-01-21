import pandas as pd
import numpy as np
import pickle
import os
from src.scripts import get_correlations
from src.scripts import prepare_drug_gene
from ml_model.ml import ensembol_gene_symbol, ensembol_gene_symbol_pd


# step 1: get the patients with RNA seq based on gene MMRF_CoMMpass_IA15a_E74GTF_HtSeq_Gene_Counts.txt
# and clinical parameter MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP.csv
def map_pateints_id_compass():
    common_patient = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP.csv', sep=';')
    dataset_gene_count = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_HtSeq_Gene_Counts.txt',
                                     sep='\t',
                                     index_col=0)

    gene_pateints = list(dataset_gene_count.columns)

    gene_patient_id = []
    for i, rows in common_patient.iterrows():
        test_gene_patient_id = rows.public_id + '_' + str(rows.line) + '_BM'
        if test_gene_patient_id in gene_pateints:
            gene_patient_id.append(test_gene_patient_id)
        else:
            test_gene_patient_id_alter = rows.public_id + '_' + str(rows.line) + '_PB'
            if test_gene_patient_id_alter in gene_pateints:
                gene_patient_id.append(test_gene_patient_id_alter)
            else:
                gene_patient_id.append(None)

    common_patient['gene_patient_id'] = gene_patient_id

    return common_patient

#step 2. extract columns from the above MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP_map.csv
# and MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL.csv
def merge_pfs_os():
    data_pfs = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP_map.csv')
    selected_col_pfs = ['public_id',
                        'line',
                        'gene_patient_id',
                        'thershnm',
                        'therclass',
                        'thersub',
                        'bestrespsh',
                        'therftrt']

    data_pfs = data_pfs[selected_col_pfs]

    # choose line = 1 and therftrt = 1 in pfs file
    data_pfs_selected = data_pfs[(data_pfs['line'] == 1) & (data_pfs['therftrt'] == 1)]

    #   read os data
    data_os = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL.csv')
    selected_column_os = ['public_id',
                           'censpfs',
                           'ttcpfs',
                           'censos',
                           'ttos']

    data_os = data_os[selected_column_os]

    # merge
    data_merge = pd.merge(data_os,
                              data_pfs_selected,
                              left_on='public_id',
                              right_on='public_id',
                              how='outer')
    data_merge = data_merge.dropna( subset= ['gene_patient_id'])
    return data_merge


#  prepare the gene that is not lower expression，
def prepare_gene_not_low_expression(data_gene):

    # only keep columns gene with at least 0.1 rpkm in 90% of your samples.
    len_genes = data_gene.shape[0]
    data_left = data_gene.loc[:, data_gene[data_gene < 0.1].count() / len_genes < 0.1]
    print('from {} genes to {} genes '.format(data_gene.shape, data_left.shape))
    # as log 2 FKPM is important, so we prepare the log2 FKPM
    data_left_log = np.log2(data_left)
    return data_left, data_left_log


def prepare_subgroups(pateint_data, gene_data):
    pateint_group_dict = {
        'Bor_Dex':{'therclass': 'Bortezomib-based',
           'thersub': 'Bor-Dex'},
        'Len_Dex':
            {'therclass': 'IMIDs-based',
            'thersub': 'Len-Dex'},
        'Bor_Cyc_Dex':
            {'therclass': 'Bortezomib-based'
            ,'thersub': 'Bor-Cyc-Dex'},
        'Bor_Len_Dex':
            {'therclass': 'combined bortezomib/IMIDs-based'
            ,'thersub': 'Bor-Len-Dex'},
        'Car_Len_Dex':
            {'therclass': 'combined IMIDs/carfilzomib-based'
            , 'thersub': 'Len-Car-Dex'},
        'Car_Cyc_Dex':
            {'therclass': 'Carfilzomib-based'
            ,'thershnm': 'Car-Cyc-Dex'},
        'Carfilzomib_based':
            {'therclass': 'Carfilzomib-based'},
        'Bor_based':{'therclass': 'Bortezomib-based'},
        'all_patients':{}
    }

    for g, v_1 in pateint_group_dict.items():
        sub_group_p_perform = pateint_data.loc[:,:]
        # select by condition of pateints
        for term, v_2 in v_1.items():
            sub_group_p_perform = sub_group_p_perform[sub_group_p_perform[term] == v_2]

        # keep only those with gene data
        sub_group_p_perform = sub_group_p_perform.loc[sub_group_p_perform.index.isin(gene_data.index),:]

        # keep gene data of corresponding patients
        sub_group_g_perform = gene_data.loc[sub_group_p_perform.index]

        # add to group dictionary for further search
        pateint_group_dict[g].update({'patient_data':sub_group_p_perform,
                                                           'gene_data':sub_group_g_perform})

    return pateint_group_dict

def prepare_ensembol_symbol_dict_commpass(sembol_genes):
    symol_dict = ensembol_gene_symbol(sembol_genes)
    with(open('results/commpass/ensembol_symbol_dict_commpass', 'wb')) as handle:
        pickle.dump(symol_dict, handle)


def ensembol_gene_symbol_pd_folder(path_selected, col_gene):
    symol_dict = pickle.load(open('results/commpass/ensembol_symbol_dict_commpass', 'rb'))

    for file in os.listdir(path_selected):
        data = pd.read_csv(path_selected + file)
        data = ensembol_gene_symbol_pd(data, col_gene, symol_dict=symol_dict)

        data.to_csv(path_selected + file, index=None)



def main():

    # map the patients ids which can be mapped in patient public id and
    def step1():
        common_patient = map_pateints_id_compass()
        common_patient.to_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP_map.csv', index=None)

    # merge pfs and os together
    def step_2():
        data_merge = merge_pfs_os()
        data_merge.to_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL_TRTRESP_merged_not_na.csv',
                          index = None)
    # DELETE LOWER EXPRESSION GENE
    def step_3():
        comm_gene_data = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM.txt',
                                     sep='\t',
                                     index_col=0)
        comm_gene_data_left, comm_gene_data_left_log = prepare_gene_not_low_expression(comm_gene_data)
        comm_gene_data_left.to_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM_processed.txt',
                                   sep='\t')
        comm_gene_data_left_log.to_csv(
            'data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_log_2_FPKM_processed.txt',
            sep='\t')

    # PREPARE SUBGROUP PATEINTS
    def step_4():
        # keep the index is pateint id
        pateint_data = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL_TRTRESP_merged_not_na.csv',
                                   index_col='gene_patient_id')

        gene_data = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM_processed.txt',
                                   sep='\t', index_col=0)
        gene_data_log = pd.read_csv(
            'data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_log_2_FPKM_processed.txt',
            sep='\t', index_col=0)

        pateint_group_dict = prepare_subgroups(pateint_data, gene_data)
        pateint_group_log_dict = prepare_subgroups(pateint_data, gene_data_log)
        with(open('results/commpass/subgroup/pateint_group_dict', 'wb')) as handle:
            pickle.dump(pateint_group_dict, handle)

        with(open('results/commpass/subgroup/pateint_group_log_dict', 'wb')) as handle:
            pickle.dump(pateint_group_log_dict, handle)


