import pandas as pd
from src.scripts import get_correlations
from src.scripts import prepare_drug_gene

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

    # choose line = 1 and therftrt = 1 in pfs filee
    data_pfs_selected = data_pfs[(data_pfs['line'] == 1) & (data_pfs['therftrt'] == 'pfs')]

    #   read os data
    data_os = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL.csv')
    selected_column_os = ['public_id',
                           'censpfs',
                           'ttcpfs',
                           'censos',
                           'ttcos']

    data_os = data_os[selected_column_os]

    # merge
    data_merge = pd.merge(data_os,
                              data_pfs,
                              left_on='public_id',
                              right_on='public_id',
                              how='outer')
    data_merge = data_merge.dropna( subset= ['gene_patient_id'])
    return data_merge


# def prepare the gene that is not lower expression
def prepare_gene_not_low_expression(gene_data):

    # only keep columns gene with at least 0.1 rpkm in 90% of your samples.
    len_genes = gene_data.shape[1]
    gene_data_left = gene_data.loc[gene_data[gene_data < 0.1].count(axis=0) / len_genes < 0.1,:]  # (921, 14680)
    print('from {} genes to {} genes '.format(gene_data.shape, gene_data_left.shape))




def main():
    def step1():
        common_patient = map_pateints_id_compass()
        common_patient.to_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP_map.csv', index=None)


    def step_2():
        data_merge = merge_pfs_os()
        data_merge.to_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL_TRTRESP_merged.csv')

    def step_3():
        comm_gene_data = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM_2.txt',
                                     sep='\t',
                                     index_col=0)


