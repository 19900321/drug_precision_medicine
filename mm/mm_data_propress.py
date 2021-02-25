import pandas as pd
import numpy as np
import pickle
import os


def prepare_groups_mm():
    pateint_data = pd.read_excel('data/mm/FIMM clinical dataset_20210224_yin.xls',
                                   sheet_name='Sheet1', index_col=0)
    pateint_data = pateint_data[pateint_data['RNA sequencing done'] == 1]
    idx = pateint_data.groupby('Patient ID')['SampleNo'].transform(min) == pateint_data['SampleNo']
    pateint_data_2 = pateint_data[idx]
    pateint_data_2['sample_id'] = pateint_data_2['sample_id'].apply(lambda x: x.replace('+', 'P'))
    pateint_data_2['map_short_id'] = pateint_data_2['shortid'].apply(lambda x:'_'.join(x.split('_')[1:]))

    # another survival data
    all_mm = pd.read_csv('data/mm/survival_allMM_censos.txt', sep='\t')
    pateint_data_2 = pd.merge(pateint_data_2, all_mm, right_on='shortid', left_on='map_short_id', how='left')

    col_select_1 = {
    'sample_id': 'sample_id',
    'OS (mo) from diagnosis': 'ttos',
    'Dead: patient is dead': 'censos',
    'Next line treatment PFS (months)': 'ttcpfs',
    'PFSCheck: Has progression occured during next line treatment': 'censpfs'}

    col_select_2 = {
                    'sample_id': 'sample_id',
                    'OS_months': 'ttos',
                    'census_OS': 'censos',
                    'Next line treatment PFS (months)': 'ttcpfs',
                    'PFSCheck: Has progression occured during next line treatment': 'censpfs'}

    # read gene data
    gene_data_log = pd.read_csv('data/mm/gene_processed.csv', index_col=0)

    n = 1
    for col_select in [col_select_1, col_select_2]:
        pateint_data_one = pateint_data_2[list(col_select.keys())]
        pateint_data_one = pateint_data_one.rename(columns=col_select)
        pateint_data_one['censos'] = pateint_data_one['censos'].astype(int)
        pateint_data_one.index = pateint_data_one['sample_id']

        # prepare
        sub_group_p_perform = pateint_data_one.loc[pateint_data_one.index.isin(gene_data_log.index), :]

        # keep gene data of corresponding patients
        sub_group_g_perform = gene_data_log.loc[sub_group_p_perform.index]

        # add to group dictionary for further search
        pateint_group_dict = {}
        pateint_group_dict['all_patients'] = {'patient_data': sub_group_p_perform,
                                      'gene_data': sub_group_g_perform}

        with(open('results/mm/subgroup/pateint_group_log_dict_{}'.format(n), 'wb')) as handle:
            pickle.dump(pateint_group_dict, handle)
        n += 1

