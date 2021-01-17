import pandas as pd
import GEOparse
import pickle
import numpy as np
import os
from ml_model.ml import ensembol_gene_symbol_pd


def prepare_pateint_info(ges_obj):
    p_info = {}
    for p_id, values in ges_obj.gsms.items():
        p_dict = {i.split(' = ')[0]:i.split(' = ')[1] for i in values.metadata['characteristics_ch1']}
        p_info[p_id] = p_dict

    info_pd = pd.DataFrame.from_dict(p_info, orient='index')
    replace_dict = {'PGx_Progression(0=No,1=Yes)': 'censpfs',
    'PGx_Days_To_Progression': 'ttcpfs',
    'Did_Patient_Die(0=No,1=Yes)': 'censos',
    'Days_Survived_From_Randomization': 'ttos'}
    info_pd = info_pd.rename(columns=replace_dict)
    info_pd.insert(0, column='sample_id', value=list(info_pd.index))
    return info_pd

def prepare_gene_data(ges_obj, platform:str):

    GPL_sample = sorted([p_id for p_id, values in ges_obj.gsms.items()
             if len(values.metadata['platform_id']) ==1
             and values.metadata['platform_id'][0] == platform])

    # elect only GPL 96
    pivoted_control_samples = ges_obj.pivot_samples('VALUE')[GPL_sample].drop_duplicates()

    # make it to log2
    pivoted_control_samples = np.log2(pivoted_control_samples)

    # filter out probes that are not expressed.  The gene is expressed (in definition here) when its average log2
    # intensity in control samples is above 0.25 quantile. I.e. we filter out worst 25% genes
    pivoted_control_samples_average = pivoted_control_samples.median(axis=1)
    print("Number of probes before filtering: ", len(pivoted_control_samples_average))

    expression_threshold = pivoted_control_samples_average.quantile(0.25)
    expressed_probes = pivoted_control_samples_average[
        pivoted_control_samples_average >= expression_threshold].index.tolist()
    print("Number of probes above threshold: ", len(expressed_probes))

    pivoted_control_samples_processed = pivoted_control_samples.loc[expressed_probes, :]

    # MAKE THE PATEINT_ID AS INDEX, THE GENEN AS COLUMNS
    pivoted_control_samples_processed_trans = pivoted_control_samples_processed.T
    return pivoted_control_samples_processed_trans


def prepare_subgroup_geo(pateint_data, gene_data):
    pateint_group_dict = {
            'PS341': {'treatment': 'PS341'},
            'Dex': {'treatment': 'Dex'}}

    for g, v_1 in pateint_group_dict.items():
        sub_group_p_perform = pateint_data.loc[:, :]
        # select by condition of pateints
        for term, v_2 in v_1.items():
            sub_group_p_perform = sub_group_p_perform[sub_group_p_perform[term] == v_2]

        # keep only those with gene data
        sub_group_p_perform = sub_group_p_perform.loc[sub_group_p_perform.index.isin(gene_data.index), :]

        # keep gene data of corresponding patients
        sub_group_g_perform = gene_data.loc[sub_group_p_perform.index]

        # add to group dictionary for further search
        pateint_group_dict[g].update({'patient_data': sub_group_p_perform,
                                      'gene_data': sub_group_g_perform})

    return pateint_group_dict


def prepare_ensembol_symbol_dict_commpass(ges_obj, platform:str):
    gene_info_pd = ges_obj.gpls[platform].table
    symol_dict = dict(zip(gene_info_pd['ID'], gene_info_pd['Gene Symbol']))
    return symol_dict


def ensembol_gene_symbol_pd_folder_geo(path_selected, col_gene):
    symol_dict = pickle.load(open('results/geo/id_symbol_dict_geo', 'rb'))

    for file in os.listdir(path_selected):
        data = pd.read_csv(path_selected + file)
        data = ensembol_gene_symbol_pd(data, col_gene, symol_dict=symol_dict)

        data.to_csv(path_selected + file, index=None)



def main():
    gse_9782 = GEOparse.get_GEO(filepath='data/geo/GSE9782_family.soft.gz')
    info_pd = prepare_pateint_info(gse_9782)
    info_pd.to_csv('data/geo/ges_8782_pateint_info.csv', index=None)
    pivoted_control_samples_processed_trans = prepare_gene_data(gse_9782, 'GPL96')
    pivoted_control_samples_processed_trans.to_csv('data/geo/ges_8782_gene_data.csv')

    # group the pateints
    pateint_data = pd.read_csv('data/geo/ges_8782_pateint_info.csv', index_col=0)
    gene_data = pd.read_csv('data/geo/ges_8782_gene_data.csv', index_col=0)
    pateint_group_dict = prepare_subgroup_geo(pateint_data, gene_data)

    with(open('results/geo/subgroup/pateint_group_dict', 'wb')) as handle:
        pickle.dump(pateint_group_dict, handle)

    symol_dict = prepare_ensembol_symbol_dict_commpass(gse_9782, 'GPL96')
    with(open('results/geo/id_symbol_dict_geo', 'wb')) as handle:
        pickle.dump(symol_dict, handle)

