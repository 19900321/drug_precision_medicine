import numpy as np
import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt


def plot_heatmap_cluster(cnv_data, gene):
    cnv_data.columns = ['_'.join(['cp', c])for c in cnv_data.columns]
    gene.columns = ['_'.join(['g', c]) for c in gene.columns]
    corr_pd = pd.concat([cnv_data, gene], axis=1, keys=['df1', 'df2']).corr(method='spearman').loc['df2', 'df1']

    # Draw the full plot
    f, ax = plt.subplots(figsize=(10, 5))
    sns.heatmap(corr_pd, cmap="vlag", linewidths=0.05, ax=ax)
    ax.set_xlabel('cnv of gene')

    # Set the font size and color for Y-axis labels
    ax.set_ylabel('gene expression of gene')
    ax.set_title('Correlation between gene cnv and gene expression')
    plt.tight_layout()
    plt.savefig('results/commpass/cnv/cnv.png')


def prepare_data(cnv_pd, pateint_group_dict):
    anno_gene_dict = pickle.load(open('results/mm/ensem_symbol_dict_mm', 'rb'))
    cnv_pd['symbol'] = cnv_pd['Gene'].apply(lambda x: anno_gene_dict.get(x))

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

    data = cnv_pd[cnv_pd['symbol'].isin(genes)]
    data = data.set_index('symbol')
    data = data.drop(columns=['Gene'])
    data = data.T


    all_patient = pateint_group_dict['all_patients']['patient_data']


    cnv_data = data.loc[data.index.isin(all_patient.index), :]
    all_patient = all_patient.loc[cnv_data.index]

    # prepare gene data

    gene = pateint_group_dict['all_patients']['gene_data']
    gene.columns = [anno_gene_dict.get(c) for c in list(gene.columns)]
    cols = [c for c in list(cnv_data.columns) if c in genes and c in list(gene.columns)]
    gene = gene.loc[cnv_data.index, cols]

    def fill_inf(x):

        x_mean = np.mean([i for i in x if i != np.float('-inf')])
        return [i if i != np.float('-inf') else x_mean for i in x]

    gene = gene.apply(lambda x: fill_inf(x), axis=0)
    return all_patient, gene, cnv_data


def comm_pre_group_dict_cnv(pateint_data, gene_data):
    pateint_group_dict = {
        'Bor_Dex': {'therclass': 'Bortezomib-based',
                    'thersub': 'Bor-Dex'},
        'Len_Dex':
            {'therclass': 'IMIDs-based',
             'thersub': 'Len-Dex'},
        'Bor_Cyc_Dex':
            {'therclass': 'Bortezomib-based'
                , 'thersub': 'Bor-Cyc-Dex'},
        'Bor_Len_Dex':
            {'therclass': 'combined bortezomib/IMIDs-based'
                , 'thersub': 'Bor-Len-Dex'},
        'Car_Len_Dex':
            {'therclass': 'combined IMIDs/carfilzomib-based'
                , 'thersub': 'Len-Car-Dex'},
        'Car_Cyc_Dex':
            {'therclass': 'Carfilzomib-based'
                , 'thershnm': 'Car-Cyc-Dex'},
        'Carfilzomib_based':
            {'therclass': 'Carfilzomib-based'},
        'Bor_based': {'therclass': 'Bortezomib-based'},
        'all_patients': {}
    }

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


def main():

    # prepare cnv data
    pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_log_dict', 'rb'))
    cnv_pd = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_CNA_Exome_PerGene_LargestSegment.txt',
                         sep='\t')
    all_patient, gene, cnv_data = prepare_data(cnv_pd, pateint_group_dict)
    all_patient.to_csv('results/commpass/cnv/all_patient.csv')
    gene.to_csv('results/commpass/cnv/gene.csv')
    cnv_data.to_csv('results/commpass/cnv/cnv.csv')

    # plot the correlation beween gene and cnv
    plot_heatmap_cluster(cnv_data, gene)

    # prepare groups by cnv
    pateint_data = pd.read_csv(
        'data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL_TRTRESP_merged_not_na.csv',
        index_col='gene_patient_id')
    cnv_pd = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_CNA_Exome_PerGene_LargestSegment.txt',
                         sep='\t', index_col=0).T
    pateint_group_dict_cnv_s100 = comm_pre_group_dict_cnv(pateint_data, cnv_pd)
    with(open('results/commpass/subgroup/pateint_group_dict_cnv_s100', 'wb')) as handle:
        pickle.dump(pateint_group_dict_cnv_s100, handle)










