import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from DEGs.files_preparation import prepare_groups_by_commonpass, get_feature_genes

# TODO: check by 10 genes
def prepare_pateints(dataset_gene_fpkm,drug_type_name, file_path, module_type:int , group_type:str):
    table_dict = {'bortezomib_1': {'table_name': 'bor_1line_1 treatment', 'drug_name': 'Bortezomib'},
                  'bortezomib_2': {'table_name': 'bor_Dex_1line_1treat', 'drug_name': 'Bortezomib'},
                  'bortezomib_3': {'table_name': 'bor_all line', 'drug_name': 'Bortezomib'},
                  'carf': {'table_name': 'carf_dex_all line'}, 'drug_name': 'Carfilzomib'}
    drug = table_dict[drug_type_name]['drug_name']

    final_columns_dict = get_feature_genes(file_path)
    selected_genes = final_columns_dict[drug][module_type]
    selected_genes = [c for c in selected_genes if c in list(dataset_gene_fpkm.index)]
    dataset_gene = dataset_gene_fpkm.loc[selected_genes, :]
    # generate the annotation of pateints, generate the gene data
    label_pd, dataset_gene_count_selected = prepare_groups_by_commonpass(table_dict[drug_type_name]['table_name'],
                                                                         ['sCR', 'CR', 'VGPR','PR'], ['SD', 'PD'],
                                                                         drug_type_name=drug_type_name,
                                                                         dataset_gene=dataset_gene,
                                                                         type='cluster')

    def get_subgroup(label_pd):
        pateints_pd = pd.read_excel('data/sample selection_compass.xlsx',
                                    sheet_name=table_dict[drug_type_name]['table_name'])
        pateints_dict_clinical = dict(zip(pateints_pd['public_id'], pateints_pd['bestrespsh']))
        pateints_dict_vivo =  {p.split('_')[0] +'_' +p.split('_')[1]:p for p in list(label_pd['file'])}
        subgroup_dict = {v: pateints_dict_clinical[k] for k,v in
                              pateints_dict_vivo.items() if v in list(label_pd['file'])}
        return subgroup_dict

    if group_type == 'two_group':
        group_dict = dict(zip(list(label_pd['file']),label_pd['condition']))
    elif group_type == 'sub_group':
        group_dict  = get_subgroup(label_pd)
    return dataset_gene_count_selected,group_dict


def cluster_pateints(dataset_gene_count_selected, label_dict,drug_type_name, module_type, group_type):
    # Create dictionary with features as keys and colors as values
    color_dict = {}
    palette = sns.color_palette()
    color_dict_dict = {g:palette[i] for i,g in enumerate(set(label_dict.values()))}
    for col in dataset_gene_count_selected.columns:
        if col in label_dict:
            color_dict[col] = color_dict_dict[label_dict[col]]

    # Convert the dictionary into a Series
    color_rows = pd.Series(color_dict)
    methods = ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
    'correlation', 'cosine','hamming','seuclidean','sqeuclidean']

    for m in  methods:
        g = sns.clustermap(dataset_gene_count_selected,
                           cmap="vlag",
                           z_score=0,
                           method='single',
                           metric=m,
                           linewidths=.75,
                           figsize=(50, 50),
                           col_colors=color_rows)
        # plt.show()

        plt.savefig('results/{}_cluster_{}_module{}_{}.png'.format(drug_type_name, m, module_type, group_type))