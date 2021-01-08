import pickle
import itertools
from ml_model.ml import uniprot_gene_symbol,ensembol_gene_symbol
import pandas as pd


def get_single_drug_fdr(drugs):
    classic = pd.read_csv('results/mm106_data_cor_pvalue_0.05.csv', sep=';')
    classic = classic[round(classic['adjust_pvalue'], 2) <0.1]
    classic_p = classic[round(classic['coefficient'], 2) >0]
    classic_n = classic[round(classic['coefficient'], 2) <0]
    drug_classic_dict_p = dict(classic_p.groupby('drug')['Gene'].apply(list))
    drug_classic_dict_n = dict(classic_n.groupby('drug')['Gene'].apply(list))
    classic_gene = {}
    for d in drugs:
        classic_gene[d] = {'final_pcm_symbol':drug_classic_dict_p[d],
                           'final_ncm_symbol':drug_classic_dict_n[d]}
    return classic_gene

# get the cinical DEGs

def get_group_fdr_gene(r_method,drug,drug_type_name):
    group_gene = pd.read_csv('DEGs/{}/{}_0.05.csv'.format(r_method, drug_type_name),index_col=0)
    group_gene_up = group_gene.index[group_gene['log2FoldChange']>0].values
    group_gene_up_symbol = list(ensembol_gene_symbol(group_gene_up).values())
    group_gene_down = group_gene.index[group_gene['log2FoldChange']<0].values
    group_gene_down_symbol = list(ensembol_gene_symbol(group_gene_down).values())
    return {drug:{'final_pcm_symbol':group_gene_up_symbol,
                           'final_ncm_symbol':group_gene_down_symbol}}


def get_common_gene_methods(drugs, data_saved, group_fdr_gene_dict, classic_gene_dict, overlap_list_str):
    common_gene_dict = {}
    overlap_list = overlap_list_str.split('_')
    for d in drugs:
        common_gene_dict_drug = {}
        for p_n_gene in ['final_pcm_symbol','final_ncm_symbol']:
            set_list_overlap = []
            set_list_overlap_merge = []
            if 'module' in overlap_list:
                a = set(data_saved[p_n_gene][d][0])
                a_merged = set(itertools.chain.from_iterable(data_saved[p_n_gene][d]))
                set_list_overlap.append(a)
                set_list_overlap_merge.append(a_merged)
            if 'clinical' in overlap_list:
                b = set(group_fdr_gene_dict[d][p_n_gene])
                set_list_overlap.append(b)
                set_list_overlap_merge.append(b)
            if 'singlecor' in overlap_list:
                c = set(classic_gene_dict[d][p_n_gene])
                set_list_overlap.append(c)
                set_list_overlap_merge.append(c)
            common_genes = set.intersection(*set_list_overlap)
            common_genes_merged = set.intersection(*set_list_overlap_merge)
            common_gene_dict_drug[p_n_gene] = {'all_module':common_genes_merged,'first_module':common_genes}
        common_gene_dict[d] = common_gene_dict_drug
    return common_gene_dict

# get the raw overlap without mudule, only clinical, single drug, raw correlation by distribution
def get_common_gene_methods_raw(drugs, data_saved, group_fdr_gene_dict, classic_gene_dict, overlap_list_str):
    common_gene_dict = {}
    overlap_list = overlap_list_str.split('_')
    chanaged_key = {'final_ncm_symbol':'dw_genes_symbol','final_pcm_symbol':'up_genes_symbol'}
    for d in drugs:
        common_gene_dict_drug = {}
        for p_n_gene in ['final_pcm_symbol','final_ncm_symbol']:
            set_list_overlap = []
            if 'module' in overlap_list:
                a = set(data_saved[chanaged_key[p_n_gene]][d])
                set_list_overlap.append(a)
            if 'clinical' in overlap_list:
                b = set(group_fdr_gene_dict[d][p_n_gene])
                set_list_overlap.append(b)
            if 'singlecor' in overlap_list:
                c = set(classic_gene_dict[d][p_n_gene])
                set_list_overlap.append(c)
            common_genes = set.intersection(*set_list_overlap)
            common_gene_dict_drug[p_n_gene] = {'all_module':common_genes,'first_module':common_genes}
        common_gene_dict[d] = common_gene_dict_drug
    return common_gene_dict


def save_as_text(r_method, drug_type_name, common_gene_dict,cor_type,overlap_list_str):
    # common_gene_dict = pickle.load(open('../results/common_genes_methods_dict','rb'))
    # Writing final modules
    with open('DEGs/{}/final_common_gene_{}_{}_{}.txt'.format(r_method, drug_type_name,cor_type,overlap_list_str), 'w') as o:
        for drug, gene_dict in common_gene_dict.items():
            for module_type, module_id_dict in gene_dict.items():
                for module_id, genes in module_id_dict.items():
                    o.write('{} \t {} \t {} \t'.format(drug, module_type, module_id) + '\t'.join(genes) + '\n')


def get_common_genes_one(r_method,drug,drug_type_name,cor_type,overlap_list_str):
    drugs  = ['Bortezomib',
              'Carfilzomib',
              'Ixazomib',
              'Oprozomib']
    classic_gene_dict = get_single_drug_fdr(drugs)
    group_fdr_gene_dict = get_group_fdr_gene(r_method,drug,drug_type_name)
    data_saved = pickle.load(open('results/result.out', 'rb'))
    if cor_type == 'raw':
        common_gene_dict = get_common_gene_methods_raw([drug], data_saved, group_fdr_gene_dict, classic_gene_dict,overlap_list_str)
    elif cor_type == 'after_module':
        common_gene_dict = get_common_gene_methods([drug], data_saved, group_fdr_gene_dict, classic_gene_dict,overlap_list_str)
    with(open('DEGs/{}/common_genes_methods_dict_{}_{}_{}'.format(r_method,drug_type_name,cor_type,overlap_list_str), 'wb')) as handle:
        pickle.dump(common_gene_dict, handle)
    save_as_text(r_method, drug_type_name, common_gene_dict,cor_type,overlap_list_str)