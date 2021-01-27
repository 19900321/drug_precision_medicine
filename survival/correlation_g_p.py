import survival.dataset as dataset
import pandas as pd
from src.scripts import get_correlations_2
import pickle
from commonpass.commpass_data_process import ensembol_gene_symbol_pd_folder

def single_correlation(data):
    cor_pd = data.corrwith(data['survive_score'])
    cor_pd_2 = pd.DataFrame(cor_pd.sort_values(ascending=False), columns=['sur_cor']).reset_index()
    cor_pd_2 = cor_pd_2.rename(columns={'index':'varibles'})
    return cor_pd_2

def correlation_death_only(data):
    data = data[data['censos'] == 1]
    cor_pd = data.corrwith(data['survive_score'])
    cor_pd_2 = pd.DataFrame(cor_pd.sort_values(ascending=False), columns=['sur_cor']).reset_index()
    cor_pd_2 = cor_pd_2.rename(columns={'index': 'varibles'})
    return cor_pd_2

def correlation_with_fdr(expr_df, drug_df, cens, tt):
    drug_df = drug_df.loc[drug_df[cens]==1, [tt]]
    expr_df = expr_df.loc[expr_df.index.isin(drug_df.index), :]
    cor_pd, cor_pd_stack,pvalue_pd,fdr_pd = get_correlations_2(expr_df, drug_df)
    return cor_pd, cor_pd_stack


def correlation_with_fdr_mutiple( pateint_group_dict, cens, tt, path_saved):

    for g, v_1 in pateint_group_dict.items():
        expr_df = v_1['gene_data']
        time_df = v_1['patient_data']
        cor_pd, cor_pd_stack = correlation_with_fdr(expr_df, time_df, cens, tt)
        pateint_group_dict[g].update({'cor_pd_{}'.format(cens): cor_pd,
                                      'cor_pd_stack_{}'.format(cens): cor_pd_stack})
        cor_pd_stack.to_csv(path_saved + 'correlation&{}&{}.csv'.format(g, cens), index=None)
    return pateint_group_dict

def map_symbol_cor():
    path_saved = 'results/commpass/correlation/'
    col_gene = 'gene'
    ensembol_gene_symbol_pd_folder(path_saved, col_gene)


def main():

    def get_cor():
        pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_log_dict', 'rb'))
        for cens_term in ['os', 'pfs']:
            cens = 'cens{}'.format(cens_term)
            if cens_term == 'pfs':
                tt ='ttc{}'.format(cens_term)
            elif cens_term == 'os':
                tt = 'tt{}'.format(cens_term)
            path_saved = 'results/commpass/correlation/'
            pateint_group_dict = correlation_with_fdr_mutiple(pateint_group_dict, cens, tt, path_saved)
        with(open('results/commpass/subgroup/pateint_group_dict', 'wb')) as handle:
                pickle.dump(pateint_group_dict, handle)



