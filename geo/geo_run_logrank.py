import pickle
import os
import pandas as pd
import multiprocessing as mp
from survival.surve_kaplanmeier import all_gene_logrank_mutiple, cal_all_gene_logrank
from geo.geo_read_samples import ensembol_gene_symbol_pd_folder_geo


def all_gene_logrank_mutiple_pipeline(g_cens_list):
    cut_percent = 0.33
    pateint_group_dict = pickle.load(open('../results/geo/subgroup/pateint_group_dict', 'rb'))
    g, cens_term = g_cens_list[0],g_cens_list[1]
    pateint_group_dict_one = {g: pateint_group_dict[g]}
    cens = 'cens{}'.format(cens_term)
    if cens_term == 'pfs':
        tt = 'ttc{}'.format(cens_term)
    elif cens_term == 'os':
        tt = 'tt{}'.format(cens_term)

    path_saved = '../results/geo/logrank/'
    all_gene_logrank_mutiple(pateint_group_dict_one, cens, tt, cut_percent, path_saved)

    # with(open('results/commpass/subgroup/pateint_group_dict', 'wb')) as handle:
    #         pickle.dump(pateint_group_dict, handle)


def only_s100_logrank():

    pateint_group_dict = pickle.load(open('results/geo/subgroup/pateint_group_dict', 'rb'))
    path_saved = 'results/geo/logrank_S100/'
    g_list = ['PS341', 'Dex']
    cens_list = ['pfs', 'os']
    cut_percent_s = [0.25, 0.33, 0.5]

    col_genes = ['202598_at',
              '217728_at',
              '205334_at',
              '202917_s_at',
              '213481_at',
              '203186_s_at',
              '208540_x_at',
              '203535_at',
              '218370_s_at',
              '200872_at',
              '205863_at',
              '206027_at',
              '204268_at',
              '218677_at',
              '200660_at']
    for g in g_list:
        expr_df = pateint_group_dict[g]['gene_data']
        time_df = pateint_group_dict[g]['patient_data']
        data = pd.merge(time_df,
                        expr_df,
                        left_index=True,
                        right_index=True,
                        how='inner')

        for cens_term in cens_list:
            cens = 'cens{}'.format(cens_term)
            if cens_term == 'pfs':
                tt = 'ttc{}'.format(cens_term)
            elif cens_term == 'os':
                tt = 'tt{}'.format(cens_term)
            for cut_percent in cut_percent_s:

                logrank_result, logrank_pd = cal_all_gene_logrank(data, tt, cens, cut_percent, col_genes)
                if not os.path.exists(path_saved):
                    os.makedirs(path_saved)
                logrank_pd.to_csv(path_saved + 'logrank&{}&{}_{}.csv'.format(g, cens_term, cut_percent * 100), index=None)
    col_gene = 'genes_logrank'
    ensembol_gene_symbol_pd_folder_geo(path_saved, col_gene)


# if __name__ == '__main__':
#
#     pool = mp.Pool(8)
#     funclist = []
#     g_list = ['PS341', 'Dex']
#     cens_list = ['pfs', 'os']
#     g_cens_list = [[g, cens] for g in g_list for cens in cens_list]
#
#     for one_case in g_cens_list:
#         f = pool.apply_async(all_gene_logrank_mutiple_pipeline,[one_case])
#         funclist.append(f)
#
#     for p in funclist:
#         p.get()
#
#     path_selected = '../results/geo/logrank/'
#     col_gene = 'genes_logrank'
#     ensembol_gene_symbol_pd_folder_geo(path_selected, col_gene)


def main():
    only_s100_logrank()