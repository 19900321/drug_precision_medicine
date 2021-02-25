import pickle
import os
import pandas as pd
import multiprocessing as mp
from survival.surve_kaplanmeier import all_gene_logrank_mutiple, cal_all_gene_logrank
from commonpass.commpass_data_process import ensembol_gene_symbol_pd_folder

def all_gene_logrank_mutiple_pipeline(g_cens_list):
    cut_percent = 0.33
    pateint_group_dict = pickle.load(open('../results/commpass/subgroup/pateint_group_log_dict', 'rb'))
    #pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_dict_cnv_s100', 'rb'))
    g, cens_term = g_cens_list[0], g_cens_list[1]
    pateint_group_dict_one = {g: pateint_group_dict[g]}
    cens = 'cens{}'.format(cens_term)
    if cens_term == 'pfs':
        tt = 'ttc{}'.format(cens_term)
    elif cens_term == 'os':
        tt = 'tt{}'.format(cens_term)

    path_saved = '../results/commpass/logrank/'
    all_gene_logrank_mutiple(pateint_group_dict_one, cens, tt, cut_percent, path_saved)

    # with(open('results/commpass/subgroup/pateint_group_dict', 'wb')) as handle:
    #         pickle.dump(pateint_group_dict, handle)


def map_symbol_logrank():
    path_saved = 'results/commpass/logrank/'
    col_gene = 'genes_logrank'
    ensembol_gene_symbol_pd_folder(path_saved, col_gene)


def only_s100_logrank():
    cut_percent_s = [0.25, 0.33, 0.5]
    #pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_log_dict', 'rb'))
    pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_dict_cnv_s100', 'rb'))
    path_saved = 'results/commpass/cnv/cnv_logrank_S100/'
    g_list = list(pateint_group_dict.keys())
    cens_list = ['pfs', 'os']
    col_genes = ['ENSG00000116497',
              'ENSG00000143546',
              'ENSG00000160678',
              'ENSG00000163221',
              'ENSG00000197747',
              'ENSG00000163191',
              'ENSG00000197956',
              'ENSG00000163220',
              'ENSG00000196154',
              'ENSG00000189171']

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

    anno_gene_dict = pickle.load(open('results/commpass/ensembol_symbol_dict_commpass', 'rb'))
    col_genes = [k for k, v in anno_gene_dict.items() if v in genes]

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
    ensembol_gene_symbol_pd_folder(path_saved, col_gene)


#def run_main():
# if __name__ == '__main__':
#
#     pool = mp.Pool(2)
#     funclist = []
#     # g_list = ['Bor_Dex', 'Bor_based', 'Carfilzomib_based']
#     # cens_list = ['pfs']
#     g_list = ['all_patients']
#     cens_list = ['pfs']
#     g_cens_list = [[g, cens] for g in g_list for cens in cens_list]
#
#     for one_case in g_cens_list:
#         f = pool.apply_async(all_gene_logrank_mutiple_pipeline,[one_case])
#         funclist.append(f)
#
#     for p in funclist:
#         p.get()


def main():
    only_s100_logrank()