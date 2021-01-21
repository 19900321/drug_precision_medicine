import pickle
import pandas as pd
import multiprocessing as mp
from survival.surve_kaplanmeier import all_gene_logrank_mutiple
from commonpass.data_process import ensembol_gene_symbol_pd_folder

def all_gene_logrank_mutiple_pipeline(g_cens_list):
    cut_percent = 0.33
    pateint_group_dict = pickle.load(open('../results/commpass/subgroup/pateint_group_log_dict', 'rb'))
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



#def run_main():
if __name__ == '__main__':

    pool = mp.Pool(2)
    funclist = []
    # g_list = ['Bor_Dex', 'Bor_based', 'Carfilzomib_based']
    # cens_list = ['pfs']
    g_list = ['all_patients']
    cens_list = ['pfs']
    g_cens_list = [[g, cens] for g in g_list for cens in cens_list]

    for one_case in g_cens_list:
        f = pool.apply_async(all_gene_logrank_mutiple_pipeline,[one_case])
        funclist.append(f)

    for p in funclist:
        p.get()
