import pickle
from survival.correlation_g_p import correlation_with_fdr_mutiple
from geo.geo_read_samples import ensembol_gene_symbol_pd_folder_geo

def get_correlation():
    pateint_group_dict = pickle.load(open('results/geo/subgroup/pateint_group_dict', 'rb'))
    for cens_term in ['os', 'pfs']:
        cens = 'cens{}'.format(cens_term)
        if cens_term == 'pfs':
            tt ='ttc{}'.format(cens_term)
        elif cens_term == 'os':
            tt = 'tt{}'.format(cens_term)
        path_saved = 'results/geo/correlation/'
        pateint_group_dict = correlation_with_fdr_mutiple(pateint_group_dict, cens, tt, path_saved)


def main():

    get_correlation()

    path_selected = '../results/geo/correlation/'
    col_gene = 'gene'
    ensembol_gene_symbol_pd_folder_geo(path_selected, col_gene)