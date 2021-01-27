import pickle
from survival.surve_kaplanmeier import compare_groups

def plot_survival_figure():
    anno_gene_dict = pickle.load(open('results/geo/id_symbol_dict_geo', 'rb'))

    pateint_group_dict = pickle.load(open('results/geo/subgroup/pateint_group_dict', 'rb'))

    g_list = ['PS341', 'Dex']
    cens_list = ['pfs', 'os']
    cut_percent_s = [0.25, 0.33, 0.5]
    gene_s = ['202598_at',
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

    path_saved = 'results/geo/geo_group_kmf_S100/'
    data_soure = 'APEX'
    for gene in gene_s:
        for g in g_list:
            for cens_term in cens_list:
                cens = 'cens{}'.format(cens_term)
                if cens_term == 'pfs':
                    tt = 'ttc{}'.format(cens_term)
                elif cens_term == 'os':
                    tt = 'tt{}'.format(cens_term)
                for cut_percent in cut_percent_s:
                    compare_groups(gene, tt, cens, g, pateint_group_dict, cut_percent, anno_gene_dict, data_soure, path_saved)


def main():
    plot_survival_figure()

