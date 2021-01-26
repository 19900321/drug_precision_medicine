import pickle
from survival.surve_kaplanmeier import compare_groups

def plot_kmf_groups_s100():

    cut_percent_s = [0.25, 0.33, 0.5]
    pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_log_dict', 'rb'))
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

    anno_gene_dict = pickle.load(open('results/commpass/ensembol_symbol_dict_commpass', 'rb'))
    path_saved = 'results/commpass/commpass_group_kmf_S100/'
    data_soure = 'CoMMpass'
    for gene in col_genes:
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
    plot_kmf_groups_s100()