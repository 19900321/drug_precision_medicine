import pickle
from survival.surve_kaplanmeier import compare_groups

def plot_kmf_groups_s100():
    # map_symbol_logrank()
    pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_log_dict', 'rb'))

    #group_s = list(pateint_group_dict.keys())
    group_s = ['all_patients']
    ttc = 'ttcpfs'
    cens = 'censpfs'

    cut_percent_s = [0.50, 0.33]
    gene_s = ['ENSG00000116497',
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
    path_saved = 'results/commpass/commpass_group_kmf/'
    data_soure = 'CoMMpass'
    for gene in gene_s:
        for g in group_s:
            for cut_percent in cut_percent_s:
                compare_groups(gene, ttc, cens, g, pateint_group_dict, cut_percent, anno_gene_dict, data_soure, path_saved)


def main():
    plot_kmf_groups_s100()