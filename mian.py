
import sys
sys.path.insert(0, 'src')
from tqdm import tqdm
from scipy import stats
import scripts as src
#src.make_output_folders()
import sys
import os
import copy as cp
import random
import numpy as np
import pandas as pd
import subprocess
import codecs
from scipy import stats
from collections import Counter
from tqdm import tqdm
import networkx as nx
from src.diamond import DIAMOnD
import multiprocessing as mp
from src.scripts import prepare_drug_gene
from DEGs.files_preparation import prepare_groups_by_commonpass,map_DEG, prepare_groups_by_commonpass_processed
from DEGs.get_common_genes import get_common_genes_one
from clusters.cluster_patients_by_genes import prepare_pateints,cluster_pateints
from ml_model import ml
from survival import correlation_g_p, surve_kaplanmeier
from geo import geo_analysis_and_plot,geo_read_samples, geo_run_logrank
from correlation import drug_response_correlation
from commonpass import commpass_run_logrank, commpass_plot_kmf

# drugs = drug_df.columns
# for direction in ['PCM','NCM']:
#     for drug in tqdm(drugs,desc='%s'%direction):
#         output = 'results/hotnet/%s/%s/'%(direction,drug)
#         src.run_iteratively_hotnet('results/hotnet_input/%s/%s.tsv'%(direction,drug),output)
#

#
# for direction in ['PCM','NCM']:
#     for drug in ['Bortezomib']:
#         output = 'results/hotnet/%s/%s/'%(direction,drug)
#         src.run_iteratively_hotnet('results/hotnet_input/%s/%s.tsv'%(direction,drug),output)


def run_one(direction_drug_list):
    direction, drug = direction_drug_list[0],direction_drug_list[1]
    input = 'results/hotnet_input/%s/%s.tsv' % (direction, drug)
    output = 'results/hotnet/%s/%s/' % (direction, drug)

    src.run_iteratively_hotnet(input, output)



def test_one(drug_df):
    drugs = drug_df.columns
    direction = 'PCM'
    drug = drugs[3]
    input_file_path = 'results/hotnet_input/%s/%s.tsv'%(direction,drug)
    output_path = 'results/hotnet/%s/%s/'%(direction,drug)
    p=100
    min_size=5
    n_iter=10

    drug = input_file_path.split('/')[-1][:-4]

    # Check if there are score for the given file
    my_input = []
    with open(input_file_path, 'r') as f:
        for line in f:
            hit = line.rstrip().split('\t')
            my_input.append(hit)

    my_genes = set([x[0] for x in my_input])

    network = 'data/string/'

    it = 1

    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    if not os.path.isdir(output_path + '1'):
        os.mkdir(output_path + '1')

    outfolder = output_path + "1"
    cmd = "python src/run_hotnet.py -p %s -n %s -s %s -o %s" % (p,network, input_file_path, outfolder)

    subprocess.Popen(cmd, shell=True).wait()



# if __name__ == '__main__':
#     pool = mp.Pool(4)
#     drug_selected = ['Bortezomib',
#                      'Carfilzomib',
#                      'Ixazomib',
#                      'Oprozomib']
#     direction_drug_list = []
#     for direction in ['PCM', 'NCM']:
#         for drug in tqdm(drug_selected, desc='%s' % direction):
#             direction_drug_list.append([direction, drug])
#
#     funclist = []
#     for d_d in direction_drug_list:
#         f = pool.apply_async(run_one,[d_d])
#         funclist.append(f)
#
#     for p in funclist:
#         p.get()


def compare_methods():
    cor_df_old = pd.read_csv('results/correlation_old.csv', index_col=0)
    cor_df_new = pd.read_csv('results/cor_pd_z.csv', index_col=0).T
    cor_df_yin= pd.read_csv('results/cor_pd.csv', index_col=0).T
    drug_selected = ['Bortezomib',
                     'Carfilzomib',
                     'Ixazomib',
                     'Navitoclax',
                     'Oprozomib',
                     'Venetoclax']
    cor_df_new = cor_df_new[drug_selected]
    cor_df_yin = cor_df_yin[drug_selected]
    cor_df_old =  cor_df_old.sort_values(by=['Bortezomib'])
    cor_df_new = cor_df_new.sort_values(by=['Bortezomib'])
    cor_df_yin = cor_df_yin.sort_values(by=['Bortezomib'])

def check_spearman():
    gene_data = pd.read_csv('data/MM_log2RPKM_expression_123_20200918.txt', sep='\t')
    patient_drug_data = pd.read_csv('data/DrugSensitivity.tsv', sep='\t')
    patient_drug_data_2 = pd.read_csv('data/MMPI_DrugSensitivity.tsv', sep='\t')
    expr_df, drug_df = prepare_drug_gene(gene_data, patient_drug_data)
    drug_selected = ['Bortezomib',
                     'Carfilzomib',
                     'Ixazomib',
                     'Navitoclax',
                     'Oprozomib',
                     'Venetoclax']
    drug_df = drug_df[drug_selected ]
    expr_df_2, drug_df_2 = prepare_drug_gene(gene_data, patient_drug_data_2)
    drug = 'Bortezomib'
    gene = 'ENSG00000205542'

    gene = 'ENSG00000163220'

    def cor_get(drug, gene, drug_df, expr_df):
        drug_vector = drug_df[drug]

        # removing nan cls
        nan_cl = pd.isnull(drug_vector)
        drug_vector = list(drug_vector[~nan_cl])
        X_expr = expr_df.loc[~nan_cl]

        gene_vector = X_expr[gene]
        # Pearson correlation
        r, pval = stats.spearmanr(gene_vector, drug_vector)

        return r, pval

    drug = 'Bortezomib'
    gene_2 = 'ENSG00000000419'
    cor_get(drug, gene_2, drug_df, expr_df)
    cor_get(drug, gene_2, drug_df_2, expr_df_2)

    drug = 'Oprozomib'

    drug = 'ENSG0000016069'

    df = pd.merge(expr_df_2, drug_df_2, left_index=True, right_index=True, how='inner')
    rho, pval = stats.spearmanr(df[drug], df[gene], nan_policy='omit')


def validation_gene():
    # get the deg pre_files from gene and annotation


    table_dict = {'bor_all_line': {'table_name': 'bor_all_line', 'drug_name': 'Bortezomib'},
                  'bor_dex_1_line': {'table_name': 'bor_dex_1_line', 'drug_name': 'Bortezomib'},
                  'bor_based_1_line': {'table_name':  'bor_based_1_line', 'drug_name': 'Bortezomib'},
                  'carf_all_line': {'table_name': 'carf_all_line', 'drug_name': 'Carfilzomib'},
                  'carf_base_all_line': {'table_name': 'carf_base_all_line', 'drug_name': 'Carfilzomib'}}

    def generate_pateint_file(table_dict):
        dataset_gene_count = pd.read_csv('data/MMRF_CoMMpass_IA15a_E74GTF_HtSeq_Gene_Counts.txt', sep='\t', index_col=0)
        for drug_type_name, values in table_dict.items():
            # prepare_groups_by_commonpass(table_dict[drug_type_name]['table_name'], ['sCR', 'CR', 'VGPR'], ['SD', 'PD'],
            #                              drug_type_name=drug_type_name,
            #                              dataset_gene=dataset_gene_count,
            #                              type='DEG')

            prepare_groups_by_commonpass_processed(table_dict[drug_type_name]['table_name'], dataset_gene_count,'DEG')

    def run_r(drug_type_name):
        # Define command and arguments
        command = 'Rscript'
        command = 'C://Program Files//R//R-3.5.2//bin//Rscript'
        command = '../../../../../Program Files/R/R-3.5.2/bin/Rscript DEGS/dge_covid.R'
        path2script = 'DEGs/dge_covid.R'

        # Variable number of args in a list
        args = [drug_type_name]

        # Build subprocess command
        cmd = [command, path2script] + args

        # check_output will run the command and store to result
        x = subprocess.check_output(cmd, universal_newlines=True)
        print('edge R done')

    # get the common genes
    def get_common_gene_final(table_dict):
        for r_method in ['deseq2','edgeR']:
            for drug_type_name, values in table_dict.items():
                for cor_type in ['raw', 'after_module']:
                    for overlap_list_str in ['module_clinical_singlecor',
                                             'clinical_singlecor',
                                             'module_clinical',
                                             'module_singlecor']:
                        drug = table_dict[drug_type_name]['drug_name']
                        #map_DEG('edgeR', drug_type_name)
                        get_common_genes_one(r_method, drug, drug_type_name,cor_type,overlap_list_str)

    get_common_gene_final(table_dict)


def cluster_patentis():
    dataset_gene_fpkm = pd.read_csv('data/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM_2.txt',
                                    sep='\t',
                                    index_col=0).T

    # drug_type_name = 'bortezomib_1'
    for drug_type_name in ['bortezomib_1', 'bortezomib_2', 'bortezomib_3','carf']:
        path_2 = 'results/'
        for group_type in ['two_group','sub_group']:
            for module_type in [0,1]:
                dataset_gene_count_selected, group_dict = prepare_pateints(dataset_gene_fpkm,drug_type_name, path_2,module_type, group_type)
                cluster_pateints(dataset_gene_count_selected, group_dict, drug_type_name, module_type, group_type)


if __name__ == "__main__":
    # validation_gene()
    # cluster_patentis()
    # ml.main()
    # correlation_g_p.main()
    # surve_kaplanmeier.main()
    # analysis_and_plot.main()
    # drug_response_correlation.main()
    # geo_read_samples.main()
    # geo_run_logrank.main()
    # geo_analysis_and_plot.main()
    # commpass_run_logrank.main()
    # commpass_plot_kmf.main()
    ml.main()

