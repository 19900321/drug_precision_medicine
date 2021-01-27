from survival.survive_transfer import cal_survive_score
import survival.dataset as dataset
import pandas as pd
from survival.correlation_g_p import single_correlation, correlation_death_only
from survival.mutiple_linear import train_linear
from survival.surve_kaplanmeier import compare_groups
from survival.surve_kaplanmeier import cal_cox_proprtion_hazard
# from neural_net import cal_basic_model
from survival.surve_kaplanmeier import cal_all_gene_logrank
from survival import get_gene_name
import pickle

# tranfer the censos and ttcos to a score
def get_survive_score():
    data = dataset.load_survive_patients()
    data['survive_score'] = data.apply(lambda x: cal_survive_score(data, x.censos, x.ttcos), axis=1)
    data.to_csv('tmp/_survive_score.csv')


# merge the gene expression with survive score by patients id
def get_merged_result(save_path):
    data_survive = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15_STAND_ALONE_SURVIVAL_TRTRESP_merged.csv',
                               index_col='gene_patient_id')
    data_gene = pd.read_csv('data/CoMMpass_IA15_FlatFiles/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM_processed.txt',
                                   sep='\t', index = None)
    data = pd.merge(data_survive,
                    data_gene,
                    left_index=True,
                    right_index=True,
                    how='inner')
    data.to_csv('results/{}/_merged_gene_ttc.csv'.format(save_path))


# calculate the correlation of genes and survive score among all the patients
def get_correlation():
    data = pd.read_csv('tmp/_merged_gene_sur.csv', index_col=0)
    cor_pd = single_correlation(data)
    cor_pd_2 = correlation_death_only(data)
    cor_all = pd.merge(cor_pd, cor_pd_2, left_on='gene', right_on='gene', how='inner')
    cor_all.to_csv('tmp/_correlation_gene_sur.csv', index=None)


# built mulitple linear regression
def built_linear():
    data = pd.read_csv('tmp/_merged_gene_sur.csv', index_col=0)
    cor_pd = pd.read_csv('tmp/_correlation_gene_sur.csv')
    col_genes, regr = train_linear(data)
    gene_dict = dict(zip(col_genes, list(regr.coef_[0]*100000)))
    cor_pd['multiple_coef'] = cor_pd['gene'].apply(lambda x: gene_dict[x] if x in gene_dict else None)
    cor_pd.to_csv('tmp/cor__multi_gene_sur.csv', index=None)

# use logrank for each gene
def cal_logrank_genes(save_path, ttc, cens):
    data = pd.read_csv('results/{}/_merged_gene_ttc.csv'.format(save_path), index_col=0)
    col_genes = [c for c in data.columns if c.startswith('ENSG')] # too larger
    # col_genes_pd = pd.read_csv('tmp/cor__multi_gene_sur.csv')
    # col_genes_pd = col_genes_pd[(col_genes_pd['p_value'] < 0.01) & ((abs(col_genes_pd['coef'])) > 0.3)]
    # col_genes = list(col_genes_pd['gene'])
    varibale_cols = [ttc, cens] + col_genes
    data = data[varibale_cols]
    logrank_result, result_pd = cal_all_gene_logrank(data)
    # save result as dictionary
    #with open('tmp/logrank_result', 'wb') as handle:
        #pickle.dump(logrank_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    result_pd.to_csv('tmp/logrank_result.csv')
    # merge with cor
    cor_pd = pd.read_csv('tmp/cor__multi_gene_sur.csv')
    cor_new = pd.merge([cor_pd, result_pd], left_on='gene', right_on='gene_logrank', how = 'left')
    cor_new.to_csv('tmp/cor__multi_gene_sur_logrank.csv')

# test gene in interest
def group_timeline_comparison_one():
    #gene = 'ENSG00000168497'# SDPR gene
    gene = 'ENSG00000163736' #PPBP
    file_path = "tmp/_merged_gene_sur.csv"
    compare_groups(gene, file_path)

# result from cox model
def cal_cox_multiple():
    data = pd.read_csv('tmp/_merged_gene_sur.csv', index_col=0)
    #col_genes = [c for c in data.columns if c.startswith('ENSG')] # too larger
    col_genes_pd = pd.read_csv('tmp/cor__multi_gene_sur.csv')
    col_genes_pd = col_genes_pd[(col_genes_pd['p_value'] < 0.01) & ((abs(col_genes_pd['coef'])) > 0.25)]
    col_genes = list(col_genes_pd['gene'])
    varibale_cols = ["ttcos", "censos"] + col_genes

    cox_m = cal_cox_proprtion_hazard(data, varibale_cols)
    cox_m.summary.to_csv('tmp/gene_varibale.csv')
    return cox_m


# use the neural model
def cal_net_model():
    data = pd.read_csv('tmp/_merged_gene_sur.csv')
    data = data[data['censos'] == 0]
    col_genes = [c for c in data.columns if c.startswith('ENSG')]
    X, Y = data[col_genes], data['ttcos']
    #cal_basic_model(X, Y)


def get_engs_gene_name():
    cor_pd = pd.read_csv('tmp/cor__multi_gene_sur_logrank.csv')
    cox_pd = pd.read_csv('tmp/gene_varibale.csv')
    col_genes = list(cor_pd['gene'])
    sym_list, sym_dict = get_gene_name(col_genes)
    cor_pd['gene_name'] = cor_pd['gene'].apply(lambda x: sym_dict[x] if x in sym_dict else None)
    cox_pd['gene_name'] = cox_pd['covariate'].apply(lambda x: sym_dict[x] if x in sym_dict else None)

    cor_pd.to_csv('results/commpass/cor__multi_gene_sur_logrank.csv')
    cox_pd.to_csv('results/commpass/gene_varibale.csv')


if __name__ == '__main__':
    # get_survive_score()
    # get_merged_result()
    # get_correlation()
    # built_linear()
    # group_timeline_comparison_one()
    # cox_m = cal_cox_multiple()
    # cal_net_model()
    cal_logrank_genes()
