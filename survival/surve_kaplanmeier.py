import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines import NelsonAalenFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from collections import defaultdict
import pickle
import tqdm
import multiprocessing as mp
import matplotlib.gridspec as gridspec
from commonpass.commpass_data_process import ensembol_gene_symbol_pd_folder

# lifetimes :https://lifelines.readthedocs.io/en/latest/

# can use 'censos' 'ttcos' or censpfs, ttcpfs

# log_rank test by two dataset from to group with censos and time
def cal_log_rank(data_h, data_l, ttc, cens):
    # Define variables :
    T_h = data_h[ttc]
    E_h = data_h[cens]
    T_l = data_l[ttc]
    E_l = data_l[cens]

    results = logrank_test(T_h,
                           T_l,
                           event_observed_A=E_h,
                           event_observed_B=E_l)
    return results


# calculate the KaplanMeier
def cal_kmf(data_h, data_l, ttc, cens):
    # The 1st arg accepts an array or pd.Series of individual survival times
    # The 2nd arg accepts an array or pd.Series that indicates if the event
    # interest (or death) occured.
    # kmf_m for male data.
    # kmf_f for female data.
    kmf_h = KaplanMeierFitter()
    kmf_l = KaplanMeierFitter()

    kmf_h.fit(durations=data_h[ttc],
              event_observed=data_h[cens],
              label="higher expression (n={})".format( data_h.shape[0]))
    kmf_l.fit(durations=data_l[ttc],
              event_observed=data_l[cens],
              label="lower expression (n={})".format( data_l.shape[0]))
    return kmf_h, kmf_l


# Cox proportional hazard model:
'''
The Cox proportional hazard model is basically a regression model generally 
used by medical researchers to find out the relationship between the survival 
time of a subject and one or more predictor variables.
'''
def cal_cox_proprtion_hazard(data, varible_col, ttc, cens):
    # Cox regression :
    data = data[varible_col]
    # Create an object :
    cph = CoxPHFitter()
    cph.fit(data, ttc, event_col=cens)
    return cph

def prepare_groups_median_way(data, gene):
    median_value = np.median(data[gene])
    data_h = data[data[gene] >= median_value]
    data_l = data.loc[data[gene] < median_value]
    return data_h, data_l

# sort by gene value for furture selection
def prepare_groups_pertage_way(data, gene,  cut_percent):
    # Organize our data
    # If gene expression larger than median, then group = higher
    # If gene expression larger than median, then group = lower
    data = data.sort_values(by=gene, ascending=False)
    data_h = np.split(data, [int(cut_percent*len(data))])[0]
    data_l = np.split(data, [int((1-cut_percent)*len(data))])[1]
    data_h['group'] = 'higher'
    data_l['group'] = 'lower'

    return data_h, data_l

# prepare the data base on selected gene(varibale) into two dataframe
def prepare_groups_for_km(gene, data, ttc, cens, cut_percent):
    data = data[[gene, cens, ttc]]
    data = data.dropna(axis=0)
    data_h, data_l = prepare_groups_pertage_way(data, gene,  cut_percent)

    return data_h, data_l


# compare group
def compare_groups(gene_1, ttc, cens, g, pateint_group_dict, cut_percent, anno_gene_dict, data_soure, path_saved):
    expr_df = pateint_group_dict[g]['gene_data']
    time_df = pateint_group_dict[g]['patient_data']
    data = pd.merge(time_df,
                    expr_df,
                    left_index=True,
                    right_index=True,
                    how='inner')
    # prepare groups
    data_h, data_l = prepare_groups_for_km(gene_1, data, ttc, cens, cut_percent)

    gene = anno_gene_dict[gene_1]
    gene = gene.split(' /// ')[0]
    # cal the kmf
    kmf_h, kmf_l = cal_kmf(data_h, data_l, ttc, cens)
    # kmf_h_perform = kmf_h.event_table
    # kmf_l_perform = kmf_l.event_table
    # kmf_h_func = kmf_h.survival_function_
    # kmf_l_func = kmf_l.survival_function_

    # Log-Rank Test
    results = cal_log_rank(data_h, data_l, ttc, cens)

    # Plot the survival_function data :
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(4, 1)
    ax1 = fig.add_subplot(gs[0:2, :])
    kmf_h.plot_survival_function(ax=ax1)
    kmf_l.plot_survival_function(ax=ax1)
    ax1.set_xlabel("Times of {}".format(ttc))
    ax1.set_ylabel(" {} (%)".format(cens))
    ax1.set_title("KMF of {} {} from {} cut off {}".format(gene, g, data_soure, cut_percent))

    # Add a table at the bottom of the axes

    merged = pd.concat([data_h, data_l], axis=0)

    merged_pd = merged.groupby(['group']).agg({ttc: "mean", gene_1: 'mean', 'group': 'count'})
    merged_pd = merged_pd.rename(columns={ttc: "Mean {}".format(ttc),
                                          gene_1: 'Mean {}'.format(gene),
                                          'group': 'The number of group'})
    merged_pd = merged_pd.apply(lambda x: round(x, 2))
    ax2 = fig.add_subplot(gs[2, :])
    ax2.axis('off')
    ax2.axis('tight')
    table = ax2.table(cellText=merged_pd.values,
                          rowLabels=list(merged_pd.index),
                          colLabels=list(merged_pd.columns),
                         cellLoc = 'center',
                          loc='center')
    table.set_fontsize(10)
    table.scale(1.2, 1.2)


    ax3 = fig.add_subplot(gs[3, :])
    ax3.axis('off')
    ax3.axis('tight')
    result_statistic = results.summary
    table_2 = ax3.table(cellText=result_statistic.values,
                 rowLabels=['logrank test'],
                 colLabels=list(result_statistic.columns),
                 cellLoc='center',
                 loc='center')
    table_2.set_fontsize(10)
    table_2.scale(1.1, 1.11)
    plt.tight_layout()
    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.2, bottom=0.2)

    # save out
    if not os.path.exists(path_saved):
        os.makedirs(path_saved)

    plt.savefig('{}/{}_{}_{}_{}_KMF_curve.png'.format(path_saved, g, gene, cens, cut_percent))

    # plot density
    # ax2 = plt.subplot(111)
    # kmf_h.plot_cumulative_density(ax=ax2)
    # kmf_l.plot_cumulative_density(ax=ax2)
    # plt.xlabel("Times passed")
    # plt.ylabel("Cumulative dead rate")
    # plt.title("Cumulative dead rate of gene {}".format(gene))
    # plt.savefig('{}/{}_dead rate.png'.format(path_saved, gene))
    #
    #
    # # Log-Rank Test
    # results.summary.to_csv('results/commpass/logrank_test by {}.csv'.format(gene))


# calcuLATE ALL LOGRANK TEST
def cal_all_gene_logrank(data, ttc, cens, cut_percent, col_genes):
    # group gene into tw groups

    logrank_result = defaultdict()

    # loop the genes, get the

    for gene in tqdm.tqdm(col_genes):
        new_data = data[[gene,
                       cens,
                       ttc]]

        data_h, data_l = prepare_groups_pertage_way( new_data, gene, cut_percent)
        # kmf_h, kmf_l = cal_kmf(data_h, data_l, ttc, cens)
        # Log-Rank Test
        results = cal_log_rank(data_h, data_l,  ttc, cens)
        logrank_result[gene] = {'log_rank_result': results}
        # logrank_result[gene] = {'higher_expression_patients': data_h,
        #                         'lower_expression_patients': data_l,
        #                         'higher_expression_kmf': kmf_h,
        #                         'lower_expression_kmf': kmf_l,
        #                         'log_rank_result': results}

    # prepare to csv
    result_pd = pd.concat([v['log_rank_result'].summary for k,v in logrank_result.items()],axis=0)
    result_pd['genes'] = list(logrank_result.keys())
    result_pd.columns = [i + '_logrank' for i in result_pd.columns]

    return logrank_result, result_pd

def all_gene_logrank_mutiple( pateint_group_dict, cens, ttc, cut_percent, path_saved):

    for g, v_1 in pateint_group_dict.items():
        expr_df = v_1['gene_data']
        time_df = v_1['patient_data']
        data = pd.merge(time_df,
                        expr_df,
                        left_index=True,
                        right_index=True,
                        how='inner')

        col_genes = list(expr_df.columns)
        logrank_result, logrank_pd = cal_all_gene_logrank(data, ttc, cens, cut_percent,  col_genes )

        # pateint_group_dict[g].update({'logrank_result_{}_{}'.format(cens, cut_percent*100): logrank_result,
        #                               'logrank_pd_{}_{}'.format(cens, cut_percent*100): logrank_pd})
        if not os.path.exists(path_saved):
            os.makedirs(path_saved)
        logrank_pd.to_csv(path_saved + 'logrank&{}&{}_{}.csv'.format(g, cens, cut_percent*100), index=None)
    return pateint_group_dict
