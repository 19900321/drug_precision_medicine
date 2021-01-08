# Check the versions of libraries
import numpy as np
import sys
import scipy
import numpy
import matplotlib.pyplot as plt
import pandas as pd
import sklearn
import pickle
import itertools
from scipy.stats import spearmanr
from operator import itemgetter, attrgetter
import mygene
from collections import defaultdict

# Load libraries
from pandas import read_csv
from pandas.plotting import scatter_matrix
from matplotlib import pyplot
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.metrics import mean_squared_error, r2_score
from src.scripts import prepare_drug_gene

from sklearn.datasets import load_boston
import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.datasets import make_friedman1
from sklearn.feature_selection import RFE
from sklearn.svm import SVR

def stepwise_selection(X, y,
                       initial_list=[],
                       threshold_in=0.01,
                       threshold_out = 0.05,
                       verbose=True):
    """ Perform a forward-backward feature selection
    based on p-value from statsmodels.api.OLS
    Arguments:
        X - pandas.DataFrame with candidate features
        y - list-like with the target
        initial_list - list of features to start with (column names of X)
        threshold_in - include a feature if its p-value < threshold_in
        threshold_out - exclude a feature if its p-value > threshold_out
        verbose - whether to print the sequence of inclusions and exclusions
    Returns: list of selected features
    Always set threshold_in < threshold_out to avoid infinite looping.
    See https://en.wikipedia.org/wiki/Stepwise_regression for the details
    """
    included = initial_list
    while True:
        changed=False
        # forward step
        excluded = list(set(X.columns)-set(included))
        new_pval = pd.Series(index=excluded)
        for new_column in excluded:
            model = sm.OLS(y, sm.add_constant(pd.DataFrame(X[included+[new_column]]))).fit()
            new_pval[new_column] = model.pvalues[new_column]
        best_pval = new_pval.min()
        if best_pval < threshold_in:
            best_feature = new_pval.argmin()
            included.append(best_feature)
            changed=True
            if verbose:
                print('Add  {:30} with p-value {:.6}'.format(best_feature, best_pval))

        # backward step
        model = sm.OLS(y, sm.add_constant(pd.DataFrame(X[included]))).fit()
        # use all coefs except intercept
        pvalues = model.pvalues.iloc[1:]
        worst_pval = pvalues.max() # null if pvalues is empty
        if worst_pval > threshold_out:
            changed=True
            worst_feature = pvalues.argmax()
            included.remove(worst_feature)
            if verbose:
                print('Drop {:30} with p-value {:.6}'.format(worst_feature, worst_pval))
        if not changed:
            break
    return included


def feature_selection(X, Y,columns_selected):
    estimator = SVR(kernel="linear")
    selector = RFE(estimator, n_features_to_select=8, step=4)
    selector = selector.fit(X, Y)
    feature_selected = np.array(columns_selected)[np.where(selector.support_ == True)]
    return feature_selected


mg = mygene.MyGeneInfo()
def uniprot_gene_symbol(uniprot_ids):
    uni_dict = {}
    for i in mg.querymany(uniprot_ids, scopes='uniprot', fields='symbol'):
        try:
            uni_dict.update({i['query']:i['symbol']})
        except:
            continue
    return uni_dict


def ensembol_gene_symbol(ensembol_ids):
    uni_dict = {}
    for i in mg.querymany(ensembol_ids, scopes='ensembl.gene', fields='symbol'):
        try:
            uni_dict.update({i['query']:i['symbol']})
        except:
            continue
    return uni_dict


# firstly prepare the drug data , gene data
def prepare_drug_gene_data():
    gene_data = pd.read_csv('data/MM_log2RPKM_expression_123_20200918.txt', sep='\t')
    patient_drug_data = pd.read_csv('data/DrugSensitivity.tsv', sep='\t')
    expr_df, drug_df = prepare_drug_gene(gene_data, patient_drug_data)
    expr_df.to_csv('data/gene_processed.csv')
    drug_df.to_csv('data/drug_processed.csv')


def get_final_columns_dict(drugs, data_saved):
    uni_sem_pd = read_csv('data/ensbl2AC.tsv', header=None, sep='\t')
    uni_sem_pd.columns = ['sem', 'uni']
    uni_sem_dict = dict(zip(uni_sem_pd['uni'], uni_sem_pd['sem']))
    final_columns_dict_d = {}
    for d in drugs:
        columns = data_saved['final_pcm'][d] + data_saved['final_ncm'][d]
        columns_1 = data_saved['final_pcm'][d][0] + data_saved['final_ncm'][d][0]

        if len(data_saved['final_pcm'][d]) >=2:
            columns_2_p = data_saved['final_pcm'][d][1]
        else:
            columns_2_p = []

        if len(data_saved['final_ncm'][d]) >= 2:
            columns_2_n = data_saved['final_ncm'][d][1]
        else:
            columns_2_n = []

        columns_2 = columns_2_p + columns_2_n

        genes = [uni_sem_dict[g] for g in set(itertools.chain.from_iterable(columns)) if g in uni_sem_dict]
        genes_1 = [uni_sem_dict[g] for g in columns_1 if g in uni_sem_dict]
        genes_2 = [uni_sem_dict[g] for g in columns_2 if g in uni_sem_dict]

        final_columns_dict_d[d] = [genes,genes_1,genes_2]

    return final_columns_dict_d


# DecisionTreeRegressor
def dmr(X_train, X_test, y_train, y_test):
    clf = DecisionTreeRegressor().fit(X_train, y_train)
    predicted = clf.predict(X_test)
    expected = y_test

    cor, p_value = spearmanr(round(expected, 2), predicted)
    print('spearman  coorelation {}, {} '.format(cor, p_value))
    # The mean squared error
    print('Mean squared error: %.2f'
          % mean_squared_error(expected, predicted))
    # The coefficient of determination: 1 is perfect prediction
    print('Coefficient of determination: %.2f'
          % r2_score(expected, predicted))

    # Plot outputs
    plt.scatter(expected, predicted)
    plt.plot(expected, predicted, color='blue', linewidth=3)
    plt.xticks(())
    plt.yticks(())
    plt.show()



def test_one_drug(dataset, drug, final_columns_dict,type_index,scale):

    def genrate_linear(X_train, X_test, y_train, y_test):
        clf = LinearRegression()
        clf.fit(X_train, y_train)
        predicted = clf.predict(X_test)
        expected = y_test

        cor, p_value = spearmanr(round(expected,2), predicted)

        print('spearman  correlation {}, {} '.format( cor, p_value))
        # The mean squared error
        MSE = mean_squared_error(expected, predicted)
        print('Mean squared error: %.2f'
              % MSE)
        # The coefficient of determination: 1 is perfect prediction
        coeffient = r2_score(expected, predicted)
        print('Coefficient of determination: %.2f'
              % coeffient)

        # Plot outputs
        plt.scatter(expected, predicted)
        plt.show()

        return clf, MSE,coeffient, cor, p_value


    def varibale_rank(linear_model,X_columns,X_symbol_dict ):
        # check the coef od each column
        gene_symbol= [X_symbol_dict[i] if i in X_symbol_dict else None for i in X_columns]
        gene_coef_s = list(zip(columns_selected, gene_symbol, linear_model.coef_[0]))
        gene_coef_s = sorted(gene_coef_s, key=itemgetter(2), reverse=True)
        return gene_coef_s


    def prepare_data(dataset, drug, columns_selected,scale):
        dataset = dataset[dataset[drug].notnull()]
        X = dataset[columns_selected]
        Y = dataset[[drug]]
        feature_selected = feature_selection(X, Y, columns_selected)
        X = X[feature_selected]
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.25, random_state=0)
        if scale == True:
            sc_X = StandardScaler()
            X_train = sc_X.fit_transform(X_train)
            X_test = sc_X.fit_transform(X_test)
        return X_train, X_test, y_train, y_test, feature_selected

    columns_selected = final_columns_dict[drug][type_index]
    X_symbol_dict = ensembol_gene_symbol(columns_selected)
    X_train, X_test, y_train, y_test,feature_selected = prepare_data(dataset, drug, columns_selected,scale)

    feature_selected_symbol = [X_symbol_dict[g] if g in X_symbol_dict else None for g in feature_selected ]
    print('the selected feature is'.format(','.join(feature_selected_symbol)))
    clf, MSE,coeffient, cor, p_value = genrate_linear(X_train, X_test, y_train, y_test)
    gene_coef_s = varibale_rank(clf,feature_selected, X_symbol_dict )
    return clf, MSE,coeffient, cor, p_value, gene_coef_s

def drugs_prediction(drugs, dataset, final_columns_dict,scale):
    result_dict_all = {}
    type_index = [0,1,2]
    for d in drugs:
        result_dict = {}
        for type_index_one in type_index:
            if len(final_columns_dict[d][type_index_one]) >5:
                clf, MSE,coeffient, cor, p_value, gene_coef_s = test_one_drug(dataset, d, final_columns_dict, type_index_one,scale)
                result_dict[type_index_one] = {'clf':clf,
                                  'MSE':MSE,
                                  'coeffient':coeffient,
                                  'cor':cor,
                                  'p_value':p_value,
                                  'gene_coef_s':gene_coef_s}
        result_dict_all[d] = result_dict
    return result_dict_all


def model_prediction(pathes, drugs, dataset, scale = True):
    for module_path in pathes:
        data_saved = pickle.load(open(module_path + 'result.out', 'rb'))
        final_columns_dict = get_final_columns_dict(drugs, data_saved)
        result_dict = drugs_prediction(drugs, dataset, final_columns_dict,scale)
        if scale ==True:
            with open(module_path + 'linear_reuslt_cut_scale','wb') as handle:
                pickle.dump(result_dict,handle)
        else:
            with open(module_path + 'linear_reuslt_cut','wb') as handle:
                pickle.dump(result_dict,handle)


# prepare ML result
def prepare_prediction_result(result_dict):
    result_pd = pd.concat([pd.DataFrame(d).T for d in list(result_dict.values())], keys = list(result_dict.keys()) ).reset_index()
    result_pd = result_pd.rename(columns={'level_0':'drug','level_1':'module'})
    result_pd['module'] = result_pd['module'].astype(str).replace('0','all_module').replace('1','module_1').replace('2','module_2')
    result_pd = result_pd.drop(columns=['clf'])

    def extend_gene_coef_s(result_pd):
        extend_coef_pd = []
        for index, rows in result_pd.iterrows():
            for record in rows['gene_coef_s']:
                extend_coef_pd.append([rows.drug, rows.module, rows.cor, rows.p_value] + list(record))

        extend_coef_pd = pd.DataFrame(extend_coef_pd, columns=['drug',
                                              'module',
                                              'cor',
                                              'p_value',
                                              'gene_sembol',
                                              'gene_symbol',
                                              'gene_coef'])
        return extend_coef_pd

    extend_coef_pd = extend_gene_coef_s(result_pd)
    return result_pd,extend_coef_pd


def merge_result(pathes):
    all_result_filter = []
    all_extend_coef_result = []
    for module_path in pathes:
        result_scale = pickle.load(open(module_path + 'linear_reuslt_cut_scale','rb'))
        result = pickle.load(open(module_path + 'linear_reuslt_cut', 'rb'))
        result_pd_scale,extend_coef_pd_scale = prepare_prediction_result(result_scale)
        result_pd_no_scale,extend_coef_pd_no_scale = prepare_prediction_result(result)
        all_result = pd.concat([result_pd_scale,result_pd_no_scale],
                               keys=['scale', 'not scale'])
        extend_coef_result = pd.concat([extend_coef_pd_scale , extend_coef_pd_no_scale],
                               keys=['scale', 'not scale'])
        all_result_filter.append(all_result)
        all_extend_coef_result.append(extend_coef_result)

    all_result_filter_pd = pd.concat(all_result_filter, keys=['no_fcg','with_fcg']).reset_index()
    all_result_filter_pd = all_result_filter_pd.rename(columns={'level_0':'fcg','level_1':'scale'})
    all_result_filter_pd = all_result_filter_pd.drop(columns='level_2')
    all_result_filter_pd.to_csv('results/all_predition_result_cut.csv')

    all_extend_coef_result_pd = pd.concat(all_extend_coef_result, keys=['no_fcg', 'with_fcg']).reset_index()
    all_extend_coef_result_pd = all_extend_coef_result_pd.rename(columns={'level_0': 'fcg', 'level_1': 'scale'})
    all_extend_coef_result_pd = all_extend_coef_result_pd.drop(columns='level_2')
    all_extend_coef_result_pd.to_csv('results/all_extend_coef_result_cut.csv')


def main():
    dataset_Y = read_csv('data/gene_processed.csv', index_col=0)
    dataset_X = read_csv('data/drug_processed.csv', index_col=0)
    dataset = pd.merge(dataset_X,dataset_Y, right_index=True, left_index=True,how='inner' )
    drugs  = ['Bortezomib',
              'Carfilzomib',
              'Ixazomib',
              'Oprozomib']

    path_1 = 'resultsv1/'
    path_2 = 'results/'

    model_prediction([path_1, path_2], drugs, dataset, scale=True)
    model_prediction([path_1, path_2], drugs, dataset, scale=False)
    merge_result([path_1, path_2])

    # linear_reuslt = pickle.load(open('results/linear_reuslt_scale'), 'rb')


