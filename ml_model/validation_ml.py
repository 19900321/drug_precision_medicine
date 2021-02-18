import pandas as pd
import pickle
from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import pickle
import matplotlib.gridspec as gridspec


def select_model(prediction_df, path):

    prediction_df_group = prediction_df.groupby('drug').apply(lambda x: x.sort_values(by='cor',
                                                                                      ascending=False).iloc[0, :])
    prediction_df_group = prediction_df_group[['fcg', 'scale', 'module', 'drug']]
    prediction_df_group_dict = prediction_df_group.to_dict(orient='index')

    best_model_dict = {}

    for drug, values in prediction_df_group_dict.items():
        module_type_dict = {'all_module': 0, 'module_1': 1, 'module_2': 2}
        module_type = module_type_dict[prediction_df_group_dict[drug]['module']]

        scale_type_dict = {'not scale': '', 'scale': '_scale'}
        scale_type = scale_type_dict[prediction_df_group_dict[drug]['scale']]

        model_result = pickle.load(open(path + values['fcg'] + '/linear_reuslt_top_fea{}'.format(scale_type), 'rb'))

        model_result[drug][module_type].update(prediction_df_group_dict[drug])
        best_model_dict[drug] = model_result[drug][module_type]

    return best_model_dict


def predict_commonpass(best_model_dict):
    pateint_group_dict = pickle.load(open('results/commpass/subgroup/pateint_group_log_dict', 'rb'))

    simple_name = {'Bortezomib': 'bor',
                    'Carfilzomib': 'car',
                   }

    drugs = ['Bortezomib', 'Carfilzomib']

    pateint_group_dict_predicted = {}
    for drug in drugs:
        groups = [g for g in list(pateint_group_dict.keys()) if simple_name[drug] in g.lower()]

        for g in groups:
            genes_selected = set(best_model_dict[drug]['feature_selected']) & set(pateint_group_dict[g]['gene_data'].columns)
            gene_data = pateint_group_dict[g]['gene_data'][list(genes_selected)]

            model = best_model_dict[drug]['clf']
            predicted = model.predict(gene_data)
            pateint_data = pateint_group_dict[g]['patient_data']
            pateint_data['predicted_response'] = predicted
            pateint_group_dict_predicted[g] = {'gene_data': gene_data,
                                               'pateint_data': pateint_data}

    return pateint_group_dict_predicted



def evaluation_prediction(pateint_group_dict):

    # use correlation
    def survival_time(pateint_data):
        pateint_data_selelct = pateint_data[(pateint_data['censpfs'] == 1) & pateint_data['censpfs'].notna()]

        expected = pateint_data_selelct['ttcpfs']
        predicted = pateint_data_selelct['predicted_response']
        cor, p_value = spearmanr(expected, predicted)
        # plt.scatter(expected, predicted)
        # plt.show()
        return cor, p_value

    def groups_cox(pateint_data, ttc, cens, cut):
        n = pateint_data.shape[0]//cut
        pateint_data = pateint_data.sort_values(by='predicted_response')
        T_h = pateint_data[ttc][-n:]
        E_h = pateint_data[cens][-n:]
        T_l = pateint_data[ttc][0:n]
        E_l = pateint_data[cens][0:n]
        results = logrank_test(T_h,
                               T_l,
                               event_observed_A=E_h,
                               event_observed_B=E_l)


        def plot_figure(T_h,  E_h, T_l, E_l):
            kmf_h = KaplanMeierFitter()
            kmf_l = KaplanMeierFitter()

            kmf_h.fit(durations=T_h,
                      event_observed=E_h,
                      label="higher expression (n={})".format(len(T_h)))
            kmf_l.fit(durations=T_l,
                      event_observed=E_l,
                      label="lower expression (n={})".format(len(T_l)))

            fig = plt.figure(figsize=(8, 8))
            gs = gridspec.GridSpec(4, 1)
            ax1 = fig.add_subplot(gs[0:2, :])
            kmf_h.plot_survival_function(ax=ax1)
            kmf_l.plot_survival_function(ax=ax1)
            ax1.set_xlabel("Times of {}".format(ttc))
            ax1.set_ylabel(" {} (%)".format(cens))
            ax1.set_title("KMF of {} from cut off {}".format(g,  n))
            plt.show()

        return results

    def repnose_level(pateint_data):
        level_dict = {'sCR': 6,
                      'CR': 5,
                      'VGPR': 4,
                      'SD': 3,
                      'PR': 2,
                      'PD': 1
                      }
        pateint_data['reponse_level'] = pateint_data['bestrespsh'].apply(lambda x: level_dict.get(x))
        pateint_data_selelct_2 = pateint_data[pateint_data['reponse_level'].notna()]

        # use correlation
        expected = pateint_data_selelct_2['reponse_level']
        predicted = pateint_data_selelct_2['predicted_response']
        cor, p_value = spearmanr(expected, predicted)
        # plt.scatter(expected, predicted)
        # plt.show()
        return cor, p_value

    for g, values in pateint_group_dict.items():
        pateint_data = values['pateint_data']

        # by split top 10, bottom 10 group
        ttc = 'ttcpfs'
        cens= 'censpfs'
        cut = 2
        results = groups_cox(pateint_data, ttc, cens, cut)
        pateint_group_dict[g].update({'group_test': results._test_statistic[0],
                                      'group_p_vlue': results.p_value})

        # by reponse level
        cor, p_value = repnose_level(pateint_data)
        pateint_group_dict[g].update({'level_cor': cor,
                                      'level_p_vlue': p_value})

        # by survival time of die patients
        cor_s, p_value_s = survival_time(pateint_data)
        pateint_group_dict[g].update({'pfs_cor': cor_s,
                                      'pfs_p_vlue': p_value_s})


    return pateint_group_dict


def main():
    path = 'results/module_ml_results/'
    # prediction_df = pd.read_csv(path + 'all_predition_result_top_fea.csv')
    # best_model_dict = select_model(prediction_df, path)
    # pickle.dump(best_model_dict, open(path+'best_model_dict', 'wb'))
    #
    # pateint_group_dict = predict_commonpass(best_model_dict)
    # pickle.dump(pateint_group_dict, open(path + 'validate/commonpass/pateint_group_log_dict_prediected', 'wb'))

    pateint_group_dict = pickle.load(open(path + 'validate/commonpass/pateint_group_log_dict_prediected', 'rb'))

    pateint_group_dict_eval = evaluation_prediction(pateint_group_dict)
    pickle.dump(pateint_group_dict_eval, open(path + 'validate/commonpass/pateint_group_dict_eval', 'wb'))

