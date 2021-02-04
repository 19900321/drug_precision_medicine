import pandas as pd
import pickle
from scipy.stats import spearmanr, pearsonr
from matplotlib import pyplot
import matplotlib.pyplot as plt

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




def validate_commonpass(best_model_dict):
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

def spearman_cor(pateint_group_dict):

    level_dict = {'sCR': 6,
                 'CR': 5,
                 'VGPR': 4,
                 'SD': 3,
                 'PR': 2,
                 'PD': 1
                    }

    # use correlation
    for g, values in pateint_group_dict.items():
        pateint_data = values['pateint_data']
        pateint_data_selelct = pateint_data[(pateint_data['censpfs'] == 1) & pateint_data['censpfs'].notna()]

        expected = pateint_data_selelct['ttcpfs']
        predicted = pateint_data_selelct['predicted_response']
        cor, p_value = spearmanr(expected, predicted)
        pateint_group_dict[g].update({'pfs_cor': cor,
                                      'pfs_p_vlue': p_value})
        # plt.scatter(expected, predicted)
        # plt.show()


    # use level correlation
    for g, values in pateint_group_dict.items():
        pateint_data = values['pateint_data']
        pateint_data['reponse_level'] = pateint_data['bestrespsh'].apply(lambda x: level_dict.get(x))
        pateint_data_selelct_2 = pateint_data[pateint_data['reponse_level'].notna()]

        # use correlation
        expected = pateint_data_selelct_2['reponse_level']
        predicted = pateint_data_selelct_2['predicted_response']
        cor, p_value = spearmanr(expected, predicted)
        pateint_group_dict[g].update({'level_cor': cor,
                                      'level_p_vlue': p_value})
        # plt.scatter(expected, predicted)
        # plt.show()

    return pateint_group_dict


def main():
    path = 'results/module_ml_results/'
    prediction_df = pd.read_csv(path + 'all_predition_result_top_fea.csv')
    best_model_dict = select_model(prediction_df, path)
    pickle.dump(best_model_dict, open(path+'best_model_dict', 'wb'))

    pateint_group_dict = validate_commonpass(best_model_dict)
    pickle.dump(pateint_group_dict, open(path + 'validate/commonpass/pateint_group_log_dict_prediected', 'wb'))

    pateint_group_dict = pickle.load(open(path + 'validate/commonpass/pateint_group_log_dict_prediected', 'rb'))

    pateint_group_dict_eval = spearman_cor(pateint_group_dict)
    pickle.dump(pateint_group_dict_eval, open(path + 'validate/commonpass/pateint_group_dict_eval', 'wb'))
