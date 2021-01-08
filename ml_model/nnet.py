import numpy as np
import pandas as pd
import pickle


from pandas import read_csv
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from ml_model.ml import get_final_columns_dict


def nerual_model(dataset, drug,final_columns_dict):
    dataset = dataset[dataset[drug].notnull()]
    X = dataset[final_columns_dict[drug]]
    length_varibale = X.shape[1]
    Y = dataset[[drug]]
    # create model
    def wider_model():
        # create model
        model = Sequential()
        model.add(Dense(length_varibale, input_dim=length_varibale, kernel_initializer='normal', activation='relu'))
        model.add(Dense(50, kernel_initializer='normal'))
        model.add(Dense(1, kernel_initializer='normal'))
        # Compile model
        model.compile(loss='mean_squared_error', optimizer='adam')
        return model

    # evaluate model with standardized dataset
    estimators = []
    estimators.append(('standardize', StandardScaler()))
    estimators.append(('mlp', KerasRegressor(build_fn=wider_model, epochs=100, batch_size=5, verbose=0)))
    pipeline = Pipeline(estimators)
    kfold = KFold(n_splits=5)
    results = cross_val_score(pipeline, X, Y, cv=kfold)
    print("Wider: %.2f (%.2f) MSE" % (results.mean(), results.std()))
    return estimators, results


data_saved = pickle.load(open('../results/result_2.out', 'rb'))
dataset_Y = read_csv('../data/gene_processed.csv', index_col=0)
dataset_X = read_csv('../data/drug_processed.csv', index_col=0)
dataset = pd.merge(dataset_X,dataset_Y, right_index=True, left_index=True,how='inner' )
drugs  = ['Bortezomib',
          'Carfilzomib',
          'Ixazomib',
          'Oprozomib']

final_columns_dict = get_final_columns_dict(drugs, data_saved)
drug = 'Bortezomib'
nerual_model(dataset, drug,final_columns_dict)