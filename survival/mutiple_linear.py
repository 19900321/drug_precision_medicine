import pandas as pd
from sklearn import linear_model

def train_linear(data, gene_col, sur_col):
    regr = linear_model.LinearRegression()
    X = data[gene_col]
    Y = data[sur_col]
    regr.fit(X, Y)
    return regr


