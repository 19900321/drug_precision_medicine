import mygene
mg = mygene.MyGeneInfo()
from collections import defaultdict
import pandas as pd


def engs_gene_name(esemb_ids):
    sym_list, sym_dict = [], defaultdict()
    for i in mg.querymany(esemb_ids, scopes='ensembl.gene', fields='symbol'):
        try:
            sym_dict[i['query']] = i['symbol']
            sym_list.append(sym_dict)
        except:
            continue
    return sym_list, sym_dict
