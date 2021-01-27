import pandas as pd

def load_genes_patients():
    file_path = 'resources/MMRF_CoMMpass_IA15a_E74GTF_Cufflinks_Gene_FPKM.txt'
    data_pd = pd.read_csv(file_path, sep='\t', index_col=0)
    data_pd = data_pd.drop(['Location'], axis=1)
    data_pd.shape
    return data_pd

def load_survive_patients():
    file_path = 'resources/OS_893_MM_RNAseqsamples.txt'
    data_pd = pd.read_csv(file_path, sep='\t', index_col=0)
    data_pd.shape
    return data_pd


