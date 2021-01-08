import pandas as pd
import numpy as np


def map_pateints_id_compass():
    common_patient = pd.read_csv('data/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP.csv', sep=';')
    dataset_gene_count = pd.read_csv('data/MMRF_CoMMpass_IA15a_E74GTF_HtSeq_Gene_Counts.txt',
                                     sep='\t',
                                     index_col=0)

    gene_pateints = list(dataset_gene_count.columns)
    gene_patient_id = []
    for i,rows in common_patient.iterrows():
        test_gene_patient_id = rows.public_id + '_' + str(rows.line) + '_BM'
        if test_gene_patient_id in gene_pateints:
            gene_patient_id.append(test_gene_patient_id)
        else:
            test_gene_patient_id_alter = rows.public_id + '_' + str(rows.line) + '_PB'
            if test_gene_patient_id_alter in gene_pateints:
                gene_patient_id.append(test_gene_patient_id_alter)
            else:
                gene_patient_id.append(None)

    common_patient['gene_patient_id'] = gene_patient_id

    return common_patient


common_patient = map_pateints_id_compass()
common_patient.to_csv('data/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP_map.csv', index=None)

