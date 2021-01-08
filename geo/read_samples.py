import pandas as pd
import GEOparse

ges = GEOparse.get_GEO(filepath='data/geo/GSE55145_family.soft.gz')
ges_2 = GEOparse.get_GEO(filepath='data/geo/geo/GSE9782_family.soft.gz')
pivoted_control_samples = ges_2.pivot_samples('VALUE')

ges.phenotype_data.to_csv('data/geo/GSE55145_family_pheno.csv')
ges_2.phenotype_data.to_csv('data/geo/GSE9782_family_pheno.csv')