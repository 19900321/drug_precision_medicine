BiocManager::install("GEOquery")
library(Biobase)
library(GEOquery)

path = 'C://Users//yinyin//Desktop//Project//GDSC-drug-modules//geo//GSE55145_family.soft.gz'
GSE55145 <- getGEO(filename=path)
GSE55145_header = GSE55145@header
GSE55145_gsms = GSE55145@gsms
GSE55145_gpls = GSE55145@gpls

one_patient_example = GSE55145_gsms$GSM1336382@header$characteristics_ch1


path_2 = 'C://Users//yinyin//Desktop//Project//GDSC-drug-modules//geo//GSE9782_family.soft.gz'
GSE9782 <- getGEO(filename=path_2)
GSE9782@header
GSE9782_gsms = GSE9782@gsms
GSE9782_gpls = GSE9782@gpls

one_patient_example_GSE9782 = GSE9782_gsms$GSM246523@header$characteristics_ch1
