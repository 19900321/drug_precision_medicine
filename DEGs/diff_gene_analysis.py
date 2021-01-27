from bioinfokit import visuz
import pandas as pd
import pickle

df = pd.read_csv('../results/deg/edgeR/carf_all_line_all.csv')
df = df.dropna(subset=['symbol'])
drug_name = 'Carfilzomib'

df = pd.read_csv('results/mm/correlation/all_drugs/extend_all_cor_pvalue_fdr.csv')
df = df[df['drug']=='Carfilzomib']
df = df.dropna(subset=['symbol'])
gene_s = ['ENSG00000116497',
                'ENSG00000143546',
                'ENSG00000160678',
                'ENSG00000163221',
                'ENSG00000197747',
                'ENSG00000163191',
                'ENSG00000197956',
                'ENSG00000163220',
                'ENSG00000196154',
                'ENSG00000189171']

anno_gene_dict = pickle.load(open('results/commpass/ensembol_symbol_dict_commpass', 'rb'))

anno_gene_dict = {g: anno_gene_dict.get(g) for g in gene_s}


visuz.gene_exp.volcano(df=df,
                       lfc="log2FoldChange",
                       pv="padj",
                       geneid="Unnamed: 0",
                       genenames=anno_gene_dict,
                       gstyle=1,
                       color=("#00239CFF", "grey", "#E10600FF"),
                       sign_line=True,
                       lfc_thr=0.5,
                       pv_thr=0.05,
                       plotlegend=True,
                       legendpos='upper right',
                       figtype='svg',
                       dotsize=40,
                       valpha=0.5,
                       axylabel='-log10(FDR)',
                       figname='The DEGs for {}'.format(drug_name),
                       axtickfontsize=12,
                       axlabelfontsize=12,
                       dim=(8,6),
                       show=True
                       )

#
visuz.gene_exp.volcano(df=df,
                       lfc="log2FoldChange",
                       pv="padj",
                       geneid="symbol",
                       genenames='deg',
                       gstyle=1,
                       color=("#00239CFF", "grey", "#E10600FF"),
                       sign_line=True,
                       lfc_thr=2,
                       pv_thr=0.05,
                       plotlegend=True,
                       legendpos='upper right',
                       figtype='svg',
                       dotsize=40,
                       valpha=0.5,
                       axylabel='-log10(FDR)',
                       figname='The DEGs for {}'.format(drug_name),
                       axtickfontsize=12,
                       axlabelfontsize=12,
                       dim=(8,6),
                       show=True
                       )


visuz.gene_exp.volcano(df=df,
                       lfc="cor",
                       pv="fdr",
                       geneid="symbol",
                       genenames='deg',
                       gstyle=1,
                       color=("#00239CFF", "grey", "#E10600FF"),
                       sign_line=True,
                       lfc_thr=0.3,
                       pv_thr=0.01,
                       plotlegend=True,
                       legendpos='upper right',
                       figtype='svg',
                       dotsize=40,
                       valpha=0.5,
                       xlm=(-1, 1, 0.5),
                       ylm=(0, 5, 1),
                       axxlabel='Drug response',
                       axylabel='-log10(FDR)',
                       figname='The DEGs for {}'.format(drug_name),
                       axtickfontsize=12,
                       axlabelfontsize=12,
                       dim=(10,5),
                       show=True
                       )