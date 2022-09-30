import benchmarkcorrelation.enrichr
from sklearn.metrics import roc_auc_score
import tqdm
import numpy as np

import benchmarkcorrelation.enrichr

def prediction(cormat):
    kegg = enrichr.get_library("KEGG_2021_Human")
    chea = enrichr.get_library("ChEA_2016")
    pheno = enrichr.get_library("MGI_Mammalian_Phenotype_Level_4_2021")
    humap = enrichr.get_library("huMAP")
    auc = predict(cormat, kegg)
    print('{}{}'.format("avg AUC (KEGG):".ljust(26), auc))
    auc = predict(cormat, chea)
    print('{}{}'.format("avg AUC (ChEA):".ljust(26), auc))
    auc = predict(cormat, pheno)
    print('{}{}'.format("avg AUC (Phenotype):".ljust(26), auc))
    auc = predict(cormat, humap)
    print('{}{}'.format("avg AUC (huMAP):".ljust(26), auc))

def predict(cormat, library):
    all_genes = []
    for k in list(library.keys()):
        library[k] = cormat.index.intersection(library[k])
        all_genes = all_genes + list(library[k])

    all_genes = sorted(list(set(all_genes)))

    scor = cormat.loc[all_genes, all_genes]
    cormatna = scor.replace(1, np.nan)

    aucs = []
    for k in tqdm.tqdm(library.keys()):
        geneset = list(library[k])
        aucs.append(roc_auc_score([x in library[k] for x in all_genes], cormatna.loc[:, geneset].mean(axis=1)))
    return np.mean(aucs)
