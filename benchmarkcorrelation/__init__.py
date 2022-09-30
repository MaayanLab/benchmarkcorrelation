from functools import update_wrapper
import pandas as pd
import numpy as np
import random
import os
import sys
import tqdm
import feather
import csv

from matplotlib import pyplot as plt

import biomart

import benchmarkcorrelation.enrichr
import benchmarkcorrelation.prediction
import benchmarkcorrelation.geneids

from importlib import reload
reload(benchmarkcorrelation.enrichr)
reload(benchmarkcorrelation.prediction)
reload(benchmarkcorrelation.geneids)


def load_correlation_file(file: str):
    mat = 0
    delimiter = ","
    try:
        if file.endswith(".f"):
            print('{}{}'.format("Matrix format:".ljust(26), "feather"))
            mat = pd.read_feather(file)
            if "index" in mat.columns:
                mat.set_index("index", inplace=True)
            else:
                mat.index = mat.columns
        elif file.endswith((".csv", ".tsv", ".txt")):
            print('{}{}'.format("Matrix format:".ljust(26), "text"))
            with open(file, 'r') as csvfile:
                dialect = csv.Sniffer().sniff(csvfile.readline())
                delimiter = dialect.delimiter
            mat = pd.read_csv(file, sep=delimiter)
    except Exception:
        sys.exit("Matrix could not be read")
    return mat

def bench_stats(cormat):
    file_size = cormat.memory_usage(deep=True).sum()/1000000
    file_dimensions = cormat.shape
    square = file_dimensions[0] == file_dimensions[1]
    same_id_order = True
    symmetric = True
    if square:
        for i in range(file_dimensions[0]):
            if cormat.index[i] != cormat.columns[i]:
                same_id_order = False
                break
        random_pos_1 = random.sample(range(file_dimensions[0]), min(file_dimensions[0], 100))
        random_pos_2 = random.sample(range(file_dimensions[0]), min(file_dimensions[0], 100))
        for i in range(len(random_pos_1)):
            if cormat.iloc[random_pos_1[i], random_pos_2[i]] != cormat.iloc[random_pos_2[i], random_pos_1[i]]:
                symmetric = False
                break
    else:
        symmetric = False
        same_id_order = False

    print('{}{}{}'.format("File size:".ljust(26), int(file_size), "MB"))

    print('{}{}'.format("Shape:".ljust(26), file_dimensions))

    print('{}{}'.format("Square:".ljust(26), square))

    print('{}{}'.format("ID order identical:".ljust(26), same_id_order))

    print('{}{}'.format("Matrix symmetric:".ljust(26), symmetric))

    print('{}{}'.format("Contains NA:".ljust(26), cormat.isnull().values.any()))

    #temp = np.array(cormat, dtype = np.float32)
    #np.fill_diagonal(temp, None)
    print('{}{}'.format("Correlation mean:".ljust(26), np.nanmean(cormat.iloc[0:min(cormat.shape[0],10000), 0:min(cormat.shape[0],10000)].astype("float32"))))
    print('{}{}'.format("Correlation STD:".ljust(26), np.nanstd(cormat.iloc[0:min(cormat.shape[0],10000), 0:min(cormat.shape[0],10000)].astype("float32"))))
    np.nanmean(cormat.iloc[0:10000, 0:10000].astype("float32"))

def gene_ids(cormat):
    upper = 0
    lower = 0
    species = ""
    for i in cormat.index:
        if i.isupper():
            upper = upper + 1
        else:
            lower = lower + 1
    if upper > lower:
        species = "human"
    else:
        species = "mouse"
    print('{}{}'.format("Predicted species:".ljust(26), species))
    gene_id_map = geneids.get_ensembl_mappings(species)
    identifiers = set.union(set(cormat.index), set(cormat.columns))
    symbol_overlap = identifiers.intersection(set(gene_id_map.iloc[:,0]))
    entrezid_overlap = identifiers.intersection(set(gene_id_map.iloc[:,1]))
    print('{}{}'.format("Total unique identifiers:".ljust(26), len(identifiers)))
    print('{}{}'.format("Gene Symbol overlap:".ljust(26), len(set(symbol_overlap))/len(identifiers)))
    print('{}{}'.format("Ensembl ID overlap:".ljust(26), len(set(entrezid_overlap))/len(identifiers)))

def get_data_path() -> str:
    path = os.path.join(
        os.path.dirname(__file__),
        'data/'
    )
    return(path)

def compare_known_cor(cormat):
    old_cor = pd.read_feather(get_data_path()+"known_cor.f")
    old_cor.index = old_cor.columns
    inter_row = cormat.index.intersection(old_cor.index)
    inter_column = cormat.index.intersection(old_cor.columns)
    flat_new = cormat.loc[inter_row, inter_column].replace(1, 1).to_numpy().flatten()
    flat_old = old_cor.loc[inter_row, inter_column].replace(1, 1).to_numpy().flatten()
    cor_match = np.corrcoef(flat_new, flat_old)[1,1]
    print('{}{}'.format("Correlation similarity:".ljust(26), cor_match))

def benchmark(file: str, identifiers=True, correlation=True, pred=True, format=True):
    # run benchmark
    print("------- Load file -------")
    cormat = load_correlation_file(file)
    # 0. Basic stats (file name / file size)

    if format:
        print("------- Matrix format -------")
        bench_stats(cormat)

    if identifiers:
        print("------- Gene Identifiers -------")
        gene_ids(cormat)

    if correlation:
        print("------- Known correlation -------")
        compare_known_cor(cormat)
    
    if pred:
        print("------- Gene function prediction -------")
        prediction.prediction(cormat)

