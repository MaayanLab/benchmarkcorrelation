import pandas as pd
import numpy as np
import biomart

def get_ensembl_mappings(species):
    # Set up connection to server
    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')
    if species == "mouse":
        mart = server.datasets['mmusculus_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'mgi_symbol', 'ensembl_gene_id', 'gene_biotype']
    else:
        mart = server.datasets['hsapiens_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'hgnc_symbol', 'ensembl_gene_id', 'gene_biotype']                                                     
    # Get the mapping between the attributes                                    
    response = mart.search({'attributes': attributes})
    data = response.raw.data.decode('ascii')
    ensembl_ids = []
    # Store the data in a dict                                                  
    for line in data.splitlines():                                              
        line = line.split('\t')                                
        ensembl_ids.append(line)
    gene_map = pd.DataFrame(ensembl_ids)
    gene_map.index = gene_map.iloc[:,0]
    nn = np.where(gene_map.iloc[:,1] == "")[0]
    gene_map.iloc[nn, 1] = gene_map.iloc[nn, 2]
    return gene_map