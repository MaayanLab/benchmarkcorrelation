# Benchmark Correlation

A suit of benchmark tools to check the correctness of a gene-gene correlation matrix. The benchmark includes the following tests:

- test gene identifier
    - are all gene identifieres known
    - are they homogeneous (e.g. only ensembl IDs) ior mixed (e.g. ensemble ID and gene symbol)
    - detect organism (human, mouse)
    - detect duplicates

- test symmetry
    - test dimensions
    - test if matrix is symmetric

- test known correlations
    - test whether known correlations are preserved

- test prediction power
    - calculate AUCs for gene-set predictions


## Support modules
 - Load gene sets from Enrichr
 - load biomart gene information / identifier mappings
 - Gene-function prediction

