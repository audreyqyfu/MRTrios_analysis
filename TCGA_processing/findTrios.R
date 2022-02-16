####################################
# function to find trios
# by matching methylation probes,
# gene expression, and copy number
# alteration (CNA)
####################################
# input
#
# methyl.data: methylation data matrix; probes in rows, individuals in columns
# expn.data: gene expression data matrix; genes in rows, individuals in columns
# cna.data: cna data matrix; genes in rows, individuals in columns
# entrez.list: master list of gene names with Entrez IDs
#
# output
# trios.rows: data matrix of 4 columns: gene name, row number
#             from each input dataset
function <- findTrios (methyl.data, expn.data, cna.data, entrez.list) {
    # identify a list of unique genes
    # for methylation probes
    
    
    # identify trios by gene name (Hugo_Symbol)
    # call function trios()
    
    
    # identify additional trios by Entrez ID
    # call function entrez()
    
    
    # find additional matches using R package org.Hs.eg.db

    
    # return a data matrix for matched trios
    # in each row: gene name, row numbers from each dataset

}
