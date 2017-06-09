This folder contains scripts for Bayesian regression.
The data used:
'golddata.csv', which is golden standard dataset from Dapseq, first column is the TFs and
the second column is the target genes.
'GSE10670_ave_TopVar.csv': which is the expression matrix, each row is a gene (including the expression
of genes and TFs) and each row is an experiment.
'tfs.csv': which is a list of the names of all the TFs.

The result will be save in a folder named 'result1'. Three folders are in 'result1',
'halfThr', 'hardThr' and 'softThr', each means a constrain method for regression.
In each of the three folders, every single recognized TFs according to the
'tfs.csv' file are save as a '.csv' file by their name. There are three columns in
each the TFs name '.csv' files. The first column is the TF names, they are all same
for each file. The second column is the inferred target genes and the third column
is the inferred weight. A positive weight means the TF might increase the expression
of the gene and a negative weight means the TF might inhibit the express of the
gene.

