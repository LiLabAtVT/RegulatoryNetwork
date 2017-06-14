This folder contain R script to calculate partial correlation from a expression matrix.

### Prerequist:
* R package: GeneNet
* R package: rshape


### Usage:
```bash
Rscript PCOR.R [InputFileName] [GoldStandardFileName] [NumberOfEdges] [OutputFileName]
```

* [InputFileName] is the name of input file, the file format should follow data/GSE10670_ave_TopVar.csv. Each row is a gene, each column is a condition.

* [GoldStandardFileName] contains gold standard interaction pairs. This is not used in partial correlation calculation. Default: data/golddata.csv

* [NumberOfEdges] specifies how many edges to write to output file. Default: 100

* [OutputFileName] is the output file. Default: ./results/pcor.csv

#### Example command line usage:
```bash
Rscript script/pcor/PCOR.R ./data/GSE10670_ave_TopVar.csv
```
### Example output format can be found in ./results/pcor.csv. Here is first few lines of this result.

```
"regulator","target","score","rank"
"AT2G35730","AT3G60140",0.00165,1.5
"AT3G60140","AT2G35730",0.00165,1.5
"AT1G01480","AT3G60140",0.00148,3.5
"AT3G60140","AT1G01480",0.00148,3.5
"AT2G35730","AT1G01480",0.00130,5.5
```

The output is a comma separated table with four columns.

|"regulator"|"target"|"score"|"rank"|
|---|---|---|---|
|"AT2G35730"|"AT3G60140"|0.00165|1.5|
|"AT3G60140"|"AT2G35730"|0.00165|1.5|
|"AT1G01480"|"AT3G60140"|0.00148|3.5|
|"AT3G60140"|"AT1G01480"|0.00148|3.5|
|"AT2G35730"|"AT1G01480"|0.00130|5.5|
