This folder contains Python script to infer gene network by hard threshold, half threshold and soft threshold (LASSO) method from a expression matrix.

### Prerequist:

* The script is run at Python 3.6.0 Anaconda
* Python package: sys, os, numpy, csv, scipy.stats, sklearn.preprocessing, math

### Usage:
```bash
Pythonscript MainLasso.py [InputFileName] [ListOfTFNames] [GoldStandardFileName] [NumberOfEdges] [OutputFileName]
```
* [InputFileName] is the name of input file,the file format should follow data /GSE10670_ave_TopVar.csv. Each row is a gene, each column is a condition.

* [ListOfTFNames] is the name of all the transcription gene names. Default: tfs.csv

* [GoldStandardFileName] contains gold standard intercation pairs. Default: data/golddata.csv

* [NumberOfEdges] spcifies how many edges to write write to output file. Default: 100

* [OutputFileName] is the output file. Default: ./results/HalfThr.csv, ./results/HardThr.csv and ./results/SoftThr.csv

The [InputFileName] and [ListOfTFNames] are required for the input.

### Example command line usage:
```bash
Pythonscript MainLasso.py ./data/GSE10670_ave_TopVar.csv 
```
### Example output format can be found in ./results/HardThr.csv, ./results/HalfThr.csv and ./results/SoftThr.csv. Here is first few lines of this result.

```
regulator,target,score,rank
AT4G34590,AT2G16660,-0.9961,9
AT4G34590,AT5G54060,0.9949,10
AT4G34590,AT4G34710,0.9945,11
AT4G34590,AT1G53910,0.9941,12
AT2G46270,AT1G64110,0.9934,13
``` 
The output is a comma separated table with four columns.

|regulator|target|score|rank|
|---|---|---|---|
|AT4G34590|AT2G16660|-0.9961|9|
|AT4G34590|AT5G54060|0.9949|10|
|AT4G34590|AT4G34710|0.9945|11|
|AT4G34590|AT1G53910|0.9941|12|
|AT2G46270|AT1G64110|0.9934|13|


