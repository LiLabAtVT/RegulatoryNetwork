This folder contains Python scripts to infer gene network by Bayesian regression from a expression matrix.

### Prerequist:

* The script is run at Python 3.6.0 Anaconda
* Python package: sys, os, warnings, numpy, random, csv, scipy.stats, sklearn.preprocessing, scipy.special, itertools, sklearn.metrics, math.

### Usage:
```bash
Pythonscript MainBayesian.py [InputFileName] [ListOfTFNames] [GoldStandardFileName] [NumberOfEdges] [OutputFileName]
```

* [InputFileName] is the name of input file,the file format should follow data /GSE10670_ave_TopVar.csv. Each row is a gene, each column is a condition.

* [ListOfTFNames] is the name of all the transcription gene names. Default: tfs.csv

* [GoldStandardFileName] contains gold standard intercation pairs. Default: data/goldata.csv

* [NumberOfEdges] spcifies how many edges to write write to output file. Default: 100

* [OutputFileName] is the output file. Default: ./results/Bayes.csv

The [InputFileName] and [ListOfTFNames] are required for the.

### Example command line usage:
```bash
Pythonscript MainBayesian.py ./data/GSE10670_ave_TopVar.csv 
```
### Example output format can be found in ./results/Bayes.csv. Here is first few lines of this result.

```
regulator,target,score,rank
AT3G51960,AT3G51960,1,1
AT5G42910,AT5G42910,1,2
AT3G51960,AT2G41230,0.9933,3
AT2G46270,AT5G53710,0.9928,4
AT3G51960,AT2G45570,0.9915,5
``` 
The output is a comma separated table with four columns.

|regulator|target|score|rank|
|---|---|---|---|
|AT3G51960|AT3G51960|1|1|
|AT5G42910|AT5G42910|1|2|
|AT3G51960|AT2G41230|0.9933|3|
|AT2G46270|AT1G01480|0.9928|4|
|AT3G51960|AT1G01480|0.9915|5|

