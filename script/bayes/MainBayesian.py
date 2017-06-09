# code transfered from paper "Inferring gene regulatory networks by integrating
# ChIP-seq/chip and  transcriptome data via Bayesian method"
# ref '1. Greenfield, A., Hafemeister, C. & Bonneau, R. Robust data-driven
# incorporation of prior knowledge into the inference of dynamic
# regulatory networks.
# Bioinformatics 29, 1060–1067 (2013).'


import numpy as np  # include the numpy pacage
import os  # include operating system
import random
import warnings
import csv
import networkx as nx
import matplotlib.pyplot as plt
from utliswh import read_input
from priorswh import getPriors
from mi_and_clrwh import mi
from mi_and_clrwh import mixedCLR
from bayesianRegressionwh import BBSR

# change the path to the folder containing the code and the data,
os.chdir("/Users/Craftsman/Documents/Wei He/Blacksburg/Study/2017\
 spring/course/Problem Solving GBCB 5874/project/RegulatoryNetwork/script/\
bayes")


''' main function '''
''' ##############################################'''
# if there is no input for the dictionary, put the value None
# warnings.filterwarnings('error')
currentpath = os.getcwd()
PAR = {'input.dir': '/Users/Craftsman/Documents/Wei He/Blacksburg/Study/2017\
 spring/course/Problem Solving GBCB 5874/project/RegulatoryNetwork/script/\
bayes/',  # change the path to the folder containing the code and the data,
                    'exp.mat.file':
                    'GSE10670_ave_TopVar.csv',
                    # change for different expression files
                    'tf.names.file': 'tfs.csv',
                    'priors.file': 'golddata.csv',
                    'gold.standard.file': 'golddata.csv',
                    'job.seed': 42,
                    'save.to.dir': currentpath,
                    'max.preds': 12,
                    'mi.bins': 10,
                    'cores': 8,
                    'delT.max': 110,
                    'delT.min': 0,
                    'tau': 45,
                    'perc.tp': np.array([100]),
                    'perm.tp': np.array([1]),
                    'perc.fp': np.array([0]),
                    'perm.fp': np.array([1]),
                    'eval.on.subset': False,
                    'method': 'BBSR',
                    'prior.weight': 2.8}
data = read_input(
    PAR['input.dir'], PAR['exp.mat.file'],
    PAR['tf.names.file'], PAR['priors.file'], PAR['gold.standard.file'])
# order genes so that TFs come before the other genes
tf_mat = np.zeros((len(data['tf.names']), len(data['exp.mat.expname'])))
priortf_mat = np.zeros((len(data['tf.names']), len(data['tf.names'])))
ind = 0
rowindTF = []
for i in data['tf.names']:
    pos1 = data['exp.mat.allgenename'].index(i)
    rowindTF.append(pos1)
    data1 = data['exp.mat.data'][pos1, :]
    if data['priors.data'] is not None:
        dataprior1 = data['priors.data'][:, pos1]
        priortf_mat[:, ind] = dataprior1
    tf_mat[ind, :] = data1
    ind = ind + 1

Maxgene = 2000
if len(data['tf.names']) > Maxgene / 2:
    raise ValueError('decrease the TF numbers')

# if gene number bigger than 2000, cut it to 2000
if len(data['exp.mat.allgenename']) > Maxgene:
    varval = []
    for i in range(0, data['exp.mat.data'].shape[0]):
        data1 = data['exp.mat.data'][i, :]
        varval.append(np.var(data1))
    varval = np.array(varval)
    var_order = np.argsort(varval)[::-1]

    expression_data = np.zeros((Maxgene, data['exp.mat.data'].shape[1]))

    expression_data[0: tf_mat.shape[0], :] = tf_mat
    genename = []
    for i in range(0, Maxgene - tf_mat.shape[0]):
        ind1 = var_order[i]
        expression_data[tf_mat.shape[0] + i, :] = data['exp.mat.data'][ind1, :]
        genename1 = data['exp.mat.allgenename'][ind1]
        genename.append(genename1)

    if data['priors.data'] is not None:
        prior_data = np.zeros((data['priors.data'].shape[0], Maxgene))
        prior_data[:, 0: tf_mat.shape[0]] = priortf_mat
        for i in range(0, Maxgene - tf_mat.shape[0]):
            ind1 = var_order[i]
            prior_data[:, tf_mat.shape[0] + i] = data['priors.data'][:, ind1]
        data['priors.data'] = prior_data
        data['gold.standard.data'] = prior_data
    genename = data['tf.names'] + genename
    data['exp.mat.data'] = expression_data
    data['exp.mat.allgenename'] = genename


# set the random seed
if PAR['job.seed'] is not None:
    random.seed(PAR['job.seed'])
    print('RNG seed has been set to ' + str(PAR['job.seed']) + '\n')
else:
    ignore = random.random()
SEED = random.getstate()

print('Output dir:' + PAR['save.to.dir'] + '\n')
if not os.path.exists(PAR['save.to.dir']):
    os.makedirs(PAR['save.to.dir'])  # create folder

# .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# create design and response matrix
# .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

X = tf_mat
Y = data['exp.mat.data']


#  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# parse priors parameters and set up priors list
#  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
priors = getPriors(
    data['exp.mat.data'], data['tf.names'], data['priors.data'],
    data['gold.standard.data'], PAR['eval.on.subset'], PAR['job.seed'],
    PAR['perc.tp'], PAR['perm.tp'], PAR['perc.fp'], PAR['perm.fp'])

#  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# main loop
#  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# for prior_name in priors.keys():
prior_name = list(priors.keys())[0]
text = 'Method: ' + PAR['method'] + '\nWeight: ' + \
    str(PAR['prior.weight'])\
    + '\nPriors: ' + prior_name + '\n'
print(text)
prior = priors[prior_name]
# set the prior weights matrix
no_pr_weight = 1
if np.sum(prior != 0) > 0:
    if PAR['prior.weight'] == no_pr_weight:
        warnings.warn(
            'Priors present, but they will not be used, because' +
            str(PAR['prior.weight']) + 'is set to' + str(no_pr_weight))
    if PAR['method'] == 'BBSR':
        no_pr_weight = 1 / PAR['prior.weight']
weight_mat = np.full(
    (len(data['tf.names']), data['exp.mat.data'].shape[0]),
    no_pr_weight)
weight_mat[prior != 0] = PAR['prior.weight']
betas = []
betas_resc = []

# fill mutual information matrices
print('Calculating MI\n')
Ms = mi(
    Y.transpose(), Y.transpose(), nbins=PAR['mi.bins'],
    cpu_n=PAR['cores'])
np.fill_diagonal(Ms, 0)


# get CLR matrix
print('Calculating CLR Matrix\n')
clr_mat = mixedCLR(Ms, Ms)
clr_matuse = np.zeros((clr_mat.shape[0], len(data['tf.names'])))
ind = 0
for i in data['tf.names']:
    pos1 = data['exp.mat.allgenename'].index(i)
    data1 = clr_mat[:, pos1]
    clr_matuse[:, ind] = data1
    ind = ind + 1
clr_mat = clr_matuse
# get the sparse ODE models
# get name
print('Calculation sparse ODE models\n')
x = BBSR(
    X, Y, clr_mat,
    PAR['max.preds'], no_pr_weight, weight_mat,
    PAR['cores'])
print('\n')

# our output will be a list holding two matrices: betas and betas.resc
bs_betas = np.zeros((Y.shape[0], X.shape[0]))

for i in range(0, len(x)):
    x1 = x[i]
    bs_betas[x1['ind'], x1['pp']] = x1['betas']

relationsall = []
for i in range(0, len(data['tf.names'])):
    tf1 = []
    tfname = data['tf.names'][i]
    for j in range(0, bs_betas.shape[0]):
        beta1 = bs_betas[j]
        gene1 = data['exp.mat.allgenename'][j]
        if beta1[i] != 0:
            relation1 = [tfname, gene1, beta1[i]]
            tf1.append(relation1)
    relationsall.append(tf1)

tfrelationsnum = np.zeros(len(relationsall))
ind = 0
for i in relationsall:
    len1 = len(i)
    tfrelationsnum[ind] = len1
    ind = ind + 1

posmax = np.where(tfrelationsnum == tfrelationsnum.max())[0]

# save the result
pathsave = PAR['save.to.dir'] + '/result1'

if not os.path.exists(pathsave):
    os.makedirs(pathsave)  # create folder

os.chdir(pathsave)

for i in range(0, len(data['tf.names'])):
    if os.path.isfile(data['tf.names'][i] + '.csv'):
        os.remove(data['tf.names'][i] + '.csv')
    with open(data['tf.names'][i] + '.csv', 'w') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
        for i in relationsall[i]:
            ii = [i[0], '    ' + i[1], '    ' + str(i[2])]
            wr.writerow(ii)
    myfile.close()


# draw one example
tf1rel = relationsall[0]
max_num = 10
if len(tf1rel) > max_num:
    tf1rel = tf1rel[0: max_num]
G = nx.Graph()
tfname = tf1rel[0][0]
G.add_node(tfname)
edge_color = []
labels = {}
labels[0] = tfname
ind = 1
for i in tf1rel:
    gene1 = i[1]
    weight = round(float(i[2]), 2)
    if weight > 0:
        colr = 'red'
    else:
        colr = 'blue'
    edge_color.append(colr)
    labels[ind] = gene1
    G.add_edge(tfname, gene1, weight=round(i[2], 2), color=colr)
    ind = ind + 1
edge_labels = dict([((u, v, ), d['weight']) for u, v, d in G.edges(data=True)])
pos = nx.spring_layout(G)
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=14)
nx.draw(G, pos, node_color='white', edge_color=edge_color,
        width=2, with_labels=True, font_size=14)
plt.show()

