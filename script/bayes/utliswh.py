import numpy as np
import csv
import scipy.stats


def read_input(
        exp_mat_file, tf_names_file,
        priors_file, gold_standard_file):
    # read the exp mat file, which is the data file contains the
    # gene expression
    IN = {}
    exp_mat = []
    with open(exp_mat_file, 'rt') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ')
        for row in spamreader:
            exp_mat.append(row)
    csvfile.close()
    # extract the name and value of the exp_mat
    exp_mat_expname = exp_mat[0]
    exp_mat_expname = exp_mat_expname[0]
    exp_mat_expname = exp_mat_expname.strip()
    exp_mat_expname = exp_mat_expname.split(',')
    # the first two should be a comma and 'GeneID', delete
    del exp_mat_expname[0:1]

    exp_mat_allgenename = []
    exp_mat_data = []
    for i in range(1, len(exp_mat)):
        exp_mati = exp_mat[i]
        exp_mati = exp_mati[0]
        exp_mati = exp_mati.strip()
        exp_mati = exp_mati.split(',')
        # the first one is a row number, delete it
        # del exp_mati[0]
        genename1 = exp_mati[0]
        del exp_mati[0]
        exp_mati = [float(ii) for ii in exp_mati]
        if sum(np.array(exp_mati) != 0) != 0:
            exp_mat_allgenename.insert(i - 1, genename1)
            exp_mat_data.append(exp_mati)
    exp_mat_data = np.array(exp_mat_data)
    if exp_mat_data.shape[1] != len(exp_mat_expname):
        raise ValueError('the column number should be same')
    if exp_mat_data.shape[0] != len(exp_mat_allgenename):
        raise ValueError('the row number should be same')
    IN['exp.mat.data'] = exp_mat_data
    IN['exp.mat.expname'] = exp_mat_expname
    IN['exp.mat.allgenename'] = exp_mat_allgenename

    # read the tf name file
    tf_names = []
    with open(tf_names_file, 'rt') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ')
        for row in spamreader:
            tf_names.append(row)
    csvfile.close()
    tf_nameuse = []
    for i in tf_names:
        if isinstance(i, list):
            i = i[0]
        if i in exp_mat_allgenename:
            tf_nameuse.append(i)
    # if there are same names in the tf, left only one
    tf_nameuse = list(set(tf_nameuse))
    IN['tf.names'] = tf_nameuse

    # read the prior data
    if priors_file is not None:
        priors_mat = []
        priors_matuse = []
        with open(priors_file, 'rt') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            for row in spamreader:
                priors_mat.append(row)
        csvfile.close()
        for i in priors_mat:
            if isinstance(i, list):
                i = i[0]
            data1 = i.split(',')
            if data1[0] in tf_nameuse and data1[1] in exp_mat_allgenename:
                priors_matuse.append(data1)
        priors_data = np.zeros((len(tf_nameuse), len(exp_mat_allgenename)))
        for i in priors_matuse:
            tf1 = i[0]
            gene1 = i[1]
            rowpos = tf_nameuse.index(tf1)
            colpos = exp_mat_allgenename.index(gene1)
            priors_data[rowpos, colpos] = 1
        IN['priors.data'] = priors_data
    else:
        IN['priors.data'] = None
        print('No priors data used')

    # read the golden standard data
    if gold_standard_file is not None:
        gold_mat = []
        gold_matuse = []
        with open(gold_standard_file, 'rt') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            for row in spamreader:
                gold_mat.append(row)
        csvfile.close()
        for i in gold_mat:
            if isinstance(i, list):
                i = i[0]
            data1 = i.split(',')
            if data1[0] in tf_nameuse and data1[1] in exp_mat_allgenename:
                gold_matuse.append(data1)
        gold_data = np.zeros((len(tf_nameuse), len(exp_mat_allgenename)))
        for i in gold_matuse:
            tf1 = i[0]
            gene1 = i[1]
            rowpos = tf_nameuse.index(tf1)
            colpos = exp_mat_allgenename.index(gene1)
            gold_data[rowpos, colpos] = 1
        IN['gold.standard.data'] = gold_data
    else:
        IN['gold.standard.data'] = None
        print('No golden standard dataset used')
    return IN

# given the condition names, create meta data data frame that assumes all
# observations are steady state measurements


def trivial_meta_data(names):
    meta_data = {}
    meta_data['condName'] = names
    meta_data['isTs'] = [False] * len(names)
    meta_data['is1stLast'] = ['e'] * len(names)
    meta_data['prevCol'] = ['NA'] * len(names)
    meta_data['del.t'] = ['NA'] * len(names)
    return meta_data

# given an expression matrix, create a trivial cluster stack - no (bi)clusters


def trivial_cluster_stack(datainp, colnameinp, rownameinp):
    clusterStack = []
    for i in range(0, len(rownameinp)):
        dic1 = {'cols': colnameinp}
        dic1['ncols'] = [len(colnameinp)] * len(colnameinp)
        dic1['rows'] = ['G1'] * len(colnameinp)
        dic1['nrows'] = [1] * len(colnameinp)
        dic1['resid'] = [float('nan')] * len(colnameinp)
        dic1['k'] = [1] * len(colnameinp)
        dic1['redExp'] = datainp[i, ]
        clusterStack.append(dic1)
    return clusterStack


def arrange_result(result_input):
    relation_total = list()
    relation_value = list()

    for i in range(0, len(result_input)):
        relationtf1 = result_input[i]
        for j in range(0, len(relationtf1)):
            relation1 = relationtf1[j]
            if relation1[0] != relation1[1]:
                relation_total.append(relation1)
                relation_value.append(relation1[-1])
    relation_value = np.asarray(relation_value)
    relation_value = np.abs(relation_value)

    orderval = relation_value.argsort()[::-1]
    relation_valueorder = relation_value[orderval]
    relation_totalorder = list(relation_total[i] for i in orderval)
    rankval = scipy.stats.rankdata(relation_valueorder, method='average')[::-1]

    return relation_totalorder, rankval
