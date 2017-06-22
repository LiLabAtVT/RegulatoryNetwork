import numpy as np
import csv
import scipy.stats
import sklearn.preprocessing


def read_input(
        exp_mat_file, tf_names_file, priors_file):
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
    exp_mat_data = sklearn.preprocessing.scale(exp_mat_data, axis=1)
    if exp_mat_data.shape[1] != len(exp_mat_expname):
        raise ValueError('the column number should be same')
    if exp_mat_data.shape[0] != len(exp_mat_allgenename):
        raise ValueError('the row number should be same')
    # exp_mat_data is matrix B, gene * experiment
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

    # construct the matrix A, TF * experiment
    tf_mat = np.zeros((len(tf_nameuse), len(exp_mat_expname)))
    ind = 0
    for i in tf_nameuse:
        pos1 = exp_mat_allgenename.index(i)
        data1 = exp_mat_data[pos1, :]
        tf_mat[ind, :] = data1
        ind = ind + 1
    IN['tf_mat'] = tf_mat

    # read the prior data
    # construct the matrix Initial, gene * TF
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
        priors_data = np.zeros((len(exp_mat_allgenename), len(tf_nameuse)))
        for i in priors_matuse:
            tf1 = i[0]
            gene1 = i[1]
            rowpostf = tf_nameuse.index(tf1)
            rowposgene = exp_mat_allgenename.index(gene1)
            tfexpress1 = tf_mat[rowpostf, :]
            geneexpress1 = exp_mat_data[rowposgene, :]
            pearsonr1 = scipy.stats.pearsonr(tfexpress1, geneexpress1)
            priors_data[rowposgene, rowpostf] = pearsonr1[0]
        IN['priors.data'] = priors_data
    else:
        IN['priors.data'] = np.zeros(
            (len(exp_mat_allgenename), len(tf_nameuse)))
        print('No golden standard dataset used')
    return IN


def arrange_result(result_input):
    relation_total = list()
    relation_value = list()

    for i in range(0, len(result_input)):
        relationtf1 = result_input[i]
        for j in range(0, len(relationtf1)):
            relation1 = relationtf1[j]
            relation_total.append(relation1)
            relation_value.append(relation1[-1])
    relation_value = np.asarray(relation_value)
    relation_value = np.abs(relation_value)

    orderval = relation_value.argsort()[::-1]
    relation_valueorder = relation_value[orderval]
    relation_totalorder = list(relation_total[i] for i in orderval)
    # rankval = scipy.stats.rankdata(relation_valueorder)[::-1]
    rankval = range(1, len(relation_valueorder))
    return relation_totalorder, rankval
