# priors.pars is list of prior parameters; every entry is a vector of four
# elements: perc.tp, tp permutation number, perc.fp, fp permutation number
import numpy as np
import math
import warnings


def getPriors(
    exp_mat, tf_names, priors_mat, gs_mat, eval_on_subset,
        job_seed, perc_tp, perm_tp, perc_fp, perm_fp):

    # Note: We require the perc.*p and perm.*p vectors to be of
    # the same length as
    # in the following code only the same positions of the vectors
    # will be paired.
    a = [perc_tp, perm_tp, perc_fp, perm_fp]
    lengths = [len(x) for x in a]
    lengths = np.array(lengths)
    priors_par = np.zeros((1, 4))
    if len(np.unique(lengths)) != 1:
        raise ValueError(
            'Error parsing prior parameters: perc.tp, perc.fp, perm.fp don''t \
            have the same length')
    for pos in range(0, np.unique(lengths)[0]):
        rmndr = max(perm_tp[pos], perm_fp[pos]) % \
            min(perm_tp[pos], perm_fp[pos])
        if rmndr != 0:
            raise ValueError(
                'Error parsing prior parameters: Larger number of permutations is \
            not multiple of smaller number')
        d = np.arange(perm_fp[pos])
        a = np.array([perc_tp[pos]] * d.shape[0])
        b = np.tile(np.arange(perm_tp[pos]), d.shape[0])
        c = np.array([perc_fp[pos]] * d.shape[0])
        priors_pos = np.column_stack((a, b, c, d))
        priors_par = np.row_stack((priors_par, priors_pos))
    priors_par = np.delete(priors_par, 0, 0)

    priors = {}
    for i in range(0, len(priors_par)):
        pp = priors_par[i]
        name = 'frac_tp_' + str(pp[0]) + '_perm_' + str(pp[1]) + \
            '--frac_fp_' + str(pp[2]) + '_perm_' + str(pp[3]) + ' '
        priors[name] = np.zeros((exp_mat.shape[0], len(tf_names)))

        if pp[0] > 0 or pp[3] > 0:
            priors[name] = getPriorMatrix(
                priors_mat, pp, gs_mat, eval_on_subset, job_seed)
    return priors


def getPriorMatrix(priors, prior_pars, gs, from_subset, seed):
    # Given the gold standard, path of the input data, and a set of prior
    # parameters, this function returns a single -1, 0, 1
    # matrix of random priors

    perc_tp = prior_pars[0]
    perm_tp = prior_pars[1]
    perc_fp = prior_pars[2]
    perm_fp = prior_pars[3]

    if not from_subset:
        p_mat = makePriorMat(
            priors, perc_tp, perm_tp + seed, false_priors=False) \
            + makePriorMat(priors, perc_fp, perm_fp + seed, false_priors=True)

    return p_mat


def makePriorMat(priors, perc, perm, false_priors):
    # Creates prior matrix based on:
    #
    # priors - matrix with priors
    # perc - percentage of priors to use
    # perm - the permutation of this priors-perc combination
    # false.priors - whether to use false priors
    #
    # Notes: By setting the RNG seed using perm, we make sure that, for a given
    # priors-perm combination, matrices with higher perc values include all the
    # priors of those with lower perc values.
    # save the state of the RNG

    prior_order = np.where(priors != false_priors)
    num = len(prior_order[0])
    num_order = np.random.permutation(np.arange(num))
    prior_ordernew = []
    prior_ordernew.append(prior_order[0][num_order])
    prior_ordernew.append(prior_order[1][num_order])

    n_priors = math.floor(num * perc / 100)
    p_mat = np.zeros((priors.shape[0], priors.shape[1]))

    if n_priors > num:
        n_priors = num

    if n_priors > 0:
        for i in range(0, n_priors):
            p_mat[prior_ordernew[0][i]][prior_ordernew[1][i]] = \
                (-1)**false_priors

    # The above code ignores whether perc is set too high, but we should warn
    # the user.
    if (n_priors > len(prior_ordernew[0])):
        warnings.warn(
            'Percent of priors set too high. Only the max used')
    return p_mat


