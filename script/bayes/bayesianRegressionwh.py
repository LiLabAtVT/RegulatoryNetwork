import sklearn.preprocessing
import scipy.special
import numpy as np
import sys
import itertools


def BBSR(X, Y, clr_mat, nS, no_pr_val, weight_mat, cores):
    # Scale and permute design response matrix
    X = sklearn.preprocessing.scale(X, axis=1)
    Y = sklearn.preprocessing.scale(Y, axis=1)
    if nS > X.shape[0]:
        nS = X.shape[0] - 1
    # number of genes
    G = Y.shape[0]
    # max number of possible predictor (number of TFs)
    K = X.shape[0]
    # predictors that will be used in the regression
    pp = np.zeros((G, K), dtype=bool)
    # keep all predictor that we have priors for
    weight_mat = weight_mat.transpose()
    pp[weight_mat != no_pr_val] = True

    # for each gene, add the nS predictors of the list to possible predictors
    for ind in range(0, G):
        clr_order = np.argsort(clr_mat[ind, ])[::-1]
        pp[ind, clr_order[0: min(K, nS)]] = True

    out = []
    for ind in range(0, G):
        out1 = BBSRforOneGene(ind, X, Y, pp, weight_mat, nS)
        out.append(out1)

    return out


def BBSRforOneGene(ind, Xin, Yin, ppin, weight_matin, nSin):

    pp_i = ppin[ind, ]
    # create BestSubsetRegression input
    y = Yin[ind, ]
    x = Xin[pp_i, ].transpose()
    g = weight_matin[ind, pp_i]

    # experimental stuff
    spp = ReduceNumberOfPredictors(y, x, g, nSin)
    pp_i[np.where(pp_i)] = spp
    x = Xin[pp_i, :].transpose()
    g = weight_matin[ind, pp_i]

    betas = BestSubsetRegression(y, x, g)
    betas_resc = PredErrRed(y, x, betas)
    out_list = {}
    out_list['ind'] = ind
    out_list['pp'] = pp_i
    out_list['betas'] = betas
    out_list['betas_resc'] = betas_resc

    return out_list


def ReduceNumberOfPredictors(yin, xin, gin, n):
    K = xin.shape[1]
    if K <= n:
        return np.ones((1, K), dtype=bool)
    aa = np.ones((K, K))
    aa = np.diag(np.diag(aa))
    combos1 = CombCols(aa)
    combos2 = (aa == 1)
    combos = np.concatenate((combos2, combos1), axis=1)
    bics = ExpBICforAllCombos(yin, xin, gin, combos)
    bicsavg = np.sum(combos * bics, axis=1)
    ret = np.zeros((K), dtype=bool)
    ret[np.argsort(bicsavg)[2: n]] = True

    return ret


def CombCols(m):
    K = m.shape[1]
    ret = np.ones((m.shape[0], int(K * (K - 1) / 2)), dtype=bool)
    ret_col = 0
    for i in range(0, (K - 1)):
        for j in range(i + 1, K):
            ret[:, ret_col] = np.logical_or(m[:, i], m[:, j])
            ret_col = ret_col + 1
    return(ret)


def AllCombinations(k):
    # Create a boolean matrix with all possible combinations of 1:k.
    # Output has k rows and 2^k columns where each column is one
    # combination. Note that the first column is all FALSE and corresponds
    # to the null model
    if k < 1:
        sys.exit('No combinations for k < 1')

    N = 2**k
    out = np.zeros((k, N), dtype=bool)
    dumm1 = list(range(k))
    result = []
    for i in range(1, len(dumm1) + 1):
        result1 = list(itertools.combinations(dumm1, i))
        result1 = [list(ele) for ele in result1]
        result.extend(result1)
    for i in range(0, len(result)):
        out[result[i], i + 1] = True

    return(out)


def ExpBICforAllCombos(y, x, g, combos):
    # For a list of combinations of predictors do Bayesian linear regression,
    # more specifically calculate the parametrization of the inverse gamma
    # distribution that underlies sigma squared using Zellner's g-prior method
    # parameter g can be a vector. The expected value of the log of
    # sigma squared is used to compute expected values of BIC.
    # returns list of expected BIC values, one for each model.
    K = x.shape[1]
    N = x.shape[0]
    C = combos.shape[1]
    bics = []
    # if the first combination the null model?
    first_combo = 0
    if sum(combos[:, 0]) == 0:
        # bics.append(float(N * np.log(np.var(y))))
        bics.append(float('Inf'))
        first_combo = 1
    # shape parameter for the inverse gamma sigma squared would be drawn from
    shape = N / 2
    # compute digamma of shape here, so we can re-use it later
    dig_shape = scipy.special.digamma(shape)
    # pre-compute the crossproducts that we will need to solve for beta
    xtx = np.dot(x.transpose(), x)
    xty = np.dot(x.transpose(), y)
    # In Zellner's formulation there is a factor in the calculation of the rate
    # parameter: 1/(g+1)
    # Here we replace the factor with the appropriate matrix since g is
    # vector now
    var_mult = np.tile(np.sqrt(1 / (g + 1)).transpose(), K)
    var_mult = np.reshape(var_mult, (K, K))
    var_mult = var_mult.transpose()
    var_mult = np.multiply(var_mult, var_mult.transpose())

    for i in range(first_combo, C):
        comb = combos[:, i]
        x_tmp = x[:, comb]
        k = sum(comb)
        try:
            xtxuse = xtx[comb, :]
            xtxuse = xtxuse[:, comb]
            xtyuse = xty[comb]
            if k == 1:
                bhat = xtyuse / xtxuse
            else:
                bhat = np.linalg.solve(xtxuse, xtyuse)

            # sum of squares of residuals
            ssr = np.sum(np.square((y - np.squeeze(np.dot(x_tmp, bhat)))))
            # rate parameter for the inverse gamma sigma squared would be drawn
            # from our guess on the regression vector beta is all 0 for
            # sparse model
            var_multuse = var_mult[comb, :]
            var_multuse = var_multuse[:, comb]
            value1 = xtxuse * var_multuse
            rate = (ssr + (0 - bhat).dot(value1).dot(0 - bhat)) / 2
            # the expected value of the log of sigma squared based on the
            # parameterization of the inverse gamma by rate and shape
            exp_log_sigma2 = np.log(rate) - dig_shape

            # expected value of BIC
            bic1 = N * exp_log_sigma2 + k * np.log(N)
            bic1 = float(bic1)
            bics.append(bic1)
        except Exception:
            print('error in solve - system is computationally singular')
            bics.append(float('Inf'))
    bics = np.array(bics)
    return bics


def BestSubsetRegression(y, x, g):
    # Do best subset regression by using all possible combinations of
    # columns of x as predictor of y. Model selection criterion is BIC
    # using results of Bayesian regression with Zellner's g-prior
    # Args:
    # y: dependent variable
    # x: independent variable
    # g: value for Zellner's g-prior; can be single value or vector
    # Returns:
    # Beta vector of best model:

    K = x.shape[1]
    g = np.zeros(g.shape)
    combos = AllCombinations(K)
    bics = ExpBICforAllCombos(y, x, g, combos)

    not_done = True
    while not_done:
        best = np.argmin(bics)
        betas = np.zeros((K))
        if best > 0:
            x_tmp = x[:, combos[:, best]]
            xtx = np.dot(x_tmp.transpose(), x_tmp)
            xty = np.dot(x_tmp.transpose(), y)
            try:
                if sum(combos[:, best]) == 1:
                    bhat = xty / xtx
                else:
                    bhat = np.linalg.solve(xtx, xty)
                betas[combos[:, best]] = bhat
                not_done = False
            except Exception:
                bics[best] = float('Inf')
        else:
            not_done = False
    return betas


def PredErrRed(yin, xin, betain):
    # calculates the error reduction (measured by variance of residuals) of
    # each predictor - compare full model to model without that predictor
    K = xin.shape[1]
    pred = betain != 0
    P = sum(pred)

    # compute sigma^2 for full model
    residuals = yin - xin.dot(betain)
    sigma_sq_full = np.var(residuals)

    # this will be the output
    err_red = np.zeros((K))
    if P == 1:
        err_red[pred] = 1 - sigma_sq_full / np.var(yin)
        return err_red

    # one by one leave out each predictor and re-compute the model
    # with the remaining ones
    index = np.where(pred)[0]
    for i in index:
        pred_tmp = pred
        pred_tmp[i] = False
        x_tmp = xin[:, pred_tmp]
        xtx = np.dot(x_tmp.transpose(), x_tmp)
        xty = np.dot(x_tmp.transpose(), yin)
        if sum(pred_tmp) == 1:
            bhat = xty / xtx
        else:
            bhat = np.linalg.solve(xtx, xty)
        residuals = yin - x_tmp.dot(bhat)
        sigma_sq = np.var(residuals)

        err_red[i] = 1 - sigma_sq_full / sigma_sq
    return err_red
