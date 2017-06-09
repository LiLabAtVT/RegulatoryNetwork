import math
import numpy as np


def HalfThr(A, b, xinput, MAX_ITERS, regulatorgene_num, s):
    # The iterative half thresholding algorithm to solve the L1/2
    # regularization model of sparse optimization
    # For the details, one can refer to Z.Xu, X.Chang, F.Xu,
    # and H.Zhang. L1/2 regularization:
    # a thresholding representation theory and a fast solver.
    # IEEE Transactions on Neural Networks and Learning Systems,
    # 23:1013"C1027,2012

    f = np.full((MAX_ITERS, 1), 1, dtype=float)
    k = 0
    x = xinput.copy()
    mu = 1  # mu is double of step-size v, which chose as 1/2
    Va1 = math.sqrt(96) / 9
# Va2 = 54**(1/3)/4
    Bu1 = mu * A.T.dot(b)
    Bu2 = mu * A.T.dot(A)
    while k < MAX_ITERS - 1:
        i = 0
        Bu = x + Bu1 - Bu2.dot(x)
        BuO = sorted(abs(Bu))
        BuV = BuO[int(regulatorgene_num - s - 1)] ** (3 / 2)
        lambda1 = Va1 * BuV
        criterion = BuO[int(regulatorgene_num - s - 1)]
        while i <= regulatorgene_num - 1:
            if abs(Bu[i]) > criterion:
                phi = math.acos(lambda1 * mu / 8 * (abs(Bu[i]) / 3) ** (-1.5))
                x[i] = 2 * Bu[i] * (1 + math.cos(2 / 3 * (math.pi - phi))) / 3
            else:
                x[i] = 0
            i = i + 1
        f[k + 1] = np.linalg.norm(A.dot(x) - b)**2
        k = k + 1
    hist = f
    return x, hist


def HardThr(A, b, x1, MAX_ITERS, regulatorgene_num, s):
    # The iterative hard thresholding algorithm to solve the L0 regularizaion
    # model of sparse optimization.
    # For the details, one can refer to T. Blumensath and M. E. Davies.
    # Iterative thresholding for sparse approximations.
    # Journal of Fourier Analysis and Applications, 14:629¨C654, 2008.
    f = np.full((MAX_ITERS, 1), 1, dtype=float)
    k = 0
    x = x1.copy()
    Ause = A
    buse = b
    mu = 1  # mu is double of step-size v, which chose as 1/2
    Bu1 = mu * Ause.T.dot(buse)
    Bu2 = mu * Ause.T.dot(Ause)
    while k < MAX_ITERS - 1:
        i = 0
        Bu = x + Bu1 - Bu2.dot(x)
        BuO = sorted(abs(Bu))
        criterion = BuO[int(regulatorgene_num - s - 1)]
        while i <= regulatorgene_num - 1:
            if abs(Bu[i]) > criterion:
                x[i] = Bu[i]
            else:
                x[i] = 0
            i = i + 1
        f[k + 1] = np.linalg.norm(Ause.dot(x) - buse) ** 2
        k = k + 1
    hist = f
    return x, hist


def SoftThr(A, b, x1, MAX_ITERS, regulatorgene_num, s):
    # The iterative soft thresholding algorithm to solve the L1
    # regularization model (LASSO) of sparse optimization.
    # For the details, one can refer to I. Daubechies, M. Defrise,
    # and C. De Mol.
    # An iterative thresholding algorithm for linear inverse problems
    # with a sparsity constraint.
    # Communications on Pure and Applied Mathematics, 57:1413¨C1457, 2004.
    f = np.full((MAX_ITERS, 1), 1, dtype=float)
    k = 0
    x = x1.copy()
    Ause = A
    buse = b
    mu = 1  # mu is double of stepsize v, which choosed as 1/2
    Bu1 = mu * Ause.T.dot(buse)
    Bu2 = mu * Ause.T.dot(Ause)
    while k < MAX_ITERS - 1:
        i = 0
        Bu = x + Bu1 - Bu2.dot(x)
        BuO = sorted(abs(Bu))
        criterion = BuO[int(regulatorgene_num - s - 1)]
        while i <= regulatorgene_num - 1:
            if abs(Bu[i]) > criterion:
                x[i] = Bu[i] - criterion * np.sign(Bu[i])
            else:
                x[i] = 0
            i = i + 1
        f[k + 1] = np.linalg.norm(Ause.dot(x) - buse) ** 2
        k = k + 1
    hist = f
    return x, hist
