import numpy as np
from sklearn.metrics import mutual_info_score


def discretize(x, nbins):
    # N = x.size
    xmin = np.min(x)
    xmax = np.max(x)
    tiny = np.finfo(np.float64).eps * (xmax - xmin)
    X_disc = np.floor((x - xmin) / (xmax - xmin + tiny) * nbins)
    X_disc.astype(int)

    return X_disc


# compute MI each column in x to each column in y
# you'll want the larger matrix to be x
def mi(x, y, nbins=10, cpu_n=1, perm_mat=None):
    if perm_mat is not None:
        xuse = np.zeros((x.shape))
        yuse = np.zeros((y.shape))
        for i in range(0, x.shape[1]):
            xuse[:, i] = x[perm_mat[:, i], i]
        for i in range(0, y.shape[1]):
            yuse[:, i] = y[perm_mat[:, i], i]
        x = xuse
        y = yuse

    # discretized the columns
    x = np.apply_along_axis(discretize, 0, x, nbins)
    y = np.apply_along_axis(discretize, 0, y, nbins)

    m = x.shape[1]
    n = y.shape[1]
    # ret = np.zeros((m, n))
    mi = np.zeros((m, n))
    # if cpu_n == 1:
    for i in range(0, m):
        for j in range(0, n):
            xx = x[:, i]
            yy = y[:, j]
            mi[i, j] = calc_MI(xx, yy, nbins)
    mi = np.round(mi, 4)
    return(mi)
    # aa = np.arange(0, m, 1)
    # out = np.array_split(aa, cpu_n)


def calc_MI(xx, yy, bins):
    c_xy = np.histogram2d(xx, yy, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi


def toZscore(x):
    out = (x - np.nanmean(x)) / np.nanstd(x)
    return out


def mixedCLR(mi_stat, mi_dyn):
    z_r_dyn = np.apply_along_axis(toZscore, 0, mi_stat)
    z_r_dyn[z_r_dyn < 0] = 0

    z_c_mix = mi_dyn
    z_c_mix = np.apply_along_axis(toZscore, 0, mi_stat)
    z_c_mix[z_c_mix < 0] = 0

    out = z_r_dyn**2 + z_c_mix**2
    out = np.sqrt(out)
    return out
