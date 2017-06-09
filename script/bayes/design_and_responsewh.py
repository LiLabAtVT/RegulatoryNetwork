#  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
# /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input: cM - colMap
#           r - expression matrix (assumed already normalized)
#           param - a vector of 2 params:
# 1-  c1 - cutoff1: represent the maximal time interval
#      allowed for time series (ts) data
# 2- 'trivial' or 'time_delayed' (param to choose the type of design matrix)
# 3- 'consecutive' or 'all_intervals' (determine if consecutive time
#       measurements [no longer than c1],
#       or all permutations of time measurements [up to c1]
#       respectively)
# 4- 'TRUE' or 'FALSE' use_t0_as_steady_state
# 5- 'TRUE' or 'FALSE' use_delt_bigger_than_cutoff_as_steady_state
#
# output:
#            steadyStateDesignMat, timeSeriesDesignMat
import numpy as np


def get_usr_chosen_design(
    cM, r, rname, delT_min, delT_max, time_delayed, all_intervals,
        use_t0_as_steady_state,
        use_delt_bigger_than_cutoff_as_steady_state):

    delT_min_vec = np.empty(len(cM['condName']))
    delT_min_vec.fill(delT_min)
    delT_max_vec = np.empty(len(cM['condName']))
    delT_max_vec.fill(delT_max)

    delT_vec = cM['del.t']
    isTs_vec = cM['isTs']
    eq_idx = [i for i, x in enumerate(isTs_vec) if x == 'FALSE' or x is False]
    ts_idx = [i for i, x in enumerate(isTs_vec) if x == 'TRUE']

    delT_vec = [delT_vec[i] for i in ts_idx]
    delT_min_vec = [delT_min_vec[i] for i in ts_idx]
    delT_max_vec = [delT_max_vec[i] for i in ts_idx]
# set delT_vec
# 0 - last time measurement in ts
# >0 - first and middle time measurements in ts
# following line make 0s indicate first time measurement in ts
    for ind, item in enumerate(delT_vec):
        if item == 'NA':
            delT_vec[ind] = 0
    delT_vec = [float(item) for item in delT_vec]
    delT_vec_trivial = delT_vec
# following 2 lines make 0s indicate last time measurement in ts
    if len(delT_vec) > 0:
        delT_vec1 = delT_vec[1:len(delT_vec)]
        delT_vec[0:-1] = delT_vec1
        delT_vec[-1] = 0
# data for steady state
    rSS = r[:, eq_idx]
    rSS_colname = [rname[i] for i in eq_idx]
# if we have time series experiments
    if len(ts_idx) > 0:
        if len(eq_idx) > 0:
            rTS = np.delete(r, eq_idx, axis=1)
            rTSname = list(set(rname) - set(rSS_colname))
        else:
            rTS = r
            rTSname = rname
        eq_idx_pseudo = []
# get ts starting conditions, we treat these as equibibrium
        if use_t0_as_steady_state:
            eq_idx_pseudo = [
                i for i,
                x in enumerate(delT_vec_trivial) if x == 0]
# get ts conditions with larger than c1 delt, we treat these as equibibrium
        if use_delt_bigger_than_cutoff_as_steady_state:
            eq_idx_pseudo = eq_idx_pseudo + \
                [i for i,
                    x in enumerate(delT_vec_trivial) if x > max(delT_max_vec)]
# create design matrix for steady state
        if len(eq_idx) > 0:
            DesignMatSS = np.concatenate((
                rSS, rTS[:, eq_idx_pseudo]), axis=1)
        else:
            DesignMatSS = rTS[:, eq_idx_pseudo]

        if time_delayed:
            ind1 = [i for i, x in enumerate(delT_vec) if (x != 0)]
            ind2 = []
            for i in range(0, len(delT_vec)):
                if delT_vec[i] <= delT_max_vec[i]:
                    ind2.append(i)
            ind = list(set(ind1).intersection(ind2))
            DesignMatTS = rTS[:, ind]
        else:
            ind1 = [i for i, x in enumerate(delT_vec_trivial) if (x != 0)]
            ind2 = []
            for i in range(0, len(delT_vec_trivial)):
                if delT_vec_trivial[i] <= delT_max_vec[i]:
                    ind2.append(i)
            ind = list(set(ind1).intersection(ind2))
            DesignMatTS = rTS[:, ind]
    else:
            DesignMatSS = rSS
            DesignMatTS = np.zeros((DesignMatSS.shape[0], 1))
    output = []
    output.append(DesignMatSS)
    output.append(DesignMatTS)
    return(output)


#  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
# /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input: cM - colMap
#           r - ratios matrix (assumed already normalized)
#           param - a vector of 5 params:
# 1-  c1 - cutoff1:
# represent the maximal time interval allowed for time series (ts) data
# 2- 'trivial','time_difference','rate','inf_1', or 'inf_1_all_intervals'
# 3- tau
# 4- TRUE or FALSE use_t0_as_steady_state
#  5- TRUE or FALSE use_delt_bigger_than_cutoff_as_steady_state
#
# output:
#            steadyStateResponseMat timeSeriesResponseMat


def get_usr_chosen_response(
    cM, r, rname, delT_min, delT_max, method, tau,
        use_t0_as_steady_state,
        use_delt_bigger_than_cutoff_as_steady_state):

    delT_min_vec = np.empty(len(cM['condName']))
    delT_min_vec.fill(delT_min)
    delT_max_vec = np.empty(len(cM['condName']))
    delT_max_vec.fill(delT_max)

    delT_vec = cM['del.t']
    isTs_vec = cM['isTs']
    eq_idx = [i for i, x in enumerate(isTs_vec) if x == 'FALSE' or x is False]
    ts_idx = [i for i, x in enumerate(isTs_vec) if x == 'TRUE']

    delT_vec = [delT_vec[i] for i in ts_idx]
    delT_min_vec = [delT_min_vec[i] for i in ts_idx]
    delT_max_vec = [delT_max_vec[i] for i in ts_idx]
# set delT_vec
# 0 - last time measurement in ts
# >0 - first and middle time measurements in ts
# following line make 0s indicate first time measurement in ts
    for ind, item in enumerate(delT_vec):
        if item == 'NA':
            delT_vec[ind] = 0
    delT_vec = [float(item) for item in delT_vec]
    delT_vec_trivial = delT_vec
# following 2 lines make 0s indicate last time measurement in ts
    if len(delT_vec) > 0:
        delT_vec1 = delT_vec[1:len(delT_vec)]
        delT_vec[0:-1] = delT_vec1
        delT_vec[-1] = 0
# data for steady state
    rSS = r[:, eq_idx]
    rSS_colname = [rname[i] for i in eq_idx]
# if we have time series experiments
    if len(ts_idx) > 0:
        if len(eq_idx) > 0:
            rTS = np.delete(r, eq_idx, axis=1)
            rTSname = list(set(rname) - set(rSS_colname))
        else:
            rTS = r
            rTSname = rname
        eq_idx_pseudo = []
# get ts starting conditions, we treat these as equibibrium
        if use_t0_as_steady_state:
            eq_idx_pseudo = [
                i for i,
                x in enumerate(delT_vec_trivial) if x == 0]
# get ts conditions with larger than c1 delt, we treat these as equibibrium
        if use_delt_bigger_than_cutoff_as_steady_state:
            eq_idx_pseudo = eq_idx_pseudo + \
                [i for i,
                    x in enumerate(delT_vec_trivial) if x > max(delT_max_vec)]

        if len(eq_idx) > 0:
            response_matrixSS = np.concatenate((
                rSS, rTS[:, eq_idx_pseudo]), axis=1)
        else:
            if len(eq_idx_pseudo) > 0:
                response_matrixSS = rTS[:, eq_idx_pseudo]
            else:
                response_matrixSS = np.array([])

        init_ind1 = [i for i, x in enumerate(delT_vec) if (x != 0)]
        init_ind2 = []
        for i in range(0, len(delT_vec)):
            if delT_vec[i] <= delT_max_vec[i]:
                init_ind2.append(i)
        init_ind = list(set(init_ind1).intersection(init_ind2))
        init_ind = np.array(init_ind)
        boundary_ind = init_ind + 1
# finished response matrices now go on to response
        if method == 'trivial':
            response_matrixTS = rTS[:, boundary_ind]
        elif method == 'time_difference':
            response_matrixTS = rTS[:, boundary_ind] - rTS[:, init_ind]
        elif method == 'rate':
            a = 1 / np.array(delT_vec)[init_ind]
            b = np.transpose((rTS[:, boundary_ind] - rTS[:, init_ind]))
            for i in range(0, len(a)):
                b[i, :] = a[i] * b[i, :]
            response_matrixTS = np.transpose(b)
        elif method == 'inf_1':
            a = tau / np.array(delT_vec)[init_ind]
            b = np.transpose((rTS[:, boundary_ind] - rTS[:, init_ind]))
            for i in range(0, len(a)):
                b[i, :] = a[i] * b[i, :]
            response_matrixTS = np.transpose(b) + rTS[:, init_ind]
        else:
            raise ValueError('unknown response ')
    else:
        response_matrixSS = rSS
        response_matrixTS = np.zeros((response_matrixSS.shape[0], 1))
    output = []
    output.append(response_matrixSS)
    output.append(response_matrixTS)
    return(output)

#  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
# /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input:
#   1- dMSS: design matrix steady state
#   2- dMTS: design matrix time series
#   3- rMSS: response matrix steady state
#   4- dMTS: response matrix time series
#   5- param: what final design matrix? choose from all, ts, or ss

# output:
#   resopnse and corresponding design matrices


def make_final_design_and_response_matrix(
        dMSS, dMTS, rMSS, rMTS, para, cS, r, tf_names, make_des_red_exp):
    if para == 'all':
        if np.any(rMTS):
            final_response_matrix = np.concatenate((
                rMSS, rMTS), axis=1)
        else:
            final_response_matrix = rMSS
        if np.any(dMTS):
            final_design_matrix = np.concatenate((
                dMSS, dMTS), axis=1)
        else:
            final_design_matrix = dMSS
    elif para == 'ts':
        final_response_matrix = rMTS
        final_design_matrix = dMTS
    elif para == 'ss':
        final_response_matrix = rMSS
        final_design_matrix = dMSS
    else:
        raise ValueError('unknown response ')

# the original r code has a part of check biclusters,
# I don't really know of that part, so I ignore them

    resp_idx = np.array(
        [range(0, final_response_matrix.shape[1])] *
        final_response_matrix.shape[0])
    output = []
    output.append(final_response_matrix)
    output.append(final_design_matrix)
    output.append(resp_idx)
    return(output)
