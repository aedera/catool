import numpy as np
import numba as nb

EPSILON = np.finfo(float).eps

#@nb.njit(parallel=True, fastmath=True)
@nb.njit(parallel=True)
def precision_and_recall_at_thr(true, pred, n_predicted_proteins, thr):
    """#
    true: (n_proteins, n_terms)
    pred: (n_proteins, n_terms)

    return:
    precision pr = true_pos / (true_pos + false_pos)
    recall    re = true_pos / (true_pos + false_neg)

    \sum_i^m pr_i
    \sum_i^N rc_N
    """
    pr = 0.
    re = 0.
    n_retrieved = 0.
    for i in nb.prange(true.shape[0]): # loop over proteins
        ### PROTEIN INDEX
        #i = 2015
        #i = 1981 # T96060010945
        # i = 2061 # T96060010935
        # i = 2126 # T96060001131
        ####################

        tp = 0. # true positives
        retrieved = 0.
        relevant = 0.
        for j in nb.prange(true.shape[1]): # loop over terms
            val = (pred[i,j] >= thr)
            tp += true[i,j] * val
            relevant += true[i,j]
            retrieved += val # count the number of pred >= thr
        pr += tp / (retrieved + EPSILON)
        re += tp / (relevant + EPSILON)
        # count proteins with at least one prediction
        if retrieved > 0:
            n_retrieved += 1.

        ######
        # break
    ############
    #print(n_retrieved)
    #print(n_predicted_proteins)

    return pr / n_retrieved, re / n_predicted_proteins

@nb.njit(parallel=True)
def precision_recall_curve(true, pred, n_predicted_proteins, out):
    for i in nb.prange(out.shape[0]):
        thr = 0.01 * i # prediction threshold
        pr, re = precision_and_recall_at_thr(true, pred, n_predicted_proteins, thr)
        # if thr == 0.01:
        #     breakpoint()

        if re > 1:
            re = 1.
        f1 = 2. * ((pr * re) / (pr + re + EPSILON)) # f1-measure
        out[i, 0] = thr
        out[i, 1] = f1
        out[i, 2] = pr
        out[i, 3] = re
