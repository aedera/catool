#!/usr/bin/env python3
import os
import gzip

import numpy as np
import numba as nb

from . import onto
from . import main

@nb.njit(parallel=True)
def propagate_terms(m, branches):
    for i in nb.prange(m.shape[0]): # loop over proteins
        for j in np.where(m[i,:])[0]: # loop over terms
            for k in np.where(branches[j,:])[0]:
                if m[i, j] > m[i, k]:
                    m[i, k] = m[i, j]

def create_ancestors_cache(go, namespace, fout):
    n_terms = go.counts[namespace]
    ancestors = np.zeros((n_terms, n_terms), dtype=np.bool)

    for i in range(n_terms):
        term = go.index2term[namespace][i]
        propagations = go.get_ancestors(term)
        # flat propagations
        propagations = set([t for b in propagations for t in b])
        # don't include root terms and leaf term
        propagations = propagations.difference({
            term,
            onto.BIOLOGICAL_PROCESS,
            onto.CELLULAR_COMPONENT,
            onto.MOLECULAR_FUNCTION})
        propagations = [go.term2index[namespace][p] for p in propagations]
        ancestors[i, propagations] = True

    np.savez_compressed(fout, a=ancestors)

def create_cache_files(log=False):
    go = main.get_go()

    for namespace in onto.NAMESPACES.values():
        fout = os.path.dirname(os.path.realpath(__file__)) + '/' + namespace + '.npz'

        if not os.path.exists(fout):
            create_ancestors_cache(go, namespace, fout)
            if log:
                print(namespace + '.npz was created.')


def predictions_into_a_matrix(pred, benchmark_prots, namespace):
    """Create matrix (n_prots, n_term) whose values are the predicted
probabilities."""
    go = main.get_go()
    # pred = { k:pred[k] for k in pred.keys() if k in common_prots  }
    # n_predicted_proteins_in_benchmark = len(pred.keys())

    mat = np.zeros((len(benchmark_prots), go.counts[namespace]), dtype=np.float64)
    # fill matrix
    n_pred_proteins_in_benchmark = 0.
    for prot_id, protein in enumerate(benchmark_prots):
        if protein in pred:
            n_pred_proteins_in_benchmark += 1.
            # loop over pairs (predicted terms, probability)
            for term, prob in pred[protein]:
                if term not in go.ont.keys(): # discard terms not present in go
                    continue
                term = go.term2index[namespace][term]
                mat[prot_id, term] = prob

    return n_pred_proteins_in_benchmark, mat


def results2string(results, mode, namespace):
    for i in range(results.shape[0]):
        print("{}\t{}\t{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}".format(
            mode,
            namespace,
            results[i,0],
            results[i,1],
            results[i,2],
            results[i,3]))


if __name__ == '__main__':
    create_caches(log=True)
