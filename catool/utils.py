#!/usr/bin/env python3
import os
import gzip

import numpy as np
import numba as nb

from . import onto

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
    OBO_FILE = '{}/go.obo'.format(os.path.dirname(os.path.realpath(__file__)))
    go = onto.Ontology(OBO_FILE, with_rels=True, include_alt_ids=False)
    for namespace in onto.NAMESPACES.values():
        fout = os.path.dirname(os.path.realpath(__file__)) + '/' + namespace + '.npz'
        create_ancestors_cache(go, namespace, fout)
        if log:
            print(namespace + '.npz was saved.')

if __name__ == '__main__':
    create_caches(log=True)
