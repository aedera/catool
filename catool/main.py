#!/usr/bin/env python3
import sys
import os

import numpy as np
import numba as nb

#from precRec import benchmark
from . import onto
from . import inou
from . import utils
from .metrics import precision_recall_curve

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# obo file used for calculations
OBO_FILE = '{}/go.obo'.format(CURRENT_DIR)

def run(true_fin, pred_fin, namespace, create_cache_files=False):
    go = onto.Ontology(OBO_FILE, with_rels=True, include_alt_ids=False)

    # obtained protein identifieres that appear in the benchmarks
    common_prots = inou.get_common_prots(true_fin, pred_fin, go, namespace)

    # read data from input files
    # true \in (n_proteins, n_onto_specific_terms)
    n_benchmarks, proteins, true = inou.read_groundtruth(
        true_fin, go, namespace, common_prots)
    # pred \in (n_proteins, n_onto_specific_terms)
    n_predicted_proteins_in_benchmark, pred = inou.read_predictions(
        pred_fin, go, proteins, namespace, common_prots)

    #
    # propagate terms based on the topological structure of the GO
    #
    if create_cache_files:
        # create cache files containing ancestor terms
        utils.create_cache_files(log=False)
    # read cache file with ancestors
    anc_fin = os.path.dirname(os.path.realpath(__file__)) + '/' + namespace + '.npz'
    ancestors = np.load(anc_fin)['a']

    # propagate terms to include their ancestors
    utils.propagate_terms(true, ancestors)
    utils.propagate_terms(pred, ancestors)
    branches = None # free memory

    # Exclude root terms
    root_term = onto.namespace2go[namespace]
    root_term_id = go.term2index[namespace][root_term]
    true[:, root_term_id] = False
    pred[:, root_term_id] = 0.

    #
    # Calculate metrics for both modes: full & partial
    #
    # Full
    results = np.empty((len(range(0, 101)), 4))
    precision_recall_curve(true, pred, n_benchmarks, results)
    for i in range(results.shape[0]):
        print("{}\t{}\t{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}".format(
            'full', namespace, results[i,0], results[i,1], results[i,2], results[i,3]))
    # Partial
    results = np.empty((len(range(0, 101)), 4))
    precision_recall_curve(true, pred, n_predicted_proteins_in_benchmark, results)
    for i in range(results.shape[0]):
        print("{}\t{}\t{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}".format(
            'partial', namespace, results[i,0], results[i,1], results[i,2], results[i,3]))
