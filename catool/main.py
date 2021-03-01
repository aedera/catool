#!/usr/bin/env python3
import sys
import os

import numpy as np
import numba as nb

#from precRec import benchmark
import onto
import inou
import utils
from metrics import precision_recall_curve

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# obo file used for calculations
OBO_FILE = '{}/go.obo'.format(CURRENT_DIR)
TRUE_FIN = '%s/true_y.tsv.gz' % (CURRENT_DIR)

global __GO__
__GO__ = [None]

def _decompress_obofile():
    import gzip
    import shutil
    with gzip.open(OBO_FILE + ".gz", 'r') as f_in:
        with open(OBO_FILE, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def get_go():
    """
    The GO is a global variable to avoid multiple reading of the obo file.
    """
    # decompress obo file if it wasn't yet
    if not os.path.exists(OBO_FILE):
        _decompress_obofile()
    # create global variable
    if __GO__[0] is None:
        __GO__[0] = onto.Ontology(OBO_FILE, with_rels=True, include_alt_ids=False)
    return __GO__[0]

def run(pred, mode, namespace):
    """This method calculates performance metrics by comparing predictions with
groundtruths.

    Input:
    =====
    - pred: a dict whose keys are protein identifier and values are tuple of
      two elements. The first element is a GO term and the second element is a
      predicted value \in [0,1]

    - mode: a string indicating "full" or "partial" mode.

    - namespace: a string indicating "biological_process",
      "cellular_component", or "molecular_function".

    Output:
    ======
    - results: a matrix (n_thresholds, 4): columns are: thr, f1, precision,
      recall.

    """
    go = get_go()

    # obtained protein identifieres that appear in the benchmarks
    common_prots = inou.get_common_prots(TRUE_FIN, pred, go, namespace)

    # read data from input files
    # true \in (n_proteins, n_onto_specific_terms)
    n_benchmarks, proteins, y_true = inou.read_groundtruth(
        TRUE_FIN, go, namespace, common_prots)

    n_predicted_proteins_in_benchmark, y_pred = utils.predictions_into_a_matrix(pred, proteins, namespace)

    # propagate terms based on the topological structure of the GO
    #
    # create cache files containing ancestor terms
    utils.create_cache_files(log=False)
    # read cache file with ancestors
    anc_fin = os.path.dirname(os.path.realpath(__file__)) + '/' + namespace + '.npz'
    ancestors = np.load(anc_fin)['a']

    # propagate terms to include their ancestors
    utils.propagate_terms(y_true, ancestors)
    utils.propagate_terms(y_pred, ancestors)
    branches = None # free memory

    # Exclude root terms
    root_term = onto.namespace2go[namespace]
    root_term_id = go.term2index[namespace][root_term]
    y_true[:, root_term_id] = False
    y_pred[:, root_term_id] = 0.

    # Calculate metrics
    results = np.empty((len(range(0, 101)), 4))
    baseline = n_benchmarks if mode == 'full' else n_predicted_proteins_in_benchmark
    precision_recall_curve(y_true, y_pred, baseline, results)

    return results

def f1max_score(pred, mode, namespace):
    """ Return maximum F1 score"""
    results = run(pred, mode, namespace)
    return max(results[:,1])
