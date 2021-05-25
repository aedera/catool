#!/usr/bin/env python3
import sys
import os

import numpy as np
import numba as nb

#from precRec import benchmark
from . import onto
from . import inou
from . import utils
from .conversion import mapper
from .metrics import precision_recall_curve


CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# obo file used for calculations
OBO_FILE = '{}/go.obo'.format(CURRENT_DIR)

TRUE_FIN1 = '%s/true_y_t1.tsv.gz' % (CURRENT_DIR) # bpo_HUMAN_type1.txt
TRUE_FIN2 = '%s/true_y_t2.tsv.gz' % (CURRENT_DIR) # bpo_HUMAN_type2.txt
TRUE_FIN3 = '%s/true_y_t3.tsv.gz' % (CURRENT_DIR) # bpo_HUMAN_type1.txt and bpo_HUMAN_type2.txt

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

def all_scores(pred,
               mode,
               namespace,
               ftype=1,
               warning_CAFA3_wrong_terms=False):
    """This method calculates performance metrics by comparing predictions
    with groundtruths.

    Args:

    pred : a dict whose keys are protein identifier and values are tuple of
    two elements.  The first element is a GO term and the second element is a
    predicted value \in [0,1]

    mode : a string indicating "full" or "partial" mode.

    namespace : a string indicating sub-ontology: "biological_process",
    "cellular_component", or "molecular_function".

    ftype : an integer indicating which type of benchmarks to use.  At the
    moment, this variable takes two values: 1 (bpo_HUMAN_type1.txt) and 2
    (bpo_HUMAN_type1.txt + bpo_HUMAN_type2.txt).

    Returns:

    results : a matrix with shape (n_thresholds, 4).  The four columns
    corresponds to: threshold, f1, precision, recall
    """
    go = get_go() # get ontology

    # convert UniProt accessions into CAFA ids
    pred = mapper.map(pred)

    # define benchmark file
    if ftype == 1:
        true_fin = TRUE_FIN1 # type1
    elif ftype == 2:
        true_fin = TRUE_FIN2 # type2
    elif ftype == 3:
        true_fin = TRUE_FIN3 # type1 + type2
    else:
        sys.exit("File type unknown.")

    # obtained protein identifiers that appear in the benchmarks
    common_prots = inou.get_common_prots(true_fin, pred, go, namespace)

    # read data from input files
    # true \in (n_proteins, n_onto_specific_terms)
    n_benchmarks, proteins, y_true = inou.read_groundtruth(
        true_fin, go, namespace, common_prots, warning_CAFA3_wrong_terms)

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

def f1max_score(pred, mode, namespace, ftype=1, warning_CAFA3_wrong_terms=False):
    """ Return maximum F1 score"""

    scores = all_scores(pred, mode, namespace, ftype, warning_CAFA3_wrong_terms)
    return max(scores[:,1])
