import gzip

import numpy as np

import main

def read_groundtruth(fin, go, namespace, prot_ids):
    """
    Ground truth is assumed to be a two-columns file.

    A boolean matrix is returned with dimensionality: (n_proteins, n_terms)
    """
    pred = {}
    with gzip.open(fin, 'rt', encoding='utf-8') as f:
        for a in f:
            prot_id, term = a.strip().split('\t')
            term_namespace = go.ont[term]['namespace']

            # if prot_id not in prot_ids:
            #     continue

            if prot_id not in pred:
                pred[prot_id] = []
            # only save term corresponding with the specified namespace
            if namespace == term_namespace:
                term_id = go.term2index[namespace][term]
                pred[prot_id].append(term_id)
    # remove proteins without terms annotated
    pred = { prot_id: terms for prot_id, terms in pred.items() if len(terms) > 0 }

    n_benchmarks = len(pred)

    reduced_pred = {}
    for k in pred.keys():
        if k in prot_ids:
            reduced_pred[k] = pred[k]
    pred = reduced_pred

    num_proteins = len(pred.keys())
    proteins = []
    #
    # boolean matrix, for low memory footprint
    mat = np.zeros((num_proteins, go.counts[namespace]), dtype=np.bool)
    for prot_id in pred:
        mat[len(proteins), pred[prot_id]] = True
        proteins.append(prot_id)

    return n_benchmarks, proteins, mat

def read_predictions(fin, go, proteins, namespace, prot_ids):
    """Predictions are reading from a tab-separated file with three columns:
    prot_id, GO_term, predicted_confidence.
    """
    pred = {}
    with open(fin) as f:
        for a in f:
            protein, term, confidence = a.strip().split('\t')

            if protein not in prot_ids:
                continue

            # exclude terms that are not found in the GO
            try:
                term_namespace = go.ont[term]['namespace']
            except KeyError:
                term_namespace = None

            if term_namespace != namespace:
                continue
            term_id = go.term2index[namespace][term]
            confidence = float(confidence)
            # only save annotations on proteins previously found in the
            # benchmark
            if protein not in pred:
                pred[protein] = [[], []]
            pred[protein][0].append(term_id)
            pred[protein][1].append(confidence)
    # high precision to avoid round-off errors
    mat = np.zeros((len(proteins), go.counts[namespace]), dtype=np.float64)
    # fill matrix
    n_pred_proteins_in_benchmark = 0.
    for protein in proteins:
        prot_id = proteins.index(protein)
        if protein in pred:
            n_pred_proteins_in_benchmark += 1.
            prot_val = pred[protein]
            mat[prot_id, prot_val[0]] = prot_val[1]

    return n_pred_proteins_in_benchmark, mat

def cast_predictions_into_dict(fin, namespace):
    """
    """
    go = main.get_go()
    pred = {}
    with gzip.open(fin, 'rt', encoding='utf-8') as f:
        for a in f:
            protein, term, confidence = a.strip().split('\t')

            # exclude terms that are not found in the GO
            try:
                term_namespace = go.ont[term]['namespace']
            except KeyError:
                term_namespace = None

            if term_namespace != namespace:
                continue

            if protein not in pred:
                pred[protein] = []
            pred[protein].append((term, float(confidence)))

    return pred


def _get_prot_ids(fin, go, namespace):
    """
    Obtain the prot_id from
    """
    prot_ids = set({})
    with gzip.open(fin, 'rt', encoding='utf-8') as f:
        for a in f:
            raw = a.strip().split('\t')
            prot_id = raw[0]
            term = raw[1]

            # exclude terms that are not found in the GO
            try :
                term_namespace = go.ont[term]['namespace']
            except KeyError:
                continue

            if term_namespace != namespace:
                continue
            prot_ids.add(prot_id)

    return prot_ids

def get_common_prots(true_fin, pred, go, namespace):
    # protein ids in the so-called benchmarks
    true_prot_ids = _get_prot_ids(true_fin, go, namespace)
    # protein ids in predictions
    pred_prot_ids = set(pred.keys())
    common_prot_ids = set.intersection(true_prot_ids, pred_prot_ids)

    return common_prot_ids
