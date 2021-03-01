import os
import gzip

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# only for human proteins
UNIP2ID = '{}/uniprot_ac_to_id_9606.map.gz'.format(CURRENT_DIR)
CAFA2ID = '{}/sp_species.9606.map.gz'.format(CURRENT_DIR)

def map(y_pred):
    """
    preds = list of UniProtKB accessions
    """
    unip2id = {}
    id2cafa = {}

    # read files
    with gzip.open(UNIP2ID, 'rt', encoding='utf-8') as f:
        for a in f:
            acc, _, prot_id = a.strip().split('\t')
            unip2id[acc] = prot_id
    #
    with gzip.open(CAFA2ID, 'rt', encoding='utf-8') as f:
        for a in f:
            cafa_id, prot_id = a.strip().split('\t')
            id2cafa[prot_id] = cafa_id

    # map ids
    new_pred = {}
    for k, v in y_pred.items():
        if k not in unip2id:
            continue

        mapped_id = unip2id[k]
        if mapped_id not in id2cafa:
            continue

        mapped_id = id2cafa[mapped_id]
        new_pred[mapped_id] = v


    return new_pred
