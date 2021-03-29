import os
import gzip

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# only for human proteins
UNIP2ID = '{}/uniprot_ac_to_id_9606.map.gz'.format(CURRENT_DIR)
CAFA2ID = '{}/sp_species.9606.map.gz'.format(CURRENT_DIR)

def map(ypred):
    """This method makes the following mapping:

    UniProtKB id -> Protein id -> CAFA id


    Args
    ----

    ypred : python dict containing UniProtKB accessions as keys and GO
    terms as values.

    Returns
    -------

    mapped_ypred: python dict containing CAFA accessions as keys and GO terms
    as values.
    """
    unip2id = {}
    id2cafa = {}

    # read files with the two files used for mappings

    # First file maps UniProtKB accessions to protein ids
    with gzip.open(UNIP2ID, 'rt', encoding='utf-8') as f:
        for a in f:
            acc, _, prot_id = a.strip().split('\t')
            unip2id[acc] = prot_id

    # Second file maps protein ids to CAFA accessions
    with gzip.open(CAFA2ID, 'rt', encoding='utf-8') as f:
        for a in f:
            cafa_id, prot_id = a.strip().split('\t')
            id2cafa[prot_id] = cafa_id

    # map protein accessions passed as argument to CAFA accessions
    mapped_ypred = {}
    for k, v in ypred.items():
        if k not in unip2id:
            continue

        mapped_id = unip2id[k]
        if mapped_id not in id2cafa:
            continue

        mapped_id = id2cafa[mapped_id]
        mapped_ypred[mapped_id] = v

    return mapped_ypred
