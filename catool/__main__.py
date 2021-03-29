#!/usr/bin/env python3
import sys
import os

from . import main
from . import inou
from . import utils

if __name__ == '__main__':
    pred_fin = sys.argv[1] # input file with predicted GO terms
    mode = sys.argv[2] # full or partial
    namespace = sys.argv[3] # biological_process
    ftype = int(sys.argv[4]) # type_1, type_2, or type_3 (=type_1 + type_2)

    go = main.get_go()

    # read input file
    pred = inou.cast_predictions_into_dict(pred_fin, namespace)

    scores = main.all_scores(pred, mode, namespace, ftype)
    utils.results2string(scores, mode, namespace)
