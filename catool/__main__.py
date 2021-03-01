#!/usr/bin/env python3
import sys
import os

from catool import main
from catool import inou
from catool import utils

if __name__ == '__main__':
    pred_fin = sys.argv[1] # predictions
    mode = sys.argv[2] # full or partial
    namespace = sys.argv[3]

    go = main.get_go()
    # read preditions from file
    pred = inou.cast_predictions_into_dict(pred_fin, go, namespace)

    #print(main.f1max_score(pred, mode, namespace))
    results = main.run(pred, mode, namespace)
    utils.results2string(results, mode, namespace)
