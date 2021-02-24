#!/usr/bin/env python3
import sys
import os

from catool import main
from catool import utils

if __name__ == '__main__':
    opname = sys.argv[1]

    if opname == 'run':
        true_fin = sys.argv[2]
        pred_fin = sys.argv[3] # predictions
        namespace = sys.argv[4]

        main.run(true_fin, pred_fin, namespace)
    else:
        utils.create_cache_files(True)
