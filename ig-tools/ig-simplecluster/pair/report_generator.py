__author__ = 'mactep'

import os
import json

from pair.data_loader import *
from pair.data_generator import make_cells
import common


def sorted_keys(cls):
    f_sort_key = lambda x: int(x[x.rfind(common.CLUSTER_SPLT) +
                                 len(common.CLUSTER_SPLT):]) if common.TRASH_NAME not in x else 100500
    return sorted(cls.keys(), key=f_sort_key)


def process(heavy_dir, light_dir, fix):
    heavy = load_clusters(os.path.abspath(heavy_dir))
    light = load_clusters(os.path.abspath(light_dir))

    if fix:
        heavy = fix_suffix(heavy)
        light = fix_suffix(light)

    hcons = load_consensus(heavy_dir)
    lcons = load_consensus(light_dir)

    heavy_info = load_info(heavy_dir)
    light_info = load_info(light_dir)

    cells, names = make_cells(heavy, light)

    # Result dict
    result = {"info": {"total": len(names), "intersect": len(cells)},
              "heavy-info": heavy_info,
              "light-info": light_info,
              "cells": {}}
    for i, (hkey, lkey) in enumerate(cells.keys()):
        result["cells"]["G%i" % (i + 1)] = {"heavy": hkey, "light": lkey,
                                            "common": cells[hkey, lkey],
                                            "hcons": str(hcons[hkey].seq) if hkey != common.TRASH_NAME else common.TRASH_NAME,
                                            "lcons": str(lcons[lkey].seq) if lkey != common.TRASH_NAME else common.TRASH_NAME}
    return result


def run(heavy_dir, light_dir, outjson, fix):
    with open(outjson, "wt") as fd:
        fd.write(json.dumps(process(heavy_dir, light_dir, fix)))