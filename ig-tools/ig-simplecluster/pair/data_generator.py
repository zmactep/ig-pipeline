__author__ = 'mactep'

import common


def get_cell(hcl, lcl):
    return list(set(hcl) & set(lcl))


def make_cells(heavy, light):
    cells = {}
    names = set()
    for hkey in heavy.keys():
        if hkey != common.TRASH_NAME:
            for n in heavy[hkey]:
                names.add(n)
        for lkey in light.keys():
            if lkey != common.TRASH_NAME:
                for n in light[lkey]:
                    names.add(n)
            cell = get_cell(heavy[hkey], light[lkey])
            if len(cell) != 0:
                cells[hkey, lkey] = cell
    return cells, names