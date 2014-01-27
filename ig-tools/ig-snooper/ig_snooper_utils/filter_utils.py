__author__ = 'pavel'

from itertools import islice, chain
from collections import Counter


def slide(seq, n):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def most_frequent(seq):
    return Counter(seq).most_common(1)[0][0]


def sliding_average(seq, size):
    return seq[:int(size / 2)] + "".join(most_frequent(window) for window in slide(seq, size)) + seq[-int(size / 2):]


def indicies_of(seq, char):
    return list(zip(*filter(lambda x: x[1] == char, enumerate(seq))))[0]


def get_centers(seq):
    centers = {}
    for i in set(seq):
        arr = indicies_of(seq, i)
        centers[i] = arr[int(len(arr) / 2)]
    return centers


def get_filled_line(length, symbol):
    return "".join(symbol for _ in range(length))


def mask_except(remains, x):
    return x if x in remains else '_'


def filter_between_fr(centers, seq):
    result = []
    for window in slide(sorted(centers.keys()), 2):
        result.append(list(map(lambda x: mask_except([window[0], window[1]], x), seq[centers[window[0]]:centers[window[-1]]])))
    return list(chain(*result))


def fill_gaps(seq):
    buffer = ['_' for _ in range(len(seq))]
    borders = list(chain(*[[seq.find(str(i)), seq.rfind(str(i)) + 1] for i in range(0, 7, 2)]))
    if list(sorted(borders)) != borders:
        return
    for index, border in enumerate(slide(borders, 2)):
        for i in range(border[0], border[-1]):
            buffer[i] = str(index)
    return "".join(buffer)


def fix(input_str, start_window=3, end_window=7, step_window=2):
    for size_window in range(start_window, end_window + 1, step_window):
        s = sliding_average(input_str, size_window)
        centers = get_centers(s)
        frs = {i: centers[i] for i in filter(lambda x: int(x) % 2 == 0, centers)}
        frs_only_string = get_filled_line(centers['0'], '0') + "".join(filter_between_fr(frs, s)) + get_filled_line(len(s) - centers['6'], '6')
        fg = fill_gaps(frs_only_string)
        if fg:
            return fg
    return ""