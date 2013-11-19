from itertools import islice


def get_kmers(sequence, k):
    if k < 1 or k > len(sequence): raise ValueError
    return map(''.join, zip(*[islice(sequence, start, None) for start in range(k)])) if 0 < k <= len(sequence) else []