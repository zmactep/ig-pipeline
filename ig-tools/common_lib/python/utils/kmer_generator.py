from itertools import islice


def get_sequence_kmers(sequence, k):
    """ Generates k-mers from any iterable sequence as a series of its subsequences """
    if k < 1 or k > len(sequence): raise ValueError("Window length out of bounds: {0}".format(k))
    return zip(*[islice(sequence, start, None) for start in range(k)])


def get_string_kmers(charsequence, k):
    """ Generates k-mers of a character sequence (a string or iterable collection of chars) """
    return map(''.join, get_sequence_kmers(charsequence, k))
