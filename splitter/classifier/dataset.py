class DomainTypeDataset(object):
    def __init__(self):
        self.data = {}

    def setData(self, d_type, data):
        self.data[d_type] = data

    def getData(self, k, window_mode=False):
        rin, rout = [], []
        for d_type in self.data:
            for seq in self.data[d_type]:
                for i in range(len(seq)):
                    if not window_mode:
                        rin.append(encode(create_region(seq, i, k)))
                    else:
                        if i > len(seq) - k:
                            break
                        rin.append(encode(seq[i:i+k]))
                    rout.append(d_type)
        return rin, rout


def encode(ns):
    nucleo = "ACGT"
    result = []
    for c in ns:
        if c == 'N':
            l = [1] * len(nucleo)
        else:
            l = [int(i) for i in reversed(bin(2 ** nucleo.index(c))[2:])]
            l = l + [0] * (len(nucleo) - len(l))
        for i in reversed(l):
            result.append(i)
    return result


def create_region(ns, pos, k):
    return 'N' * (k - pos) + ns[pos - k if pos - k > 0 else 0: pos + k + 1] +\
           'N' * (pos + k + 1 - len(ns))
