import copy
from functools import reduce
from itertools import starmap, product
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATA_DIR = 'data'

class Sequence:
    def __init__(self, name, sequence, nomenclature):
        self.name = name
        self.sequence = sequence
        self.nomenclature = nomenclature

    def __str__(self):
        return '{0}\n{1}'.format(self.sequence, ''.join([str(i)*(reg[1] - reg[0] + 1) for i, reg in self.nomenclature.items()]))

    def merge(self, other):
        result = copy.deepcopy(self)
        my_last_region = max(self.nomenclature)
        others_first_unmerged_region = min(other.nomenclature)
        if my_last_region == others_first_unmerged_region:
            my_tuple = self.nomenclature[my_last_region]
            result.nomenclature[my_last_region] = (my_tuple[0], my_tuple[1] + other.nomenclature[others_first_unmerged_region][1])
            others_first_unmerged_region += 1
        for r in range(others_first_unmerged_region, max(other.nomenclature) + 1):
            if r in other.nomenclature:
                assert r not in self.nomenclature
                others_tuple = other.nomenclature[r]
                result.nomenclature[r] = (len(self.sequence) + others_tuple[0], len(self.sequence) + others_tuple[1])
        result.sequence += other.sequence
        result.name += '_' + other.name
        return result



class SegmentParser:
    REGION_ID = {'FR1': 0, 'CDR1': 1, 'FR2': 2, 'CDR2': 3, 'FR3': 4, 'CDR3': 5, 'FR4': 6}

    def __init__(self, path):
        self.sequences = []
        self._parse(path)

    def __iter__(self):
        return iter(self.sequences)

    def _parse(self, path):
        with open(path) as file:
            prefix = file.readline().strip()
            region_names = file.readline().split()

            for line in file:
                name = prefix + line.split()[0]
                region_sequences = line.split()[1:]

                seq = ''
                nomenclature = {}
                for region_seq, region_name in zip(region_sequences, region_names):
                    region_seq = region_seq.replace('-', '')
                    start_index = len(seq) + 1
                    seq += region_seq
                    nomenclature[self.REGION_ID[region_name]] = (start_index, len(seq))

                self.sequences.append(Sequence(name, seq, nomenclature))


def generate_sequences(files):
    for seqs in product(*[SegmentParser(f) for f in files]):
        yield reduce(Sequence.merge, seqs)


def write_data(fasta_filename, kabat_filename, *inputs):
    with open(fasta_filename, 'w') as fasta, open(kabat_filename, 'w') as kabat:
        for s in generate_sequences(inputs):
            #fasta.write('>{0}\n{1}\n'.format(s.name, str(s)))
            SeqIO.write(SeqRecord(Seq(s.sequence), id=s.name, description=''), fasta, 'fasta')
            kabat.write('{0}\t{1}\n'.format(s.name, '\t'.join(
                str(d) for r in sorted(s.nomenclature) for d in s.nomenclature[r])))


def main():
    write_data('vh.fasta', 'vh.kabat', *(os.path.join(DATA_DIR, f) for f in ['vh.txt', 'dh.txt', 'jh.txt']))
    write_data('vl.fasta', 'vl.kabat', *(os.path.join(DATA_DIR, f) for f in ['vl.txt', 'jl.txt']))
    write_data('vk.fasta', 'vk.kabat', *(os.path.join(DATA_DIR, f) for f in ['vk.txt', 'jk.txt']))

if __name__ == '__main__':
    main()