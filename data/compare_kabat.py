#!/usr/bin/python

import argparse
from itertools import zip_longest

def main():
    parser = argparse.ArgumentParser(description='Calculateprediction score for predicted kabat.')
    parser.add_argument('--ref', nargs=1, help='reference kabat')
    parser.add_argument('--input', nargs=1, help='test kabat')
    args = parser.parse_args()

    with open(args.ref.pop(0), 'rU') as reference_file:
        ref_dict = dict([parse_kabat_line(line) for line in reference_file])

    with open(args.input.pop(0), 'rU') as input_file:
        input_dict = dict([parse_kabat_line(line) for line in input_file])

    common_keys = set(ref_dict.keys()) & set(input_dict.keys())
    abs_score = sum([get_abs_score(ref_dict[key], input_dict[key]) for key in common_keys])
    print('Absolute score: %s' % abs_score)

# return: key = sequence_name, value = list of regions len
def parse_kabat_line(line):
    tokens = line.split()
    key = tokens.pop(0).strip()
    value = [int(tokens[i + 1]) - int(tokens[i]) + 1 for i in range(0, len(tokens), 2)]
    return key, value

def get_abs_score(first, second):
    return sum([abs(fst - snd) for fst, snd in zip_longest(first, second, fillvalue=0)])

if __name__ == "__main__":
   main()