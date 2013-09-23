#!/usr/bin/python

import sys

def main(argv):
    input_file = argv[1]
    result = open('result.arff', 'w')
    header_file = open('header.txt', 'w')
    header = []
    is_header_end = False
    with open(input_file, 'rU') as input_file:
        for line in input_file:
            if line.startswith('@') or line.startswith('\n'):
                header += [line]
            else:
                if not is_header_end:
                    fixed_header = get_fixed_header(header)
                    result.writelines(fixed_header)
                    header_file.writelines(fixed_header)
                    is_header_end = True
                result.write(line)
    input_file.close()
    header_file.close()
    result.close()

def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""

def get_fixed_header(header):
    features_set = set()
    num_of_features = 0
    result = header[0:2]
    header.pop(0)
    header.pop(0)
    for line in header:
        if line.startswith('@attribute att_'):
            num_of_features += 1
            features_set = features_set.union({x for x in find_between(line, '{', '}').split(',')})
        if line.startswith('@attribute class'):
            result += ['@attribute att_' + str(i + 1) + ' {' + ','.join(features_set) + '}\n' for i in range(num_of_features)]
            result += [line]
    result += ['\n', '@data\n', '\n']
    return result


if __name__ == "__main__":
   main(sys.argv)
