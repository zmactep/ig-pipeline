import os


# '*' is a cap marker near the sequence edge. Reads are extended with '*' on the left and on the right sides
# to deal with sliding window near read's edges
nucleo_bases = ['A', 'C', 'G', 'T', 'N', '*']
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


def fix_header(input_file, output_file, output_dir):
    result = open(os.path.join(output_dir, output_file), 'w')
    with open(input_file, 'rU') as input_file:
        for line in input_file:
            if line.startswith('@attribute'):
                result.write(fix_line(line) + '\n')
            else:
                result.write(line)

    input_file.close()
    result.close()


def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


def fix_line(line):
    if line.startswith('@attribute class'):
        return '@attribute class {0,1,2,3,4,5,6}'
    else:
        if len(find_between(line, '{', '}').split(',')) <= len(nucleo_bases):
            return line[:line.index('{') + 1] + ','.join([str(ord(x)) for x in nucleo_bases]) + '}'
        else:
            return line[:line.index('{') + 1] + ','.join([str(ord(x)) for x in amino_acids]) + '}'


