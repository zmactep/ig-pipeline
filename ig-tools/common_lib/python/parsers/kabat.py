""" Parser of KABAT file format """

import re


def parse(file):
    """
    Parses kabat file into a {read id : [(region start, region end)]} dictionary, making resulting region start and end
    indices zero-based. Throws exception if some line does not conform to file format.
    """

    line_pattern = re.compile("^[\s]*[^\s]+(([\s]+[\d]+){2})+[\s]*$")
    # possible whitespace at the start, followed by record_id of anything-non-whitespace,
    # followed by pairs of integer numbers, possibly terminated with whitespace at line end

    def fail():
        raise ParseException("Error parsing kabat file")

    def check_line(line):
        return True if line_pattern.match(line) else fail()

    try:
        return {pieces[0]: [(int(start) - 1, int(end) - 1) for start, end in zip(pieces[1::2], pieces[2::2])]
                for pieces in [line.split() for line in file if check_line(line)]}
    except ValueError:
        fail()


class ParseException(Exception):
    pass



