#!/usr/bin/python

import sys
from itertools import groupby, islice
from collections import Counter

def main(argv):
    input_file = argv[1]
    sliding_window_size = int(argv[2])
    svm_window_size = int(argv[3])
    svm_head_offset = svm_window_size // 2
    svm_tail_offset = svm_window_size - svm_head_offset

    debug_prediction = open('debug_prediction.txt', 'w')
    debug_prediction_avg = open('debug_prediction_avg.txt', 'w')
    results = open('results.txt', 'w')
    with open(input_file) as file:
        name2prediction = [line.rstrip().split('\t')[:2] for line in file]
        # name2prediction is supposed to be sorted here -> no need to sort in groupby
        for key, group in groupby(name2prediction, lambda x: x[0]):
            prediction = [x[1] for x in group] # list
            averaged_prediction = average_prediction(prediction, sliding_window_size) # list
            debug_prediction.write("%s\n" % (key + '\t' + ''.join(prediction)))
            debug_prediction_avg.write("%s\n" % (key + '\t' + ''.join(averaged_prediction)))
            results.write("%s\n" % kabat_range(key, averaged_prediction, svm_head_offset, svm_tail_offset))
    debug_prediction.close()
    debug_prediction_avg.close()
    results.close()

def average_prediction(prediction_list, window_size):
    avg = [str(Counter(window).most_common()[0][0]) for window in slide_window(prediction_list, window_size)]
    # restore cut-offed head and tail from prediction_list
    missing_edge_count = (window_size - 1) // 2
    head = avg[0: missing_edge_count]
    tail = avg[-missing_edge_count:]
    return head + avg + tail

def slide_window(iterable, window_size):
    shiftedStarts = [islice(iterable, s, None) for s in range(window_size)]
    return zip(*shiftedStarts)

def kabat_range(name, prediction, svm_head_offset, svm_tail_offset):
    array = []
    offset = 0
    for key, group in groupby(prediction):
        range = sum([1 for x in group])
        start = offset + 1
        stop = start + range - 1
        array += start, stop
        offset = stop
    result = [x + svm_head_offset for x in array]
    result[0] = 1
    result[-1] += svm_tail_offset
    return name + '\t' + '\t'.join([str(x) for x in result])

if __name__ == "__main__":
   main(sys.argv)