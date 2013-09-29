#!/usr/bin/python

import argparse
from itertools import groupby, islice
from collections import Counter

def main():
    parser = argparse.ArgumentParser(description='Calculateprediction score for predicted kabat.')
    parser.add_argument('--input_file', nargs=1, help='input kabat')
    parser.add_argument('--sliding_window_size', nargs=1, type=int, help='sliding window size')
    parser.add_argument('--merge_threshold', nargs=1, type=int, help='threshold len for merge with neighbour')
    parser.add_argument('--output', nargs=1, help='output_dir')
    args = parser.parse_args()

    input_file = args.input_file.pop(0)
    output_dir = args.output.pop(0)
    sliding_window_size = args.sliding_window_size.pop(0)
    merge_threshold = args.merge_threshold.pop(0)

    debug_prediction = open(output_dir + 'debug_prediction.txt', 'w')
    debug_prediction_avg = open(output_dir + 'debug_prediction_avg.txt', 'w')
    results = open(output_dir + 'results.txt', 'w')

    with open(input_file) as file:
        #name2prediction is a list of key value pairs: key is seq_name, value is region prediction for single nucleotide
        #Nucleotide predictions are sorted by position and listed as one position per line in file (per tuple in list).
        #We need to concatenate predictions for each read.
        name2prediction = [line.rstrip().split('\t')[:2] for line in file]
        # name2prediction is supposed to be sorted here -> no need to sort in groupby
        for key, group in groupby(name2prediction, lambda x: x[0]):
            if not key:
                continue
            # key is a seq_name, and prediction is a list of region predictions for each position in sequence
            prediction = [x[1] for x in group]
            # Post processing step: averaging, merging small regions, etc...
            averaged_prediction = average_prediction(prediction, sliding_window_size) # list
            debug_prediction.write("%s\n" % (key + '\t' + ''.join(prediction)))
            debug_prediction_avg.write("%s\n" % (key + '\t' + ''.join(averaged_prediction)))
            results.write("%s\n" % kabat_range(key, averaged_prediction, merge_threshold))

    debug_prediction.close()
    debug_prediction_avg.close()
    results.close()

def average_prediction(prediction_list, window_size):
    return [str(Counter(window).most_common()[0][0]) for window in slide_window(prediction_list, window_size)]

def slide_window(iterable, window_size):
    shiftedStarts = [islice(iterable, s, None) for s in range(window_size)]
    return zip(*shiftedStarts)

def kabat_range(name, prediction, merge_threshold):
    result = []
    previous_stop = 0
    for key, group in groupby(prediction):
        sum_range = sum([1 for x in group])
        start = previous_stop + 1
        stop = start + sum_range - 1
        result += start, stop
        previous_stop = stop

    #merge regions with len less than threshold
    #ex merge_threshold = 7:
    #[(1, 2), (3, 5), (6, 20)] -> [(1, 20)]
    #[(10, 20), (30, 50)] -> [(10, 20), (30, 50)]
    previous = (result[0], result[1])
    result.pop(0)
    result.pop(0)
    merged_ranges = []
    for x, y in [(result[i], result[i + 1]) for i in range(0, len(result), 2)]:
        if previous[1] - previous[0] + 1 < merge_threshold:
            # region is small, merge with previous one and wait until merged
            # region len will be sufficient to add to result (in else clause)
            previous = previous[0], y
        else:
            merged_ranges.append(previous)
            previous = (x, y)
    merged_ranges.append(previous)
    return name + '\t' + '\t'.join([str(x) + '\t' + str(y) for x, y in merged_ranges])

if __name__ == "__main__":
   main()