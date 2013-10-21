#!/usr/bin/python

import os
import re
from itertools import groupby, islice
from collections import Counter
from Bio import SeqIO
from math import floor


def parse(input_fasta, input_prediction, read_names, output_dir, avg_window_size, merge_threshold):
    debug_prediction = open(os.path.join(output_dir, 'debug_prediction.txt'), 'w')
    debug_prediction_avg = open(os.path.join(output_dir, 'debug_prediction_avg.txt'), 'w')
    results_kabat = open(os.path.join(output_dir, 'results.kabat'), 'w')
    results_pic = open(os.path.join(output_dir, 'results_pic.txt'), 'w')
    input_fasta_iterator = SeqIO.parse(input_fasta, "fasta")

    with open(input_prediction) as pred, open(read_names) as names:
        #name2prediction is a list of key value pairs: key is seq_name, value is region prediction for single nucleotide
        #Nucleotide predictions are sorted by position and listed as one position per line in file (per tuple in list).
        #We need to concatenate predictions for each read.
        name2prediction = []
        for i in range(5):  # skip header
            pred.readline()

        regex = re.compile(":(\d+)")
        for name, val in zip(names, pred):
            val = regex.findall(val)[1]
            name2prediction += [(name.rstrip(), val)]

        # name2prediction is supposed to be sorted here -> no need to sort in groupby
        for key, group in groupby(name2prediction, lambda x: x[0]):
            if not key:
                continue
            # key is a seq_name, and prediction is a list of region predictions for each position in sequence
            prediction = [x[1] for x in group]
            # Post processing step: averaging, merging small regions, etc...
            averaged_prediction = average_prediction(prediction, avg_window_size) # list
            debug_prediction.write("%s\n" % (key + '\t' + ''.join(prediction)))
            debug_prediction_avg.write("%s\n" % (key + '\t' + ''.join(averaged_prediction)))
            regions = kabat_range(averaged_prediction, merge_threshold)
            results_kabat.write("%s\n" % (key + '\t' + '\t'.join([str(x) + '\t' + str(y) for x, y in regions])))
            results_pic.write("%s" % print_regions(input_fasta_iterator.__next__(), regions))

    debug_prediction.close()
    debug_prediction_avg.close()
    results_kabat.close()
    results_pic.close()


def average_prediction(prediction_list, window_size):
    return [str(Counter(window).most_common()[0][0]) for window in slide_window(prediction_list, window_size)]


def slide_window(iterable, window_size):
    shiftedStarts = [islice(iterable, s, None) for s in range(window_size)]
    return zip(*shiftedStarts)


def kabat_range(prediction, merge_threshold):
    result = []
    previous_stop = 0
    first_region_number = prediction[0]
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
    # fill leading missing regions with (0, 0).
    # Ex: if prediction starts with region 2, we add (0,0) twice - for region 0 and 1
    for i in range(int(first_region_number)):
        merged_ranges.insert(0, (0, 0))
    return merged_ranges


def print_regions(seq_record, regions):
    names = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    min_region_len = 6  # minimal printable region len: <CDR1>
    result = []
    region_num = 0
    unknown_region_num = 0
    for (start, stop) in regions:
        region_len = stop - start + 1
        if region_len < min_region_len:
            result += [('*' if unknown_region_num % 2 == 0 else '+') * region_len]
            unknown_region_num += 1
        else:
            region_name = names[region_num] if region_num < len(names) else 'N/A'
            num_of_dashes = region_len - len(region_name) - len('<>')
            result += ['<', '-' * floor(num_of_dashes / 2), region_name,
                       '-' * (num_of_dashes - floor(num_of_dashes / 2)), '>']
            region_num += 1

    return seq_record.id + '\n' + seq_record.seq + '\n' + ''.join(result) + '\n'