__author__ = 'Kos'

import argparse
import shutil
import multiprocessing
import json
import sys
import os
from pipeline import train_pipeline
from pipeline import predict_pipeline
from diff_info import get_diff
from compare_marking import compare_marking


def run_single_prediction(working_dir, train_fasta, train_kabat, test_fasta, test_kabat, config_path):
    print('Running test for (%s, %s, %s, %s)' % (train_fasta, train_kabat, test_fasta, test_kabat))

    params = json.load(open(config_path))
    params.update({'fasta': train_fasta, 'kabat': train_kabat})
    params['model_path'] = working_dir
    params['outdir'] = working_dir
    train_pipeline(params)

    params.update({'fasta': test_fasta, 'kabat': test_kabat})
    predict_pipeline(params)


def run_diff_info(result_kabat, test_kabat):
    d, error_list = get_diff(result_kabat, test_kabat, False, 0)
    metrics_line = []
    with open(os.path.join(os.path.dirname(result_kabat), 'diff_info.txt'), 'a') as f:
        print("REGIONS:\terrors\t(signed\t/ unsigned)\t(signed\t/ unsigned)\t error rate", file=f)
        for name, (v, s, u, se, ue, er) in d:
            print("%s:\t%d\t(%.3f\t/ %.3f)\t(%.3f\t/ %.3f)\t%.4f" % (name, v, s, u, se, ue, er), file=f)
            metrics_line += [v, s, u, se, ue, er]
        print("Double-sequence list:", file=f)
        for e in error_list:
            print(e, file=f)
    return metrics_line


def run_compare_metrics(result_kabat, test_kabat):
    total_error, seq_error,  misfit_error, total = compare_marking(result_kabat, test_kabat, False, 0)
    d = {
        'wrong_rate': (total_error / total) if total else 0,
        'error_rate': (seq_error / total) if total else 0,
        'misfit_rate': (misfit_error / total) if total else 0,
        'error_wrong_rate': (seq_error / total_error) if total_error else 0,
        'misfit_wrong_rate': (misfit_error / total_error) if total_error else 0
    }
    with open(os.path.join(os.path.dirname(result_kabat), 'compare_info.txt'), 'a') as f:
        print("Wrong sequence annotation rate: %.4f" % d['wrong_rate'], file=f)
        print("Average error rate per sequence: %.4f" % d['error_rate'], file=f)
        print("Average region misfit per sequence: %.4f" % d['misfit_rate'], file=f)
        print("Average error rate per wrong sequence: %.4f" % d['error_wrong_rate'], file=f)
        print("Average region misfit per wrong sequence: %.4f" % d['misfit_wrong_rate'], file=f)
    return d


def get_data_paths(test_reference_data_dir, train_sub_dir_abs):
    test_fasta = os.path.join(test_reference_data_dir, 'vh-test.fasta')
    test_kabat = os.path.join(test_reference_data_dir, 'vh-test.kabat')
    train_fasta = os.path.join(train_sub_dir_abs, 'vh-train.fasta')
    train_kabat = os.path.join(train_sub_dir_abs, 'vh-train.kabat')
    return train_fasta, train_kabat, test_fasta, test_kabat


def print_stat_header(test_dir_abs):
    rates = ['test_dir', 'wrong_rate', 'error_rate', 'misfit_rate', 'error_wrong_rate', 'misfit_wrong_rate']
    regions = ["%s%d %s %s" % (r, n, se, err) for n in range(1, 5) for r in ["FR", "CDR"] for se in ["start", "end"]
               for err in ["errors", "signed", "unsigned", "signed", "unsigned", "error rate"]]
    header = rates + regions

    with open(os.path.join(test_dir_abs, 'stat.txt'), 'a') as f:
        f.write('\t'.join(header) + '\n')


def dump_stat(test_dir_abs, train_sub_dir, metrics_diff, metrics_compare):
    data = [train_sub_dir, metrics_compare['wrong_rate'], metrics_compare['error_rate'],
            metrics_compare['misfit_rate'], metrics_compare['error_wrong_rate'],
            metrics_compare['misfit_wrong_rate']]
    data += metrics_diff

    with open(os.path.join(test_dir_abs, 'stat.txt'), 'a') as f:
        f.write(data[0] + '\t' + '\t'.join(map(lambda x: '%.4f' % x, data[1:])) + '\n')


def run_job(data):
    print('Job params: %s %s %s %s' % data)
    index, test_dir, data_dir, config_path = data
    test_dir_abs = os.path.join(data_dir, test_dir)
    test_reference_data_dir = os.path.join(test_dir_abs, 'test-data')

    print_stat_header(test_dir_abs)

    working_dir_base = os.path.join('/tmp', test_dir)
    if not os.path.exists(working_dir_base):
        os.mkdir(working_dir_base)

    for train_sub_dir in filter(lambda x: x.startswith('train-'), os.listdir(test_dir_abs)):
        train_sub_dir_abs = os.path.join(test_dir_abs, train_sub_dir)
        print('Running tests from %s' % train_sub_dir_abs)

        train_fasta, train_kabat, test_fasta, test_kabat = get_data_paths(test_reference_data_dir, train_sub_dir_abs)
        try:
            working_dir = os.path.join(working_dir_base, train_sub_dir)
            if not os.path.exists(working_dir):
                os.mkdir(working_dir)

            run_single_prediction(working_dir, train_fasta, train_kabat, test_fasta, test_kabat, config_path)
            result_kabat = os.path.join(train_sub_dir_abs, 'prediction.kabat')
            shutil.copyfile(os.path.join(working_dir, 'results.kabat'), result_kabat)
            dump_stat(test_dir_abs, train_sub_dir, run_diff_info(result_kabat, test_kabat), run_compare_metrics(result_kabat, test_kabat))
        except:
            (type, value, traceback) = sys.exc_info()
            print("Unexpected error: ", value)


def run_tests(data_dir, config_path):
    test_dirs = [d for d in filter(lambda x: x.startswith('test-'), os.listdir(data_dir))]
    jobs = [(index, job, data_dir, config_path) for (index, job) in enumerate(test_dirs)]

    po = multiprocessing.Pool()
    po.map(run_job, jobs)
    po.close()
    po.join()


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('data_dir', help='path to dir with ig-pipeline/data/test/germline like structure')
    parser.add_argument('config_path', help='path to ig-snooper-config like ig-pipeline/ig-config/ig-snooper.conf')

    return parser.parse_args()


def main():
    args = parse_args()

    if not (os.path.exists(args.data_dir) and os.path.isdir(args.data_dir)):
        print('%s doesn\'t exist or it is not a directory' % args.data_dir)
        sys.exit()

    if not (os.path.exists(args.config_path) and os.path.isfile(args.config_path)):
        print('%s doesn\'t exist or it is a directory' % args.config_path)
        sys.exit()

    run_tests(args.data_dir, args.config_path)

if __name__ == "__main__":
    main()
