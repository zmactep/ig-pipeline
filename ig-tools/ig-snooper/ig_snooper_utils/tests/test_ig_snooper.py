__author__ = 'Kos'

import argparse
import json
import sys
import os
from pipeline import train_pipeline
from pipeline import predict_pipeline


def run_single_test(train_fasta, train_kabat, test_fasta, test_kabat, config_path):
    print('Running test for (%s, %s, %s, %s)' % (train_fasta, train_kabat, test_fasta, test_kabat))

    params = json.load(open(config_path))
    params.update({'fasta': train_fasta, 'kabat': train_kabat})
    train_pipeline(params)

    params.update({'fasta': test_fasta, 'kabat': test_kabat})
    predict_pipeline(params)


def run_tests(data_dir, config_path):
    for test_dir in filter(lambda x: x.startswith('test-'), os.listdir(data_dir)):
        test_dir_abs = os.path.join(data_dir, test_dir)
        test_reference_data_dir = os.path.join(test_dir_abs, 'test-data')
        test_fasta = os.path.join(test_reference_data_dir, 'vh-test.fasta')
        test_kabat = os.path.join(test_reference_data_dir, 'vh-test.kabat')

        for train_sub_dir in filter(lambda x: x.startswith('train-'), os.listdir(test_dir_abs)):
            train_sub_dir_abs = os.path.join(test_dir_abs, train_sub_dir)
            print('Running tests from %s' % train_sub_dir_abs)
            train_fasta = os.path.join(train_sub_dir_abs, 'vh-train.fasta')
            train_kabat = os.path.join(train_sub_dir_abs, 'vh-train.kabat')

            run_single_test(train_fasta, train_kabat, test_fasta, test_kabat, config_path)


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
