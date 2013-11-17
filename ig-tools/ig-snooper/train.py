import argparse
import json
import sys
import os
from subprocess import Popen, PIPE
import datetime
from ig_snooper_utils import fix_weka_header
import svm_data_generator as sdg

svm_data_generator = os.path.join(sdg.__path__[0], 'bin/svm_data_generator')
weka = 'common_lib/third_party/weka-3.6.10/weka.jar'


def main():
    try:
        start_time = datetime.datetime.now()
        params = get_params()
        check_files(params)

        # prepare data: convert .fasta to .libsvm
        print('Generating train data in libsvm format...')
        output, error = Popen([os.path.join(params['tools_root'], svm_data_generator), 'train', params['fasta'],
                               params['kabat'], str(params['ml_window_size']), '1', params['outdir']],
                              stdout=PIPE, stderr=PIPE).communicate()

        train_libsvm_path = os.path.join(params['outdir'], 'train.libsvm')
        if not (os.path.exists(train_libsvm_path) and os.path.getsize(train_libsvm_path) > 0):
            raise FileNotFoundError('Failed to convert .fasta to .libsvm: %s; %s' % (output, error))

        # run weka conversion
        print('Done. Applying NumericToNominal conversion...')
        os.environ["CLASSPATH"] = os.path.join(params['tools_root'], weka)
        output, error = Popen(['java', '-Xmx4096M', 'weka.filters.unsupervised.attribute.NumericToNominal', '-i',
                               os.path.join(params['outdir'], 'train.libsvm'), '-o',
                               os.path.join(params['outdir'], 'train_nominal.arff')],
                              stdout=PIPE, stderr=PIPE).communicate()

        train_nominal_path = os.path.join(params['outdir'], 'train_nominal.arff')
        if not (os.path.exists(train_nominal_path) and os.path.getsize(train_nominal_path) > 0):
            raise FileNotFoundError('Failed to convert Numeric to Nominal: %s; %s' % (output, error))

        fix_weka_header.fix_header(os.path.join(params['outdir'], 'train_nominal.arff'), 'train_nominal_fixed.arff',
                                   params['outdir'])

        # run weka training
        print('Done. Train...')
        output, error = Popen(['java', '-Xmx4096M', 'weka.classifiers.trees.RandomForest', '-I', '10',
                               '-K', '0', '-S', '1', '-no-cv', '-p', '0',
                               '-t', os.path.join(params['outdir'], 'train_nominal_fixed.arff'),
                               '-d', os.path.join(params['outdir'], 'model.model')],
                              stdout=PIPE, stderr=PIPE).communicate()

        model_path = os.path.join(params['outdir'], 'model.model')
        if not (os.path.exists(model_path) and os.path.getsize(model_path) > 0):
            raise FileNotFoundError('Failed to train: %s; %s' % (output, error))

        # final cleanup if necessary
        if 'clean_up' in params:
            print('Deleting unnecessary files...')
            for f in ['train.libsvm', 'train_nominal.arff', 'train_nominal_fixed.arff']:
                os.remove(os.path.join(params['outdir'], f))

        print('Done in %s' % str(datetime.datetime.now() - start_time))
        print('Your model is in %s' % model_path)

    except Exception as detail:
        print('%s' % detail)
        sys.exit()


def get_params():
    parser = argparse.ArgumentParser(description='Train model given train dataset')
    parser.add_argument('--config_path', nargs='?', help='Path to config')
    parser.add_argument('--fasta', help='Train data in fasta')
    parser.add_argument('--kabat', help='Train data in kabat')
    parser.add_argument('--ml_window_size', nargs='?', help='ML windows size')
    parser.add_argument('--tools_root', nargs='?', help='Tools dir')
    parser.add_argument('--outdir', nargs='?', help='Output directory')

    params = {k: v for k, v in vars(parser.parse_args()).items() if v is not None}

    if 'config_path' in params and not params['config_path'] is None:
        config = json.load(open(params['config_path']))
        config.update(params)  # commandline args has higher priority
        params = config

    required = ('fasta', 'kabat', 'ml_window_size', 'tools_root', 'outdir')
    if not all(key in params and params[key] is not None for key in required):
        raise KeyError('Missing required params. Make sure that you have all of these: %s' % ' '.join(required))

    return params


def check_files(params):
    svm_data_generator_path = os.path.join(params['tools_root'], svm_data_generator)
    weka_path = os.path.join(params['tools_root'], weka)

    if not os.path.exists(svm_data_generator_path):
        raise FileNotFoundError('svm_data_generator not found in %s' % svm_data_generator_path)

    if not os.path.exists(weka_path):
        raise FileNotFoundError('weka.jar not found in %s' % weka_path)

    if not os.path.exists(params['fasta']):
        raise FileNotFoundError('fasta file not found in %s' % params['fasta'])

    if not os.path.exists(params['kabat']):
        raise FileNotFoundError('kabat not found in %s' % params['kabat'])

    if not os.path.exists(params['outdir']):
        raise FileNotFoundError('output dir not found in %s' % params['outdir'])


if __name__ == "__main__":
   main()