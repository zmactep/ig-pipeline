import argparse
import json
import sys
import os
from subprocess import Popen, PIPE
import datetime
import fix_weka_header, parse_svm_output, compare_kabat

svm_data_generator = 'ig-snooper/svm_data_generator/bin/svm_data_generator'
weka = 'common_lib/third_party/weka-3.6.10/weka.jar'


def main():
    try:
        start_time = datetime.datetime.now()
        params = get_params()
        check_files(params)

        # prepare data: convert .fasta to .libsvm
        print('Generating predict data in libsvm format...')
        output, error = Popen([os.path.join(params['tools_root'], svm_data_generator), 'predict', params['fasta'],
            str(params['ml_window_size']), params['outdir']], stdout=PIPE, stderr=PIPE).communicate()

        predict_libsvm_path = os.path.join(params['outdir'], 'predict.libsvm')
        if not (os.path.exists(predict_libsvm_path) and os.path.getsize(predict_libsvm_path) > 0):
            raise FileNotFoundError('Failed to convert .fasta to .libsvm: %s; %s' % (output, error))

        # run weka conversion
        print('Done. Applying NumericToNominal conversion...')
        os.environ["CLASSPATH"] = os.path.join(params['tools_root'], weka)
        output, error = Popen(['java', '-Xmx4096M', 'weka.filters.unsupervised.attribute.NumericToNominal', '-i',
            os.path.join(params['outdir'], 'predict.libsvm'), '-o', os.path.join(params['outdir'], 'predict_nominal.arff')],
            stdout=PIPE, stderr=PIPE).communicate()

        predict_nominal_path = os.path.join(params['outdir'], 'predict_nominal.arff')
        if not (os.path.exists(predict_nominal_path) and os.path.getsize(predict_nominal_path) > 0):
            raise FileNotFoundError('Failed to convert Numeric to Nominal: %s; %s' % (output, error))

        fix_weka_header.fix_header(os.path.join(params['outdir'], 'predict_nominal.arff'), 'predict_nominal_fixed.arff',
                                   params['outdir'])

        # run weka training
        print('Done. Predict...')
        prediction_file_name = os.path.join(params['outdir'], 'prediction.txt')
        with open(prediction_file_name, 'w') as out:
            Popen(['java', '-Xmx4096M', 'weka.classifiers.trees.RandomForest', '-no-cv', '-p', '0', '-l',
                params['model_path'], '-T', os.path.join(params['outdir'], 'predict_nominal_fixed.arff')],
                stdout=out).wait()

        if not (os.path.exists(prediction_file_name) and os.path.getsize(prediction_file_name) > 0):
            raise FileNotFoundError('Failed to predict.')

        parse_svm_output.parse(params['fasta'], prediction_file_name, os.path.join(params['outdir'], 'read_names.txt'),
                               params['outdir'], params['avg_window_size'], params['merge_threshold'])

        if 'kabat' in params and os.path.exists(params['kabat']):
            print('Calculating quality metrics...')
            compare_kabat.compare(params['kabat'], os.path.join(params['outdir'], 'results.kabat'), params['outdir'])


        # final cleanup if necessary
        if params['clean_up']:
            print('Deleting unnecessary files...')
            for f in ['predict.libsvm', 'predict_nominal.arff', 'predict_nominal_fixed.arff', 'read_names.txt',
                      'prediction.txt', 'debug_prediction.txt', 'debug_prediction_avg.txt']:
                os.remove(os.path.join(params['outdir'], f))

        print('Done in %s' % str(datetime.datetime.now() - start_time))

    except (KeyError, FileNotFoundError) as detail:
        print('%s' % detail)
        sys.exit()


def get_params():
    parser = argparse.ArgumentParser(description='Train model given train dataset')
    parser.add_argument('--config_path', nargs='?', help='Path to config')
    parser.add_argument('--fasta', help='Predict data in fasta')
    parser.add_argument('--kabat', help='Optional kabat for benchmarking')
    parser.add_argument('--model_path', help='input model path')
    parser.add_argument('--ml_window_size', nargs='?', help='ML windows size')
    parser.add_argument('--avg_window_size', nargs='?', help='averaging windows size')
    parser.add_argument('--merge_threshold', nargs='?', help='merge small regions if len < merge_threshold')
    parser.add_argument('--tools_root', nargs='?', help='Tools dir')
    parser.add_argument('--outdir', nargs='?', help='Output directory')

    params = {k: v for k, v in vars(parser.parse_args()).items() if v is not None}

    if 'config_path' in params and not params['config_path'] is None:
        config = json.load(open(params['config_path']))
        config.update(params) # commandline args has higher priority
        params = config

    required = ('fasta', 'model_path', 'ml_window_size', 'avg_window_size', 'merge_threshold', 'tools_root', 'outdir')
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

    if not os.path.exists(params['model_path']):
        raise FileNotFoundError('model not found in %s' % params['model_path'])

    if not os.path.exists(params['outdir']):
        raise FileNotFoundError('output dir not found in %s' % params['outdir'])


if __name__ == "__main__":
   main()