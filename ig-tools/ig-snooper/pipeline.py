"""
This module holds boilerplate code for construction of a generic pipeline of actions, and a complete pipeline for
running training and prediction of antibody regions using WEKA. It serves as a backend for two client-facing scripts,
'train.py' and 'predict.py'.
"""
from datetime import datetime
import logging
import os
from subprocess import Popen, PIPE
from ig_snooper_utils import svm_data_generator, fix_weka_header, parse_svm_output, compare_kabat


WEKA_PATH = 'common_lib/third_party/weka-3.6.10/weka.jar'


PREDICTOR_FIXED_ARGS = {'RandomForest': {'train': ['weka.classifiers.trees.RandomForest', '-I', '10', '-K', '0', '-S', '1'],
                                         'predict': ['weka.classifiers.trees.RandomForest']}}


# The following functions produce pipeline stages as closures:


def _get_check_files_action(pathlist):
    def action():
        for path in pathlist:
            if not os.path.exists(path):
                logging.error('Path not exists: %s' % path)
                raise FileNotFoundError(path)

    return action


def _get_convert_to_libsvm_action(args, generate_comments=False):
    def action():
        try:
            svm_data_generator.process(args['ml_window_size'], True, args['outdir'], 'dataset.libsvm',
                                       args['fasta'], args['kabat'] if 'kabat' in args else None,
                                       'read_names.txt' if generate_comments else None)
        except Exception:
            logging.error('Failed to convert .fasta to .libsvm, converter reported error')
            raise

        libsvm_path = os.path.join(args['outdir'], 'dataset.libsvm')
        if not (os.path.exists(libsvm_path) and os.path.getsize(libsvm_path) > 0):
            logging.error('Failed to convert .fasta to .libsvm, file is missing or empty')
            raise RuntimeError()

    return action


def _get_weka_conversion_action(args):
    def action():
        os.environ['CLASSPATH'] = os.path.join(args['tools_root'], WEKA_PATH)
        try:
            Popen(['java', '-Xmx4096M', 'weka.filters.unsupervised.attribute.NumericToNominal', '-i',
                   os.path.join(args['outdir'], 'dataset.libsvm'), '-o',
                   os.path.join(args['outdir'], 'dataset_nominal.arff')],
                  stdout=PIPE, stderr=PIPE).communicate()
        except Exception:
            logging.error('Failed to convert numeric to nominal: exception in calling weka converter')
            raise

        nominal_path = os.path.join(args['outdir'], 'dataset_nominal.arff')
        if not (os.path.exists(nominal_path) and os.path.getsize(nominal_path) > 0):
            logging.error('Failed to convert numeric to nominal: output file missing or empty')
            raise RuntimeError()

        try:
            fix_weka_header.fix_header(os.path.join(args['outdir'], 'dataset_nominal.arff'),
                                       'dataset_nominal_fixed.arff', args['outdir'])
        except Exception:
            logging.error('Error in fixing .arff header')
            raise

    return action


def _get_weka_train_action(args):
    def action():
        try:
            command = ['java', '-Xmx4096M'] + PREDICTOR_FIXED_ARGS[args['predictor']]['train'] + \
                      ['-no-cv', '-p', '0', '-t', os.path.join(args['outdir'], 'dataset_nominal_fixed.arff'),
                       '-d', os.path.join(args['outdir'], 'model.model')]
            logging.debug('Train command is: %s' % '\t'.join(command))
            Popen(command, stdout=PIPE, stderr=PIPE).communicate()
        except Exception:
            logging.error('Error in calling weka training')
            raise

        model_path = os.path.join(args['outdir'], 'model.model')

        if not (os.path.exists(model_path) and os.path.getsize(model_path) > 0):
            logging.error('Failed to train: output file missing or empty')
            raise RuntimeError()

    return action


def _get_weka_predict_action(args):
    def action():
        prediction_file_name = os.path.join(args['outdir'], 'prediction.txt')

        try:
            with open(prediction_file_name, 'w') as out:
                command = ['java', '-Xmx4096M'] + PREDICTOR_FIXED_ARGS[args['predictor']]['predict'] + \
                          ['-no-cv', '-p', '0', '-l', os.path.join(args['model_path'], 'model.model'), '-T',
                           os.path.join(args['outdir'], 'dataset_nominal_fixed.arff')]
                logging.debug('Predict command is: %s' % '\t'.join(command))
                Popen(command, stdout=out).wait()
        except Exception:
            logging.error('Error in running weka prediction')
            raise

        if not (os.path.exists(prediction_file_name) and os.path.getsize(prediction_file_name) > 0):
            logging.error('Failed to predict: output file %s is missing or empty' % prediction_file_name)
            raise RuntimeError()

        try:
            parse_svm_output.parse(args['fasta'], prediction_file_name, os.path.join(args['outdir'], 'read_names.txt'),
                                   args['outdir'], int(args['avg_window_size']), int(args['merge_threshold']))
        except Exception:
            logging.error('Error in parsing prediction output')
            raise

    return action


def _get_kabat_quality_metrics_action(args):
    def action():
        if not os.path.exists(args['kabat']):
            logging.error('File does not exist: %s' % args['kabat'])
            raise FileNotFoundError()
        try:
            compare_kabat.compare(args['kabat'], os.path.join(args['outdir'], 'results.kabat'), args['outdir'])
        except Exception:
            logging.error('Error in comparing resulting .kabat')
            raise RuntimeError()

    return action


def _get_cleanup_action(pathlist):
    def action():
        for f in pathlist:
            try:
                os.remove(f)
            except Exception:
                logging.error('Error removing file {0}'.format(f.name))
                raise

    return action


class _Stage():
    """
    Represents a stage in processing pipeline. Holds 'desc', a message printed before execution, and 'action', a closure
    doing the job, returning a tuple of (bool, string): the success flag and outcome message. On their meaning, read
        the _pipeline method description.
    """

    def __init__(self, desc, action):
        self._desc = desc
        self._action = action

    def run(self):
        logging.info(self._desc)
        try:
            self._action()
        except Exception:
            logging.error('"{0}" stage failed.'.format(self._desc))
            raise


def _pipeline(tasks):
    """
    Executes a series of tasks by iterating over provided collection of _Task objects. For every task, first prints its
    description, then calls the 'run' method. It is expected to return a tuple of (success flag, outcome message).
    Processing continues over to the next task only if success flag is true, otherwise the outcome message gets
    printed and execution terminates. Return value is a boolean indicating successful completion of the pipeline overall
    """
    start_time = datetime.now()

    try:
        for task in tasks:
            task.run()
    except Exception:
        logging.error("Pipeline failed")
        raise

    logging.info('Done in %s' % str(datetime.now() - start_time))
    return True


def train_pipeline(args):
    """
    The pre-assembled pipeline for training of antibody region markup using WEKA. Args is a dictionary of command line
    arguments.
    """
    chain = [
        _Stage('Checking input files', _get_check_files_action(
            [os.path.join(args['tools_root'], WEKA_PATH), args['fasta'], args['kabat'], args['outdir']])),
        _Stage('Converting data to libsvm format', _get_convert_to_libsvm_action(args)),
        _Stage('Applying NumericToNominal conversion', _get_weka_conversion_action(args)),
        _Stage('Training', _get_weka_train_action(args))]
    if 'clean_up' in args:
        chain += [_Stage('Deleting temporary files', _get_cleanup_action(
            [os.path.join(args['outdir'], f)
             for f in ['dataset.libsvm', 'dataset_nominal.arff', 'dataset_nominal_fixed.arff']]))]
    return _pipeline(chain)


def predict_pipeline(args):
    """
    The pre-assembled pipeline for prediction of antibody region markup using WEKA.
    """
    chain = [
        _Stage('Checking input files', _get_check_files_action(
            [os.path.join(args['tools_root'], WEKA_PATH), args['fasta'], args['model_path'], args['outdir']])),
        _Stage('Converting data to libsvm format', _get_convert_to_libsvm_action(args, generate_comments=True)),
        _Stage('Applying NumericToNominal conversion', _get_weka_conversion_action(args)),
        _Stage('Predicting', _get_weka_predict_action(args))]
    if 'kabat' in args:
        chain += [_Stage('Evaluating quality metrics', _get_kabat_quality_metrics_action(args))]
    if 'clean_up' in args:
        chain += [_Stage('Deleting unnecessary files', _get_cleanup_action(
            [os.path.join(args['outdir'], f)
             for f in ['dataset.libsvm', 'dataset_nominal.arff', 'dataset_nominal_fixed.arff', 'read_names.txt',
                       'prediction.txt', 'debug_prediction.txt', 'debug_prediction_avg.txt']]))]
    return _pipeline(chain)