"""
This module holds boilerplate code for construction of a generic pipeline of actions, and a complete pipeline for
running training and prediction of antibody regions using WEKA. It serves as a backend for two client-facing scripts,
"train.py" and "predict.py".
"""
from datetime import datetime
import os
from subprocess import Popen, PIPE
from ig_snooper_utils import svm_data_generator, fix_weka_header, parse_svm_output, compare_kabat

WEKA_PATH = 'common_lib/third_party/weka-3.6.10/weka.jar'


def _format_exception(e):
    return "{0}: {1}\n{2}".format(e.__class__.__name__, str(e), _format_exception(e.__cause__) if e.__cause__ else "")


# The following functions produce pipeline stages as closures:


def _get_check_files_action(pathlist):
    def action():
        for path in pathlist:
            if not os.path.exists(path):
                return False, "Path not exists: %s" % path
        return True, ""

    return action


def _get_convert_to_libsvm_action(args, generate_comments=False):
    def action():
        try:
            svm_data_generator.process(args['ml_window_size'], True, args['outdir'], 'dataset.libsvm',
                                       args['fasta'], args['kabat'] if 'kabat' in args else None,
                                       "read_names.txt" if generate_comments else None)
        except (ValueError, RuntimeError) as e:
            return False, "Failed to convert .fasta to .libsvm, exception detail:\n%s" % _format_exception(e)

        libsvm_path = os.path.join(args['outdir'], 'dataset.libsvm')
        if not (os.path.exists(libsvm_path) and os.path.getsize(libsvm_path) > 0):
            return False, 'Failed to convert .fasta to .libsvm, file is missing or empty'

        return True, ""

    return action


def _get_weka_conversion_action(args):
    def action():
        os.environ["CLASSPATH"] = os.path.join(args['tools_root'], WEKA_PATH)
        output, error = Popen(['java', '-Xmx4096M', 'weka.filters.unsupervised.attribute.NumericToNominal', '-i',
                               os.path.join(args['outdir'], 'dataset.libsvm'), '-o',
                               os.path.join(args['outdir'], 'dataset_nominal.arff')],
                              stdout=PIPE, stderr=PIPE).communicate()

        nominal_path = os.path.join(args['outdir'], 'dataset_nominal.arff')
        if not (os.path.exists(nominal_path) and os.path.getsize(nominal_path) > 0):
            return False, 'Failed to convert Numeric to Nominal: %s; %s' % (output, error)

        fix_weka_header.fix_header(os.path.join(args['outdir'], 'dataset_nominal.arff'),
                                   'dataset_nominal_fixed.arff', args['outdir'])

        return True, ""

    return action


def _get_weka_train_action(args):
    def action():
        output, error = Popen(['java', '-Xmx4096M', 'weka.classifiers.trees.RandomForest', '-I', '10',
                               '-K', '0', '-S', '1', '-no-cv', '-p', '0',
                               '-t', os.path.join(args['outdir'], 'dataset_nominal_fixed.arff'),
                               '-d', os.path.join(args['outdir'], 'model.model')],
                              stdout=PIPE, stderr=PIPE).communicate()

        model_path = os.path.join(args['outdir'], 'model.model')

        if not (os.path.exists(model_path) and os.path.getsize(model_path) > 0):
            return False, 'Failed to train: %s; %s' % (output, error)

        return True, ""

    return action


def _get_weka_predict_action(args):
    def action():
        prediction_file_name = os.path.join(args['outdir'], 'prediction.txt')
        with open(prediction_file_name, 'w') as out:
            Popen(['java', '-Xmx4096M', 'weka.classifiers.trees.RandomForest', '-no-cv', '-p', '0',
                   '-l', args['model_path'], '-T', os.path.join(args['outdir'], 'dataset_nominal_fixed.arff')],
                  stdout=out).wait()

        if not (os.path.exists(prediction_file_name) and os.path.getsize(prediction_file_name) > 0):
            return False, 'Failed to predict: output file is empty'

        parse_svm_output.parse(args['fasta'], prediction_file_name, os.path.join(args['outdir'], 'read_names.txt'),
                               args['outdir'], int(args['avg_window_size']), int(args['merge_threshold']))

        return True, ""

    return action


def _get_kabat_quality_metrics_action(args):
    def action():
        if not os.path.exists(args['kabat']):
            return False, "File does not exist: %s" % args['kabat']
        compare_kabat.compare(args['kabat'], os.path.join(args['outdir'], 'results.kabat'), args['outdir'])
        return True, ""
    return action


def _get_cleanup_action(pathlist):
    def action():
        for f in pathlist:
            os.remove(f)
        return True, ""

    return action


class _Stage():
    """
    Represents a stage in processing pipeline. Holds "desc", a message printed before execution, and "action", a closure
    doing the job, returning a tuple of (bool, string): the success flag and outcome message. On their meaning, read
    the _pipeline method description.
    """

    def __init__(self, desc, action):
        self._desc = desc
        self._action = action

    def run(self):
        print(self._desc)
        success, error_message = self._action()
        if not success:
            print(error_message)
        return success


def _pipeline(tasks):
    """
    Executes a series of tasks by iterating over provided collection of _Task objects. For every task, first prints its
    description, then calls the "run" method. It is expected to return a tuple of (success flag, outcome message).
    Processing continues over to the next task only if success flag is true, otherwise the outcome message gets
    printed and execution terminates. Return value is a boolean indicating successful completion of the pipeline overall
    """
    start_time = datetime.now()

    for task in tasks:
        if not task.run():
            print("Processing aborted on error")
            return False

    print('Done in %s' % str(datetime.now() - start_time))
    return True


def train_pipeline(args):
    """
    The pre-assembled pipeline for training of antibody region markup using WEKA. Args is a dictionary of command line
    arguments, such as the output of argparse.
    """
    chain = [
        _Stage("Checking input files ...", _get_check_files_action(
            [os.path.join(args['tools_root'], WEKA_PATH), args['fasta'], args['kabat'], args['outdir']])),
        _Stage("Converting data to libsvm format ...", _get_convert_to_libsvm_action(args)),
        _Stage("Applying NumericToNominal conversion ...", _get_weka_conversion_action(args)),
        _Stage("Training ...", _get_weka_train_action(args)) ]
    if 'clean_up' in args:
        chain += _Stage("Deleting unnecessary files ...", _get_cleanup_action(
            [os.path.join(args['outdir'], f)
             for f in ['dataset.libsvm', 'dataset_nominal.arff', 'dataset_nominal_fixed.arff']]))
    return _pipeline(chain)


def predict_pipeline(args):
    """
    The pre-assembled pipeline for prediction of antibody region markup using WEKA.
    """
    chain = [
        _Stage("Checking input files ...", _get_check_files_action(
            [os.path.join(args['tools_root'], WEKA_PATH), args['fasta'], args['model_path'], args['outdir']])),
        _Stage("Converting data to libsvm format ...", _get_convert_to_libsvm_action(args, generate_comments=True)),
        _Stage("Applying NumericToNominal conversion ...", _get_weka_conversion_action(args)),
        _Stage("Predicting ...", _get_weka_predict_action(args))]
    if 'kabat' in args:
        chain += _Stage("Evaluating quality metrics ...", _get_kabat_quality_metrics_action(args))
    if 'clean_up' in args:
        chain += _Stage("Deleting unnecessary files ...", _get_cleanup_action(
            [os.path.join(args['outdir'], f)
             for f in ['dataset.libsvm', 'dataset_nominal.arff', 'dataset_nominal_fixed.arff', 'read_names.txt',
                       'prediction.txt', 'debug_prediction.txt', 'debug_prediction_avg.txt']]))
    return _pipeline(chain)