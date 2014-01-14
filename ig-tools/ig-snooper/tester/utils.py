import os
import subprocess
from compare_marking import compare_marking
from parsers import kabat
from test_data_generator import split, get_out_filenames

DATA_DIR = 'data'
TEMP_DIR = 'tmp'
SNOOPER_DIR = 'snooper'
TRAIN_COMMAND = './train.sh'
PREDICT_COMMAND = './predict.sh'
# TRAIN_COMMAND = 'train.bat'
# PREDICT_COMMAND = 'predict.bat'


def train(fasta_path, kabat_path, ml_window_size, io_dir):
    with open('snooper.log', 'a') as snooper_log:
        subprocess.Popen(
            ' '.join([TRAIN_COMMAND, os.path.abspath(fasta_path), os.path.abspath(kabat_path), str(ml_window_size),
             os.path.abspath(io_dir)]),
            cwd=SNOOPER_DIR, stdout=snooper_log, stderr=snooper_log, shell=True).wait()


def test(fasta_path, kabat_path, ml_window_size, avg_window_size, merge_threshold, io_dir):
    with open('snooper.log', 'a') as snooper_log:
        subprocess.Popen(' '.join([PREDICT_COMMAND, os.path.abspath(fasta_path), str(ml_window_size), str(avg_window_size),
                          str(merge_threshold), os.path.abspath(io_dir)]),
                         cwd=SNOOPER_DIR, stdout=snooper_log, stderr=snooper_log, shell=True).wait()

    with open(kabat_path) as ref_kabat, open(os.path.join(io_dir, 'results.kabat')) as pred_kabat:
        reference = kabat.parse(ref_kabat)
        prediction = kabat.parse(pred_kabat)

    errs = compare_marking(reference, prediction)
    return errs['wrong_rate'], errs['error_rate'], errs['misfit_rate'], errs['error_wrong_rate'], errs[
        'misfit_wrong_rate']


def experiment(fasta_path, kabat_path, ml_window_size, avg_window_size, merge_threshold, count, out_dir):
    split(fasta_path, kabat_path, out_dir, .3, count)
    train_fasta, test_fasta = get_out_filenames(out_dir, fasta_path)
    train_kabat, test_kabat = get_out_filenames(out_dir, kabat_path)
    train(train_fasta, train_kabat, ml_window_size, out_dir)
    # return test(test_fasta, test_kabat, ml_window_size, avg_window_size, merge_threshold, out_dir)
    return test('vh.fasta', 'vh.kabat', ml_window_size, avg_window_size, merge_threshold, out_dir)

def log_result(file, ml_window, avg_window, merge_threshold, *errs):
    file.write('{0} {1} {2} {3}\n'.format(ml_window, avg_window, merge_threshold, ' '.join(str(e) for e in errs)))
