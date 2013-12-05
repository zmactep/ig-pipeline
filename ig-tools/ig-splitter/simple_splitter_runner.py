import logging
import shlex
import subprocess
from IPython.parallel import Client, require, error
import argparse
import os
from Bio import SeqIO
from itertools import chain
import time

def parse_args():
    parser = argparse.ArgumentParser(description='Split dataset sequences according to mids')
    parser.add_argument('filename', help='input filename')
    parser.add_argument('type', help='input file type, a string recognized by SeqIO.parse')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('--threads', type=int, default=4, help='number of threads to use (default = 4)')
    return parser.parse_args()


def reduce_dicts(ds):
    r = {}
    for d in ds:
        for k, v in d.items():
            if k in r:
                r[k].append(iter(v))
            else:
                r[k] = [iter(v)]
    for k in r:
        r[k] = chain(*r[k])
    return r


class ClusterHandle:
    def __init__(self, num_threads):
        self.client = None
        self._num_threads = num_threads

        logging.info('Launching ipcluster')

        try:
            subprocess.Popen(shlex.split('ipcluster3 start --n={0}'.format(num_threads)))
        except Exception:
            logging.error('Error in launching ipcluster')
            raise

        self._wait()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            subprocess.call(shlex.split('ipcluster3 stop'))
        except Exception:
            logging.error('Error in shutting down ipcluster')
            raise

    def _wait_controller(self):
        logging.info('Probing the controller ...')
        for i in range(5):
            try:
                self.client = Client()
                return
            except (ValueError, error.TimeoutError):
                time.sleep(5)
        logging.error('Could not connect to controller')
        raise RuntimeError()

    def _wait_engines(self):
        logging.info('Probing %i engines ...' % self._num_threads)
        running = len(self.client)
        while running < self._num_threads:
            time.sleep(1)
            previous = running
            running = len(self.client)

    def _wait(self):
        #  based on code from http://mail.scipy.org/pipermail/ipython-dev/2013-August/012189.html
        self._wait_controller()
        self._wait_engines()


def main():
    args = parse_args()
    
    with ClusterHandle(args.threads) as cluster:
        view = cluster.client.direct_view()

        @view.parallel()
        @require('simple_splitter')
        def split_dataset_parallel(recs):
            return simple_splitter.split_dataset(recs)

        @view.remote()
        @require('simple_splitter')
        def cleanup():
            simple_splitter.cleanup()

        try:
            records = SeqIO.parse(args.filename, args.type)
        except Exception:
            logging.error('Error opening input file')
            raise

        logging.info('Splitting dataset ...')
        try:
            async = split_dataset_parallel(records)
        except Exception:
            logging.error('Error in splitting dataset')
            raise

        async.wait_interactive()
        buckets = reduce_dicts(async.result)

        for k in buckets.keys():
            try:
                with open(os.path.join(args.outdir, k + '.fasta'), 'w') as file:
                    SeqIO.write(buckets[k], file, 'fasta')
            except Exception:
                logging.error('Error writing output file')
                raise

        try:
            cleanup()
        except Exception:
            logging.error('Error in cleaning up ipcluster engines')
            raise


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s %(asctime)s %(message)s')
    try:
        main()
    except Exception:
        logging.critical('Terminated.')