import argparse
import json
import logging
import sys
from pipeline import predict_pipeline


def main():
    predict_pipeline(get_params())

def get_params():
    parser = argparse.ArgumentParser(description='Predict data given a model')
    parser.add_argument('--config_path', nargs='?', help='Path to config')
    parser.add_argument('--fasta', help='Predict data in fasta')
    parser.add_argument('--kabat', help='Optional kabat for benchmarking')
    parser.add_argument('--model_path', help='input model path')
    parser.add_argument('--ml_window_size', type=int, nargs='?', help='ML windows size')
    parser.add_argument('--avg_window_size', type=int, nargs='?', help='averaging windows size')
    parser.add_argument('--merge_threshold', type=int, nargs='?', help='merge small regions if len < merge_threshold')
    parser.add_argument('--tools_root', nargs='?', help='Tools dir')
    parser.add_argument('--outdir', nargs='?', help='Output directory')

    params = {k: v for k, v in vars(parser.parse_args()).items() if v is not None}

    if 'config_path' in params and not params['config_path'] is None:
        config = json.load(open(params['config_path']))
        config.update(params)  # commandline args has higher priority
        params = config

    required = ('fasta', 'model_path', 'ml_window_size', 'avg_window_size', 'merge_threshold', 'tools_root', 'outdir')
    if not all(key in params and params[key] is not None for key in required):
        print('Missing required params. Make sure that you have all of these: %s' % ' '.join(required))
        sys.exit()

    return params

if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.critical("Terminated.")