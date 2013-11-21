import argparse
import json
from pipeline import train_pipeline
import sys

def main():
    train_pipeline(get_params())

def get_params():
    parser = argparse.ArgumentParser(description='Train model given train dataset')
    parser.add_argument('--config_path', nargs='?', help='Path to config')
    parser.add_argument('--fasta', help='Train data in fasta')
    parser.add_argument('--kabat', help='Train data in kabat')
    parser.add_argument('--model_name', help='Output model name')
    parser.add_argument('--ml_window_size', type=int, nargs='?', help='ML windows size')
    parser.add_argument('--tools_root', nargs='?', help='Tools dir')
    parser.add_argument('--outdir', nargs='?', help='Output directory')

    params = {k: v for k, v in vars(parser.parse_args()).items() if v is not None}

    if 'config_path' in params and not params['config_path'] is None:
        config = json.load(open(params['config_path']))
        config.update(params)  # commandline args has higher priority
        params = config

    required = ('fasta', 'kabat', 'model_name', 'ml_window_size', 'tools_root', 'outdir')
    if not all(key in params and params[key] is not None for key in required):
        print('Missing required params. Make sure that you have all of these: %s' % ' '.join(required))
        sys.exit()

    return params

if __name__ == "__main__":
   main()