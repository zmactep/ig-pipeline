import sys
from splitter.splitter import run

input_file = sys.argv[1]
results_dir = sys.argv[2]
run(input_file, results_dir)
