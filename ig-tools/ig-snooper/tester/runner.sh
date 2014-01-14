#!/bin/bash
export PYTHONPATH=/home/mactep/DEV/production/snooper-tester/dependencies/common_lib/python:/home/mactep/DEV/production/snooper-tester/dependencies/ig-utils:$PYTHONPATH
rm snooper.log
python3 runner.py
