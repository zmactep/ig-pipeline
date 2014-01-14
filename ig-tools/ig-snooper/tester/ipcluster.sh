#!/bin/bash
export PYTHONPATH=/home/mactep/DEV/production/snooper-tester/dependencies/common_lib/python:/home/mactep/DEV/production/snooper-tester/dependencies/ig-utils:$PYTHONPATH
ipcluster3 start --n=4
