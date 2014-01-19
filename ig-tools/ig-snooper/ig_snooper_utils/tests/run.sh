#!/bin/bash

export PYTHONPATH=$PYTHONPATH:/opt/ig-pipeline/ig-tools/common_lib/python:/opt/ig-pipeline/ig-tools/ig-snooper:/opt/ig-pipeline/ig-tools/ig-snooper/ig_snooper_utils
time python test_ig_snooper.py /opt/ig-pipeline/data/test/germline /opt/ig-pipeline/ig-config/ig-snooper.conf
