#!/bin/bash
set -xe

python setup.py bdist_wheel
python2 setup.py bdist_wheel

rm -rf build/ src/asydo.egg-info