#!/bin/bash
testdir="$(dirname $0)"
if [ -n "$1" ]; then
  python=$1
else
  python=python
fi
root="../"
export PYTHONPATH="$root/src:$PYTHONPATH"
cd $testdir
$python -m unittest discover 
