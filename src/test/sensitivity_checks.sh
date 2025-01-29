#!/bin/bash
if [ ! -f "$1" ]; then
  echo "Error, first argument should be a file!"
  exit
fi
if [ ! -f "$2" ]; then
  echo "Error, second argument should be a file!"
  exit
fi

MYPATH=`dirname "$0"` # path to this script, used when we calling other scripts from the same directory

$MYPATH/ADdiff.py --tolerance $4 --ndata $3 $1 $2

STATUS=$? # exit status of ADdiff.py
# The exit status tells whether the test were successfull
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS
