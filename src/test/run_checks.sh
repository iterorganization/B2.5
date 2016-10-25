#!/bin/bash
if [ ! -d "$1" ]; then
  echo "Error, first argument should be a directory!"
  exit
fi
if [ ! -d "$1/b2mn.exe.dir" ]; then
  echo "Error, $1/b2mn.exe.dir not found!"
  exit
fi
if [ ! -d "$2" ]; then
  echo "Error, second argument should be a directory!"
  exit
fi
if [ ! -d "$2/b2mn.exe.dir" ]; then
  echo "Error, $2/b2mn.exe.dir not found!"
  exit
fi

MYPATH=`dirname "$0"` # path to this script, used when we calling other scripts from the same directory

for f in $1/b2mn.exe.dir/b2f*; do 
   filenames="$f $2/b2mn.exe.dir/`basename $f`"; 
   $MYPATH/check_b2_output $filenames  
done  > compare_results.log

# For most of the variables we compare the maximum error to the average array value
# (and this average is the average of abs(var)).
# We have a large tolerance for most of the variables (1e-7).
# For basic quantities like te, ti, na, ni, ne, po, ua we use more strick tolerance level.
# Except for the velocity, the we check the maximum relative error for the basic quantities'
$MYPATH/b2diff.py --tolerance 1e-7 --maxerr 'te ti na ni ne po' --specific-tolerance 'ua 1e-12 ne 1e-12 ni 1e-12 na 1e-10 te 1e-12 ti 1e-12 po 1e-12' -i 'del*|res*' -v compare_results.log

STATUS=$? # exit status of b2diff.py
# The exit status tells whether the test were successfull
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS