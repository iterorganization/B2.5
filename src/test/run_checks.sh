#!/bin/bash
if [ ! -d "$1" ]; then
  echo "Error, first argument should be a directory!"
  exit
fi
if [ ! -d "$2" ]; then
  echo "Error, second argument should be a directory!"
  exit
fi

MYPATH=`dirname "$0"` # path to this script, used when we calling other scripts from the same directory

files="b2fmovie b2fparam b2fplasma b2fstate b2ftrace b2ftrack"

#for f in $files; do 
#   filenames="$1/$f $2/$f"; 
#   check_b2_output $filenames
#done  > compare_results.log

missing=0
for f in $files; do
   if [ ! -f $1/$f ]; then
     echo "Error, file not found $1/$f";
     missing=$((missing+1))
   elif [ ! -f $2/$f ]; then
     echo "ERROR, file not found $2/$f";
     missing=$((missing+1))
   else
     echo "Checking file $f"
     check_b2_output $1/$f $2/$f
   fi
done > compare_results.log

if [ $missing -gt 0 ]; then
  echo Results are NOT correct, files missing.
  exit 1
fi

# For most of the variables we compare the maximum error to the average array value
# (and this average is the average of abs(var)).
# We have a large tolerance for most of the variables (1e-6).
# For basic quantities like te, ti, na, ni, ne, po, ua, we use a stricter tolerance level.
# Except for the velocity, we check the maximum relative error for the basic quantities.
$MYPATH/b2diff.py --tolerance 1e-6 --maxerr 'te ti na ni ne po' --specific-tolerance 'ua 1e-11 ne 1e-11 ni 1e-11 na 1e-8 te 1e-11 ti 1e-11 po 1e-10' -i 'time|del*|res*|b2stbc_smo|b2stbc_sna|b2stbc_she|b2stbc_shi|fllim0fna' -v compare_results.log

STATUS=$? # exit status of b2diff.py
# The exit status tells whether the test were successfull
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS
