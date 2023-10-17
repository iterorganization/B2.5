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
     check_b2_output $1/$f $2/$f
   fi
done  > compare_results.log

if [ $missing -gt 0 ]; then
  echo Results are NOT correct, files missing.
  exit 1
fi

# For most of the variables we compare the maximum error to the average array value
# (and this average is the average of abs(var)).
# To avoid statistical variance, this test must be done using correlated sampling and use the APCAS Eirene MPI parallelization strategy.
# Except for the velocity, we check the maximum relative error for the basic quantities.
$MYPATH/b2diff.py --tolerance 0.01 --maxerr 'te ti na ni ne po' --specific-tolerance 'fhe 0.02 fhe_eir 0.02 fhep 0.02 fhe0 0.02 fhe_mdf 0.02 fht 0.02 fhm 0.02 rqahe 0.05 rqrad 0.05 ti 0.02 ne2 0.02 hce0 0.05 hci0 0.05 kinrgy 0.05' -i 'time|data|b2stb*|res*|del*|sm*|na*|ne0|nep|ni0|ua*|rcx*|rsa*|cdpa|csigin|cvsa*|dpa*|fllimvisc|b2sihs_visa|b2sihs_divu*|b2npmo_sm*' -v compare_results.log

STATUS=$? # exit status of b2diff.py
# The exit status tells whether the test were successful
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS
