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
$MYPATH/b2diff.py --tolerance 1e-6 --maxerr 'te ti na ni ne po' --specific-tolerance 'ua 0.01 ne 1e-11 ni 1e-11 na 0.02 ne 0.02 te 0.02 ti 0.01 tn 0.01 po 0.02 calf 0.01 ceqp 0.01 chce 0.01 chci 0.01 chvemx 0.01 chvimx 0.01 csig 0.01 csigin 0.01 dpa0 0.01 fhe 0.01 fhp 0.02 fhi 0.01 fhi_mdf 0.01 fht 0.02 floe_noc 0.01 floi_noc 0.01 fna_52 0.01 fna_he 0.01 fne 0.01 fne_52 0.025 fni 0.01 fni_52 0.025 fna_eir 0.01 fne_eir 0.01 hce0 0.01 hci0 0.01 kinrgy 0.01 sig0 0.01 she 0.01 vsa0 0.01' -i 'time|data|del*|res*' -v compare_results.log

STATUS=$? # exit status of b2diff.py
# The exit status tells whether the test were successful
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS
