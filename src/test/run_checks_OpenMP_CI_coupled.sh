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
# To avoid statistical variance, this test must be done using correlated sampling.
# The NLIDENT switch in Eirene must be turned on.
# For each stratum, the NPTSDEL number must be equal to the NPTS value divided by some
# integer multiple of the number of threads used.
# Except for the velocity, we check the maximum relative error for the basic quantities.
$MYPATH/b2diff.py --tolerance 0.01 --maxerr 'te ti na ni ne po' --specific-tolerance 'fhe 0.1 fhe_eir 0.1 fhep 0.10 fhe0 0.05 fhe_mdf 0.1 fhi 0.05 fhi0 0.05 fhi_eir 0.05 fhet 0.25 fhit 0.05 fhi_mdf 0.05 fhip 0.05 fht 0.05 fhm 0.15 fna0 0.10 fmo 0.10 fne_eir 0.1 te 0.1 ti 0.1 ti0 0.1 tip 0.1 po 0.25 pop 0.25 ne2 0.02 calf 0.10 chce 0.25 hce0 0.10 chci 0.10 hci0 0.10 fne_32 0.1 fne_52 0.1 fch_32 0.02 fch_52 0.025' -i 'time|data|b2stb*|res*|del*|sm*|po0|na*|ne0|nep|ni0|ua*|kinrgy|fna_he|fna_52|fch|fch0|fchp|fch_p|fhj|fne|fni*|rcx*|rqahe|rqrad|rra*|rsa*|alf*|calf_an|cdpa|csig*|dpa*|fllim*|b2sihs_*|b2npmo_sm*' -v compare_results.log

STATUS=$? # exit status of b2diff.py
# The exit status tells whether the test were successful
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS
