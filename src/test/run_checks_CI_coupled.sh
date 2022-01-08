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
# To allow for statistical variance, we have a large tolerance for most of the variables (0.125).
# Except for the velocity, we check the maximum relative error for the basic quantities.
$MYPATH/b2diff.py --tolerance 0.15 --maxerr 'te ti na ni ne po' -i 'time|data|na|ni0|nep|ne0|ua|ua0|kinrgy|fne*|fmo|fhi_eir|fch_p|del*|res*|rrahi|rsahi|b2stbr_*|b2stbc_*|b2sihs_visa|b2sihs_divua|b2sihs_divue|b2sihs_fraa|b2npmo_*|rcx*|smfr|smq|sch|she|shi|csig*|calf|alf0|hce0|hci0|sig0|dpa*|cdpa|cvsa*|fllim*|floe_noc|floi_noc|ne2' -v compare_results.log

STATUS=$? # exit status of b2diff.py
# The exit status tells whether the test were successful
# if [ "$STATUS" == 0 ]; then
#   echo "Results agree"
# else
#   echo "Error, results do NOT agree"
# fi

exit $STATUS
