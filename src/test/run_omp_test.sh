#!/bin/bash
# Script to run B2 OpenMP test using 1 and [max_threads] threads, compare the results.
#
# Usage:
# run_omp_test benchmark rundir [steps] [max_threads]
#
# Arguments:
# benchmark -- the name of a benchmark like AUG_16151_D+C+He, ITER_535_D+He+Ar, etc.
# rundir    -- is the subdirectory specific for the test case, a few examples:
#
#   ./run_omp_test.sh AUG_16151_D+C+He 16151_1.6MW_2.0e19_D=0.4_chi=1.6_standalone
#   ./run_omp_test.sh ITER_535_D+He+Ar run_as_5.2
#   ./run_omp_test.sh ITER_Be-W_D+T+He+Ne testrun_coupled
#
# Optional arguments:
# steps       -- number of time steps to take (default 1)
# max_threads -- number of threads used for the second test (default 48)
#
# Returns 0 if the results agree, otherwise return 1.

set -e    # exit on error

if [ -z "$SOLPSTOP" ]; then
  echo "Error, SOLPSTOP is not set"
  exit 1
fi

cd $SOLPSTOP/runs/examples

if [ -z "$1" ]; then
  echo "Give a testcase name as an argument. The available testcases are"
  make help
  exit 1
fi

if [ ! -d $1 ]; then
  make $1
fi

cd $1

if [ "$1" = "ITER_2588_D+He+N" ]; then
  ln -s baserun_N_atom baserun
fi

correct_baserun_timestamps

if [ ! -d "$2" ]; then
  echo "Error, $2 not found, available directories for $1 are"
  ls -d -1 */
  exit 1
else
  SOURCEDIR=$2
fi

if [ -n "$3" ]; then
  STEPS=$3
else
  STEPS=1
fi

if [ -n "$4" ]; then
  MAX_THREADS=$4
else
  MAX_THREADS=48
fi

echo We will run B2.5 for $STEPS steps.

function run_test {
  # run_test sourcedir targetdir nsteps nthreads
  # run test based on sourcedir using nthreads
  if [ ! -d $1 ]; then
    echo "Error, $1 is not a directory"
    exit 1
  fi
  if [ -z "$2" ]; then
    echo "error, give target directory name"
    exit 1
  fi
  rsync -avP $1/. $2
  cd $2
  touch b2fstati
  # change the number of time steps
  sed -i "s/\('b2mndr_ntim'\s*\)'[0-9]\+'/\1'$3'/" b2mn.dat

  # disable eirene
  sed -i "s/\('b2mndr_eirene'\s*\)'1'/\1'0'/" b2mn.dat

  # Turn off excessive io for the 98 species test case
  sed -i "s/\('b2npmo_iout' *\)'.'/\1'0'/" b2mn.dat

  #set number of threads, compact pinning, and stack size
  export OMP_NUM_THREADS=$4
  export KMP_AFFINITY=verbose,norespect,compact
  export KMP_STACKSIZE=32MB

  echo running B2 using $4 threads
  make -f $SOLPSTOP/runs/Makefile STAND_ALONE=yes  b2mn.prt 2>&1 | tee b2run.log
  cd ..
}

DATESTAMP=$(datestamp)

# run B2 with 1 and 48 cores
run_test $SOURCEDIR test1_$DATESTAMP $STEPS 1
run_test $SOURCEDIR test2_$DATESTAMP $STEPS $MAX_THREADS

# compare the results
# run_checks will call check_b2_output, make sure it is compiled
# cd $SOLPSTOP/modules/B2.5/src/test/
# ifort -g -O2 check_b2_output.F90 -o check_b2_output
$SOLPSTOP/modules/B2.5/src/test/run_checks.sh test1_$DATESTAMP test2_$DATESTAMP
