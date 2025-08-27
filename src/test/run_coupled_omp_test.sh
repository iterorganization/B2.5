#!/bin/bash
# Script to run B2.5-Eirene OpenMP test using 1 and [max_threads] threads, compare the results.
#
# Usage:
# run_omp_test benchmark rundir [steps] [max_threads]
#
# Arguments:
# benchmark -- the name of a benchmark like AUG_16151_D+C+He, ITER_535_D+He+Ar, etc.
# rundir    -- is the subdirectory specific for the test case, a few examples:
#
#   ./run_omp_test.sh AUG_16151_D+C+He run_coupled_OpenMP_for_CI
#   ./run_omp_test.sh ITER_535_D+He+Ar standalone
#   ./run_omp_test.sh ITER_2171_D+He+Be+Ne standalone
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
  if [ -s b2.feedback_save.parameters.i ]; then
    mv b2.feedback_save.parameters.i b2.feedback_save.parameters
  fi
  if [ -s b2.sputter_save.parameters.i ]; then
    mv b2.sputter_save.parameters.i b2.sputter_save.parameters
  fi
  if [ -s b2.wall_save.parameters.i ]; then
    mv b2.wall_save.parameters.i b2.wall_save.parameters
  fi
  if [ -s b2time.nc.i ]; then
    mv b2time.nc.i b2time.nc
  fi
  if [ -s b2batch.nc.i ]; then
    mv b2batch.nc.i b2batch.nc
  fi
  if [ -s b2tallies.nc.i ]; then
    mv b2tallies.nc.i b2tallies.nc
  fi
  if [ -s fort.10.i ]; then
    mv fort.10.i fort.10
  fi
  if [ -s fort.11.i ]; then
    mv fort.11.i fort.11
  fi
  if [ -s fort.13.i ]; then
    mv fort.13.i fort.13
  fi
  if [ -s fort.14.i ]; then
    mv fort.14.i fort.14
  fi
  if [ -s fort.15.i ]; then
    mv fort.15.i fort.15
  fi
  if [ -s fort.44.i ]; then
    mv fort.44.i fort.44
  fi
  if [ -s fort.46.i ]; then
    mv fort.46.i fort.46
  fi
  touch b2fstati
  setup_baserun_eirene_links
  # Change the number of time steps
  sed -i "s/\('b2mndr_ntim'\s*\)'[0-9]\+'/\1'$3'/" b2mn.dat

  # Turn off excessive I/O for the 98 species test case
  sed -i "s/\('b2npmo_iout' *\)'.'/\1'0'/" b2mn.dat

  # Set number of threads, compact pinning, and stack size
  export OMP_NUM_THREADS=$4
  export KMP_AFFINITY=verbose,norespect,compact
  if [ -z "$KMP_STACKSIZE" ]; then
    export KMP_STACKSIZE=128MB
  fi
  if [ -z "$OMP_STACKSIZE" ]; then
    if [ "$COMPILER" != "ifort64" ]; then
      export OMP_STACKSIZE=128MB
    fi
  fi
  if [ -n "$SOLPS_MPI" ]; then
    echo "Found MPI mode"
    MPI_OPTS='-m "mpiexec -np 1"'
    echo b2run $MPI_OPTS b2mn
  fi

  if [ $4 -eq 1 ]; then
    echo Running B2.5-Eirene in serial mode
    b2run $MPI_OPTS b2mn 2>&1 | tee run.log
  else
    echo Running B2.5-Eirene using $4 threads
    b2run $MPI_OPTS -t $4 b2mn 2>&1 | tee run.log
  fi
  cd ..
}

DATESTAMP=$(datestamp)

# Run B2.5-Eirene with 1 and MAX_THREADS cores
run_test $SOURCEDIR test1_$DATESTAMP $STEPS 1
run_test $SOURCEDIR test2_$DATESTAMP $STEPS $MAX_THREADS

# Compare the results
# run_checks will call check_b2_output, make sure it is compiled
# cd $SOLPSTOP/modules/B2.5/src/test/
# ifort -g -O2 check_b2_output.F90 -o check_b2_output
$SOLPSTOP/modules/B2.5/src/test/run_checks_OpenMP_CI_coupled.sh $SOURCEDIR test1_$DATESTAMP
$SOLPSTOP/modules/B2.5/src/test/run_checks_OpenMP_CI_coupled.sh test1_$DATESTAMP test2_$DATESTAMP
