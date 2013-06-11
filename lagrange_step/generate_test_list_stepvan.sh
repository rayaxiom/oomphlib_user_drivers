#!/bin/bash

# How to use:
# (./generate_test_list_stepvan.sh 2>&1) | tee generate_test_list_stepvan.dat
# First generates a list of tests into the file TEST_LIST i.e. test_list_step_van.dat
# The stdout from all the tests are stored in test_list_step_van.out
# The results I'm interested in should be in ./runsStepVan/ (this is set in the program)
# The results from this script is in generate_test_list_stepvan.dat,
#   in case we need to re-run some tests if they fail, or if we need to stop the
#   script at any point. We simply find which test it is currently doing in test_list_step_van.dat
#   and start from there.

PROGRAM="step_van"
OOMPHTOPDIR=$(echo `grep abs_top_srcdir Makefile` | sed -e "s/abs_top_srcdir = //g")
CURRENTDIR=`pwd`
TESTLIST_FILE="test_list_$PROGRAM.dat"
OUT_FILE="test_list_$PROGRAM.out"

gen_tests()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - Exact
# 1 - LSC: P SuperLU, F SuperLU
# 2 - LSC: P AMG, F AMG

PRECLIST="0 1 2"

# The precs are set according to the PRECLIST above.
NSPRECLIST="0" # 0 - Exact, 1 - LSC
PPRECLIST="0" # Only for LSC, 0 - Exact, 1 - AMG
FPRECLIST="0" # Only for LSC, 0 - Exact, 1 - AMG

VISLIST="0 1"
ANGLIST="0"
RELIST="-1"
NOELLIST="8 16 32 64"

for PREC  in $PRECLIST
do
  case "$PREC" in
    0)
    WPRECLIST="0"
    NSPRECLIST="0"
    ;;
    1)
    NSPRECLIST="1"
    PPRECLIST="0"
    FPRECLIST="0"
    ;;
    2)
    NSPRECLIST="1"
    PPRECLIST="1"
    FPRECLIST="1"
    ;;
    esac
  for NSPREC in $NSPRECLIST
  do
    if [ $NSPREC == 1 ]; then
      for PPREC in $PPRECLIST 
      do
        for FPREC in $FPRECLIST 
        do
          for VIS in $VISLIST
          do
            for ANG in $ANGLIST
            do
              for RE in $RELIST
              do
                for NOEL in $NOELLIST
                do
echo "mpirun -np 1 ./$PROGRAM --ns_solver $NSPREC --p_solver $PPREC --f_solver $FPREC --visc $VIS --ang $ANG --rey $RE --noel $NOEL" >> $TESTLIST_FILE
                done
              done
            done
          done
        done
      done
    else
      for VIS in $VISLIST 
      do
        for ANG in $ANGLIST
        do
          for RE in $RELIST
          do
            for NOEL in $NOELLIST
            do
echo "mpirun -np 1 ./$PROGRAM --ns_solver $NSPREC --visc $VIS --ang $ANG --rey $RE --noel $NOEL" >> $TESTLIST_FILE
            done
          done
        done
      done
    fi
  done
done

} # gen_tests function

run_tests()
{
  while read RUNCOMMAND; do
    echo "Doing $RUNCOMMAND on $(date)"
    bash -c "$RUNCOMMAND" </dev/null >> $OUT_FILE
  done < $TESTLIST_FILE
}

#cd $OOMPHTOPDIR/src/ && make && make install && \
#cd $CURRENTDIR && make $PROGRAM && \

rm -rf $TESTLIST_FILE && rm -rf $OUT_FILE && gen_tests && run_tests

