#!/bin/bash
FILE="sq_van"

run_tests()
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
RELIST="0 100 200"
NOELLIST="4 8 16 32 64 128"

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
./$FILE --ns_solver $NSPREC --p_solver $PPREC --f_solver $FPREC --visc $VIS --ang $ANG --rey $RE --noel $NOEL
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
./$FILE --ns_solver $NSPREC --visc $VIS --ang $ANG --rey $RE --noel $NOEL
            done
          done
        done
      done
    fi
  done
done

} # run_tests function


#OOMPHPATH="/home/mly/oomphlib_current"
#cd $OOMPHPATH/src/ && make && make install && \
#cd $OOMPHPATH/user_drivers/lagrange_square/ && \
#make square2 && run_tests

run_tests

