#!/bin/bash

# SET THIS.
PROGRAM="sqcurve_van"

# On interrupt INT, i.e. Ctrl c, we exit.
# Without this, we cannot exit the loop, it keeps on continuing.
# if we "set -e", the script will exit if the PROGRAM crashes. We do not want this.
trap 'exit' INT

# How to use:
# (./loop_sqcurve_van.sh 2>&1) | tee -a loop_sqcurve_van.output
#
# What this does:
#
# 1) The function gen_tests generates a list of tests to the file TEST_LIST 
#    i.e. loop_sq_po.testlist
#
# 2) The function run_tests will read in each test from loop_sq_po.testlist
#    into the variable "RUNCOMMAND" and run them. 
#    Output from each test is redirected to TESTOUT_FILE
# 
# 3) Once the test is completed, we echo the "$RUNCOMMAND" into TESTLISTFIN_FILE
#    This file is used to determine where we continue from.
#   
#   List of files:
#
#   loop_xxx.sh: this current run script.
#   loop_xxx.output: output from this run script. If any test fails, we can see them here.
#   ./res_loop_xxx/ the results I'm interested, within this dir we also have:
#     loop_xxx.testlist: All the tests to run. This is generated by the function gen_tests
#     loop_xxx.testlistfin: A list of finished tests so far. Used for continuing tests later.
#     loop_xxx.testoutput: the stdout from all the tests.
#
# IMPORTANT:
#
#   * On the first run, the loop_xxx.output file should be removed.
#
#   * THE FIRST RUN: If the res_loop_xxx directory does not exist, we start off fresh.
#       The directory is created and the executable is re-built, 
#       The three files .testlist, .testlistfin, and .testoutput are regenerated.
#       The function run_tests is called.
#   
#   * CONTINUATION: If the res_loop_xxx directory exists, 
#     we assume that continuation is wanted.
#       The executable is NOT re-built, if it does not exist we stop.
#       If the executable exists, then the function run_tests_con is called.
#       

FILE=$0 # This contains "./", which we do not want.
FILE=${FILE:2} # Gets rid of "./"
FILEBASE=${FILE%%.*} # Get rid of the extension (in this case, ".sh")
CURRENTDIR=`pwd`
OOMPHROOT=$(make -s --no-print-directory print-top_builddir)

ITSTIMEDIR="res_$FILEBASE"
TESTLIST_FILE="$ITSTIMEDIR/$FILEBASE.testlist"
TESTLISTFIN_FILE="$ITSTIMEDIR/$FILEBASE.testlistfin"
TESTOUT_FILE="$ITSTIMEDIR/$FILEBASE.testoutput"

# The output file for this script (if it exits)
SHELLOUTPUT_FILE="${FILEBASE}.output"

gen_tests()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F AMG

PRECLIST="0 1 2"
# The precs are set according to the PRECLIST above.
NSPRECLIST="0" # 0 - Exact, 1 - LSC
PPRECLIST="0" # Only for LSC, 0 - Exact, 1 - AMG
FPRECLIST="0" # Only for LSC, 0 - Exact, 1 - AMG

VISLIST="0 1"
ANGLIST="0"
RELIST="0 100 200"
NOELLIST="4 8 16 32 64 128 256"

for PREC  in $PRECLIST
do
  case "$PREC" in
    0)
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
echo "mpirun -np 1 ./$PROGRAM --ns_solver $NSPREC --p_solver $PPREC --f_solver $FPREC --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $ITSTIMEDIR" >> $TESTLIST_FILE
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
echo "mpirun -np 1 ./$PROGRAM --ns_solver $NSPREC --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $ITSTIMEDIR" >> $TESTLIST_FILE
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
  # State that this is the initial test run (as opposed to run_tests_con)
  # Total number of tests
  NUMTESTS=$(wc $TESTLIST_FILE | awk '{print $1;}')
  echo "RAYRAY: initial test run, total tests: $NUMTESTS"
  
  # Check if SHELLOUTPUT_FILE contains previous runs.
  LINESINSHELLOUT=$(grep "RAYRAY: initial" $SHELLOUTPUT_FILE | wc | awk '{print $1;}')

  if [ -f $SHELLOUTPUT_FILE ]; then
    if [ "$LINESINSHELLOUT" -ne "1" ]; then
      echo "This is the first run but you have not removed the previous $SHELLOUTPUT_FILE file."
      echo "Please delete both the file $SHELLOUTPUT_FILE and the directory $ITSTIMEDIR"
      echo "Then re-run this script."
      exit 1
    fi
  fi


  # Current test number.
  TESTNUM=1
  
  while read RUNCOMMAND; do
    echo "Test $TESTNUM/$NUMTESTS: \"$RUNCOMMAND\" on $(date)"
    bash -c "$RUNCOMMAND" </dev/null >> $TESTOUT_FILE

    # For some reason I have to echo to stdout, otherwise it will execute
    # the rest of the loop, updating TESTLISTFIN_FILE even when I ctrl+c
    echo ""
    TESTNUM=$[$TESTNUM+1]
    echo "$RUNCOMMAND" >> $TESTLISTFIN_FILE
  done < $TESTLIST_FILE
}

run_tests_con()
{
  # Current test number.
  TESTNUM=1

  # How many tests have been done so far?
  # Count the number of files in the ITSTIMEDIR directory
  NUMTESTFIN=$(wc $TESTLISTFIN_FILE | awk '{print $1;}')
  # Total number of tests
  NUMTESTS=$(wc $TESTLIST_FILE | awk '{print $1;}')

  echo "RAYRAY: continuing, done $NUMTESTFIN/$NUMTESTS tests so far..."

  while read RUNCOMMAND; do
    if [ "$TESTNUM" -le "$NUMTESTFIN" ]; then
      TESTNUM=$[$TESTNUM+1]
    else
      echo "Test $TESTNUM/$NUMTESTS: \"$RUNCOMMAND\" on $(date)"
      bash -c "$RUNCOMMAND" </dev/null >> $TESTOUT_FILE

      # For some reason I have to echo to stdout, otherwise it will execute
      # the rest of the loop, updating TESTLISTFIN_FILE even when I ctrl+c
      echo ""
      TESTNUM=$[$TESTNUM+1]
      echo "$RUNCOMMAND" >> $TESTLISTFIN_FILE
    fi
  done < $TESTLIST_FILE
}

###############################################################################
############### "main function"-like code starts here #########################
###############################################################################

# IF the results directory exists, 
#   try to continue the tests in TESTLIST_FILE, 
#   assuming that the executable PROGRAM exists.
#   IF the executable does not exist, we echo out a message and do nothing.
# ELSE (the results directory does not exist)
#   build the executable
#   create the results directory
#   generate the xxx.testlist
#   run the tests (output in xxx.testoutput)
if [ -d "$ITSTIMEDIR" ]; then
  if [ ! -f $PROGRAM ]; then
    echo "The executable $PROGRAM does not exist but the result directory"
    echo "$ITSTIMEDIR does. Please check that the .cc file is the same before"
    echo "compiling and CONTINUING the tests."
    echo " "
    echo "Otherwise, delete the result directory $ITSTIMEDIR" and re-run this
    echo "script to compile the executable and run all the tests."
    
    exit 1
  else
    run_tests_con
    
    exit 0
  fi
else
  make $PROGRAM && rm -rf $ITSTIMEDIR && mkdir $ITSTIMEDIR &&
  gen_tests && run_tests
  
  exit 0
fi

