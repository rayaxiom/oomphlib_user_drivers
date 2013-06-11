#!/bin/bash
# Used to test step_po
# Usage: ./run_steppo.sh 1 4
# 1 core, 4 elements in 1D

NPROC=$1
NOEL=$2
PROGRAM="step_po"

RAYMACHINE=$HOSTNAME # wulfling.maths.manchester.ac.uk, rice, cake...
OOMPHTOPDIR=$(echo `grep abs_top_srcdir Makefile` | sed -e "s/abs_top_srcdir = //g")
#OOMPHROOT=$(grep "abs_top_srcdir" Makefile)
#OOMPHROOT=${OOMPHROOT:17}
CURRENT_DIR=`pwd`

# Flags:
# --doc_soln - if present, doc the solution.
# --doc_prec");
# --w_solver", &myvar.W_solver);
# --ns_solver", &myvar.NS_solver);
# --p_solver", &myvar.P_solver);
# --f_solver", &myvar.F_solver);
# --visc", &myvar.Vis);
# --ang", &myvar.Ang);
# --rey", &myvar.Rey);
# --noel", &myvar.Noel);
# --sigma",
# --bdw");
# --amg_str", &RayParam::amg_strength);
# --amg_damp", &RayParam::amg_damping);


WS="0"
NS="1"
F="0" # 0, 69
P="0" # 0,1
V="0"
A="30"
R="-1"

function run_code {
if [ "$NS" = "0" ]; then
./$FILE --w_solver $WS --ns_solver $NS --visc $V \
          --ang $A --rey $R --noel $NOEL
else
  if [ "$NS" = "1" ]; then
    mpirun -np $NPROC ./$FILE \
      --w_solver $WS --ns_solver $NS --p_solver $P --f_solver $F \
      --visc $V --ang $A --rey $R --noel $NOEL
  else
    mpirun -np $NPROC ./$FILE \
      --w_solver $WS --ns_solver $NS \
      --visc $V --ang $A --rey $R --noel $NOEL
  fi
fi
}

cd $OOMPHTOPDIR/src/ && make && make install && \
cd $CURRENT_DIR && \
make $PROGRAM && \
run_code



