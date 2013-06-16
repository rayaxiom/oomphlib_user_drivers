#!/bin/bash
# Used to test step_van
# Usage: ./run_stepvan.sh 1 4
# 1 core, 4 elements in 1D

NPROC=$1
NOEL=$2
PROGRAM="step_van"

RAYMACHINE=$HOSTNAME # wulfling.maths.manchester.ac.uk, rice, cake...
OOMPHTOPDIR=$(echo `grep abs_top_srcdir Makefile` | sed -e "s/abs_top_srcdir = //g")
#OOMPHROOT=$(grep "abs_top_srcdir" Makefile)
#OOMPHROOT=${OOMPHROOT:17}
CURRENT_DIR=`pwd`

# Flags:
# --doc_soln - if present, doc the solution.
# --doc_prec");
# --ns_solver", &myvar.NS_solver);
# --p_solver", &myvar.P_solver);
# --f_solver", &myvar.F_solver);
# --visc", &myvar.Vis);
# --ang", &myvar.Ang);
# --rey", &myvar.Rey);
# --noel", &myvar.Noel);
# --amg_str", &RayParam::amg_strength);
# --amg_damp", &RayParam::amg_damping);

NS="0"
F="1"
P="1"
V="0"
A="0"
R="-1"

function run_code {
if [ "$NPROC" = "0" ]; then
./$PROGRAM --ns_solver $NS --f_solver $F --p_solver $P \
        --visc $V --ang $A --rey $R --noel $NOEL
else
  if [ "$NS" = "1" ]; then
    mpirun -np $NPROC ./$PROGRAM --ns_solver $NS --f_solver $F --p_solver $P \
                              --visc $V --ang $A --rey $R --noel $NOEL
  else
    mpirun -np $NPROC ./$PROGRAM --ns_solver $NS \
                              --visc $V --ang $A --rey $R --noel $NOEL --doc_soln
  fi
fi
}

cd $OOMPHTOPDIR/src/ && make && make install && \
cd $CURRENT_DIR && \
make $PROGRAM && \
run_code


