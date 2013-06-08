#!/bin/bash

NPROC=$1
NOEL=$2

RAYMACHINE=$HOSTNAME # wulfling.maths.manchester.ac.uk, rice, cake...
OOMPHBASE="/home/ray/oomphlib/clean_checkout"
FILE="sq_van"

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


WS="0"
NS="1"
F="1"
P="1"
V="0"
A="0"
R="100"

function run_code {
if [ "$NPROC" = "0" ]; then
./$FILE --ns_solver $NS --f_solver $F --p_solver $P \
        --visc $V --ang $A --rey $R --noel $NOEL
else
mpirun -np $NPROC ./$FILE --ns_solver $NS --f_solver $F --p_solver $P \
                          --visc $V --ang $A --rey $R --noel $NOEL
fi
}

#cd $OOMPHBASE/src/ && make && make install && \
#cd $OOMPHBASE/user_drivers/lagrange_square/ && \
make $FILE && \
run_code



