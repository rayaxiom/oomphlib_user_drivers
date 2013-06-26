#!/bin/bash



NPROC=$1
NOEL=$2
FILE="sq_po"

RAYMACHINE=$HOSTNAME # wulfling.maths.manchester.ac.uk, rice, cake...
OOMPHBASE=$(make -s --no-print-directory print-top_builddir)
CURRENTDIR=`pwd`

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
V="1"
A="0"
R="200"

ITSTIMEDIR="res_ray"

rm -rf $ITSTIMEDIR && mkdir $ITSTIMEDIR

function run_code {
if [ "$NPROC" = "0" ]; then
./$FILE --w_solver $WS --ns_solver $NS --visc $V \
          --ang $A --rey $R --noel $NOEL
else
mpirun -np $NPROC ./$FILE \
  --w_solver $WS --ns_solver $NS --p_solver 1 --f_solver 69 \
  --visc $V --ang $A --rey $R --noel $NOEL --itstimedir $ITSTIMEDIR

#mpirun -np $NPROC ./square0 \
#  --w_solver $WS --ns_solver $NS --p_solver 1 --f_solver 69 \
#  --visc $V --ang $A --rey $R --noel $NOEL
#NS="1"
#F="0"
#P="0"
#mpirun -np $NPROC ./square0 \
#  --w_solver $WS --ns_solver $NS --f_solver $F --p_solver $P \
#  --visc $V --ang $A --rey $R --noel $NOEL
fi
}

cd $OOMPHBASE/src/ && make && make install && \
cd $CURRENTDIR && \
make $FILE && \
run_code



