#!/bin/bash

NPROC=$1
NOEL=$2

RAYMACHINE=$HOSTNAME # wulfling.maths.manchester.ac.uk, rice, cake...
OOMPHBASE="/home/ray/oomphlib/clean_checkout"

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


WS=0
NS=0
V=0
A=30
R=100
function run_code {
if [ "$NPROC" = "0" ]; then
./square0 --w_solver $WS --ns_solver $NS --visc $V \
          --ang $A --rey $R --noel $NOEL
else
mpirun -np $NPROC ./square0 --w_solver $WS --ns_solver $NS \
                            --visc $V --ang $A --rey $R --noel $NOEL
fi
}

cd $OOMPHBASE/src/ && make && make install && \
cd $OOMPHBASE/user_drivers/lagrange_square/ && \
make square0 && \
run_code



#./square0 --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0 \
#          --visc 0 --ang 30 --rey 100 --noel 4

#./square0 --w_solver 0 --ns_solver 1 --visc 0 \
#          --ang 30 --rey 100 --noel 4


#if [ "$RAYTARGET" == rice]; then
#OOMPHPATH="/home/ray/learning/phd/wulfling/oomphlib_current/"
#  cd $OOMPHPATH/src/ && make && make install && \
#  cd $OOMPHPATH/user_drivers/lagrange_square/ && \
#  make square0 && \
#  ./square0 --w_solver 0 --ns_solver 1 --f_solver 1 --p_solver 1 \
#            --visc 0 --ang 30 --rey 100 --noel 4
#elif [ "$RAYTARGET" == "wulfling.maths.manchester.ac.uk" ]; then
#  OOMPHPATH="/home/mly/oomphlib_current"
#  cd $OOMPHPATH/src/ && make && make install && \
#  cd $OOMPHPATH/user_drivers/lagrange_square/ && \
#  make square0 && \
#  ./square0 --w_solver 0 --ns_solver 1 --f_solver 2 --p_solver 1 \
#            --visc 0 --ang 30 --rey 100 --noel 4
#else
#  echo "No such target for $RAYTARGET"
#  read -p "Press [Enter] key to quit."
#  exit
#fi




#./square0 --w_solver 0 --ns_solver 0 --visc 0 \
#          --ang 30 --rey 100 --noel 4


#./square0 --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 3 --visc 0 \
#          --ang 30 --rey 100 --noel 4

#cd /home/mly/v327/src/navier_stokes/
#make && make install \
#&& cd /home/mly/v327/user_drivers/lagrange_square/ \
#&& make square3 && ./square3 --w_solver 0 --ns_solver 0 --visc Sim --ang 30 --rey 100 --noel 8 --diagw --doc_soln

#./square0.sh > square0.dat \
#&& ./square1.sh > square1.dat \
#&& ./square3.sh > square3.dat





