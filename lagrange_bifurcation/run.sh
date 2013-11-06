#!/bin/bash

# Use
# tetgen -a0.2 fsi_bifurcation_fluid.poly
# to generate the mesh.
# 0.2 is the maximum area/volumne of the elements.
# decrease it to get a finer mesh.

## PROGRAM NAME
PROGRAM="unstructured_three_d_fluid"
RESLTDIR="RESLT_TH"

touch ./$RESLTDIR && rm -rf ./$RESLTDIR && mkdir $RESLTDIR && \
make $PROGRAM && \
./$PROGRAM --w_solver 0 --ns_solver 0 --rey 200 --doc_soln $RESLTDIR && \
cd $RESLTDIR && \
../plot_it.bash && \
paraview soln.pvd


