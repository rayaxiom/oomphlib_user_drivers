#!/bin/bash

# Use
# tetgen -a0.2 fsi_bifurcation_fluid.poly
# to generate the mesh.
# 0.2 is the maximum area/volumne of the elements.
# decrease it to get a finer mesh.

touch ./RESLT_TH && rm -rf ./RESLT_TH && mkdir RESLT_TH && \
make unstructured_three_d_fluid && \
./unstructured_three_d_fluid && \
cd RESLT_TH && \
./../plot_it.bash && \
paraview fluid_soln.pvd


