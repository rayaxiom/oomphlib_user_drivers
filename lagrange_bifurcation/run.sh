#!/bin/bash

touch ./RESLT_TH && rm -rf ./RESLT_TH && mkdir RESLT_TH && \
make unstructured_three_d_fluid && \
./unstructured_three_d_fluid && \
cd RESLT_TH && \
./../plot_it.bash && \
paraview fluid_soln.pvd


