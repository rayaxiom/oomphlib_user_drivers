#!/bin/bash

PROGRAM="unstructured_three_d_fluid"

# Current directory
CURRENTDIR=`pwd`

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

tetgen -a0.2 fsi_bifurcation_fluid.poly && \
make $PROGRAM && \
mpirun ./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0



