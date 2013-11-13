#!/bin/bash

## PROGRAM NAME
PROGRAM="unstructured_three_d_fluid"
AUTOGENSCRIPT="autogen_wipebuild_noselftest.sh"

## Directories
RESLTDIR="RESLT"
DOCPRECDIR="prec_data"
CURRENTDIR=`pwd`
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


CMD_INPUT="--w_solver 0 --ns_solver 1 --visc 0 --rey 0 --mesh_type 0 --mesh_area 0.0125 --ntsteps 4 --doc_prec $DOCPRECDIR"

# Setup directories
function setup_dir
{
  touch ./$RESLTDIR && rm -rf ./$RESLTDIR && mkdir $RESLTDIR && \
  touch ./$DOCPRECDIR && rm -rf ./$DOCPRECDIR && mkdir $DOCPRECDIR
}

function make_src
{
  cd $OOMPH_ROOT_DIR/src && \
  make && make install && \
  cd $CURRENTDIR

#  cd $OOMPH_ROOT_DIR && \
#  ./$AUTOGENSCRIPT --jobs=4 --rebuild && \
#  cd $CURRENTDIR
}

function make_run_program
{
  make $PROGRAM && \
  mpirun ./$PROGRAM $CMD_INPUT
}

###############################################################################
###############################################################################
#setup_dir && make_src && make_run_program
make_run_program



###############################################################################
###############################################################################

#make $PROGRAM && \
#./$PROGRAM --w_solver 0 --ns_solver 0 --rey 100 --doc_soln $RESLTDIR && \
#cd $RESLTDIR && \
#../plot_it.bash && \
#paraview soln.pvd


