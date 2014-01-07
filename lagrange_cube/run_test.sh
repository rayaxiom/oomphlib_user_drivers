#!/bin/bash

###############################################################################
# Current directory
CURRENTDIR=`pwd`

# Where the static validata is.
VALIDATADIR="validata"

# This is where the validation is performed. 
# This will be removed at the beginning of every validation.
VALIDATEDIR="Validate"
TEMPVALIDATADIR="temp_validata"



# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
###############################################################################
## EDIT THIS
############

PROGRAM="cube_lgr"

###############################################################################
#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 2 ./$PROGRAM $RUNARGUMENTS


#cd $OOMPH_ROOT_DIR/src && make && make install && \
#cd $CURRENTDIR && \
#make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

#touch $VALIDATEDIR
#rm -rf $VALIDATEDIR
#mkdir $VALIDATEDIR
#cd $VALIDATEDIR
#mkdir $TEMPVALIDATADIR

#RUNARGUMENTS="--dist_prob --doc_soln soln_dir --doc_prec prec_dir --prob_id 11 --visc 1 --angx 0 --angy 0 --angz 0 --rey 100 --noel 4"
RUNARGUMENTS="--dist_prob --doc_soln soln_dir --prob_id 11 --visc 1 --angx 30 --angy 30 --angz 30 --rey 100 --noel 4" #--bdw

PRECARGUMENTS="--w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 "

FINALRUNCOMMAND="mpirun -np 1 $PROGRAM $PRECARGUMENTS $RUNARGUMENTS"
echo "Doing $FINALRUNCOMMAND"

make $PROGRAM && $FINALRUNCOMMAND



