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

PROGRAM="sq_lgr"

###############################################################################
#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 2 ./$PROGRAM $RUNARGUMENTS


#cd $OOMPH_ROOT_DIR/src && make && make install && \
#cd $CURRENTDIR && \
#make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

touch $VALIDATEDIR
rm -rf $VALIDATEDIR
mkdir $VALIDATEDIR
cd $VALIDATEDIR
mkdir $TEMPVALIDATADIR

RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS

RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS




RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS


RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS

# Now I need to compare files... grep for RAYITS
files=("SqPoWeNeFePeSimA42R100N16WdNP1R0" \
"SqPoWeNeFePeSimA42R100N16WdNP2R0" \
"SqPoWeNeFePeSimA42R100N16WdNP2R1" \
"SqPoWeNeFePeSimA42R100N16WdNP3R0" \
"SqPoWeNeFePeSimA42R100N16WdNP3R1" \
"SqPoWeNeFePeSimA42R100N16WdNP3R2" \
"SqPoWeNeFePeSimA42R100N16WdNP4R0" \
"SqPoWeNeFePeSimA42R100N16WdNP4R1" \
"SqPoWeNeFePeSimA42R100N16WdNP4R2" \
"SqPoWeNeFePeSimA42R100N16WdNP4R3" \
"SqPoWeNeFePeSimA42R100N32WdNP1R0" \
"SqPoWeNeFePeSimA42R100N32WdNP2R0" \
"SqPoWeNeFePeSimA42R100N32WdNP2R1" \
"SqPoWeNeFePeSimA42R100N32WdNP3R0" \
"SqPoWeNeFePeSimA42R100N32WdNP3R1" \
"SqPoWeNeFePeSimA42R100N32WdNP3R2" 
"SqPoWeNeFePeSimA42R100N32WdNP4R0" \
"SqPoWeNeFePeSimA42R100N32WdNP4R1" \
"SqPoWeNeFePeSimA42R100N32WdNP4R2" \
"SqPoWeNeFePeSimA42R100N32WdNP4R3" \
"SqPoWeNeFePeSimA42R100N8WdNP1R0" \
"SqPoWeNeFePeSimA42R100N8WdNP2R0" \
"SqPoWeNeFePeSimA42R100N8WdNP2R1" \
"SqPoWeNeFePeSimA42R100N8WdNP3R0" \
"SqPoWeNeFePeSimA42R100N8WdNP3R1" \
"SqPoWeNeFePeSimA42R100N8WdNP3R2" \
"SqPoWeNeFePeSimA42R100N8WdNP4R0" \
"SqPoWeNeFePeSimA42R100N8WdNP4R1" \
"SqPoWeNeFePeSimA42R100N8WdNP4R2" \
"SqPoWeNeFePeSimA42R100N8WdNP4R3" \
"SqPoWeNeFePeStrA42R100N16WdNP1R0" \
"SqPoWeNeFePeStrA42R100N16WdNP2R0" \
"SqPoWeNeFePeStrA42R100N16WdNP2R1" \
"SqPoWeNeFePeStrA42R100N16WdNP3R0" \
"SqPoWeNeFePeStrA42R100N16WdNP3R1" \
"SqPoWeNeFePeStrA42R100N16WdNP3R2" \
"SqPoWeNeFePeStrA42R100N16WdNP4R0" \
"SqPoWeNeFePeStrA42R100N16WdNP4R1" \
"SqPoWeNeFePeStrA42R100N16WdNP4R2" \
"SqPoWeNeFePeStrA42R100N16WdNP4R3" \
"SqPoWeNeFePeStrA42R100N32WdNP1R0" \
"SqPoWeNeFePeStrA42R100N32WdNP2R0" \
"SqPoWeNeFePeStrA42R100N32WdNP2R1" \
"SqPoWeNeFePeStrA42R100N32WdNP3R0" \
"SqPoWeNeFePeStrA42R100N32WdNP3R1" \
"SqPoWeNeFePeStrA42R100N32WdNP3R2" \
"SqPoWeNeFePeStrA42R100N32WdNP4R0" \
"SqPoWeNeFePeStrA42R100N32WdNP4R1" \
"SqPoWeNeFePeStrA42R100N32WdNP4R2" \
"SqPoWeNeFePeStrA42R100N32WdNP4R3" \
"SqPoWeNeFePeStrA42R100N8WdNP1R0" \
"SqPoWeNeFePeStrA42R100N8WdNP2R0" \
"SqPoWeNeFePeStrA42R100N8WdNP2R1" \
"SqPoWeNeFePeStrA42R100N8WdNP3R0" \
"SqPoWeNeFePeStrA42R100N8WdNP3R1" \
"SqPoWeNeFePeStrA42R100N8WdNP3R2" \
"SqPoWeNeFePeStrA42R100N8WdNP4R0" \
"SqPoWeNeFePeStrA42R100N8WdNP4R1" \
"SqPoWeNeFePeStrA42R100N8WdNP4R2" \
"SqPoWeNeFePeStrA42R100N8WdNP4R3" \
"SqPoWeNlFePeSimA42R100N16WdNP1R0" \
"SqPoWeNlFePeSimA42R100N16WdNP2R0" \
"SqPoWeNlFePeSimA42R100N16WdNP2R1" \
"SqPoWeNlFePeSimA42R100N16WdNP3R0" \
"SqPoWeNlFePeSimA42R100N16WdNP3R1" \
"SqPoWeNlFePeSimA42R100N16WdNP3R2" \
"SqPoWeNlFePeSimA42R100N16WdNP4R0" \
"SqPoWeNlFePeSimA42R100N16WdNP4R1" \
"SqPoWeNlFePeSimA42R100N16WdNP4R2" \
"SqPoWeNlFePeSimA42R100N16WdNP4R3" \
"SqPoWeNlFePeSimA42R100N32WdNP1R0" \
"SqPoWeNlFePeSimA42R100N32WdNP2R0" \
"SqPoWeNlFePeSimA42R100N32WdNP2R1" \
"SqPoWeNlFePeSimA42R100N32WdNP3R0" \
"SqPoWeNlFePeSimA42R100N32WdNP3R1" \
"SqPoWeNlFePeSimA42R100N32WdNP3R2" \
"SqPoWeNlFePeSimA42R100N32WdNP4R0" \
"SqPoWeNlFePeSimA42R100N32WdNP4R1" \
"SqPoWeNlFePeSimA42R100N32WdNP4R2" \
"SqPoWeNlFePeSimA42R100N32WdNP4R3" \
"SqPoWeNlFePeSimA42R100N8WdNP1R0" \
"SqPoWeNlFePeSimA42R100N8WdNP2R0" \
"SqPoWeNlFePeSimA42R100N8WdNP2R1" \
"SqPoWeNlFePeSimA42R100N8WdNP3R0" \
"SqPoWeNlFePeSimA42R100N8WdNP3R1" \
"SqPoWeNlFePeSimA42R100N8WdNP3R2" \
"SqPoWeNlFePeSimA42R100N8WdNP4R0" \
"SqPoWeNlFePeSimA42R100N8WdNP4R1" \
"SqPoWeNlFePeSimA42R100N8WdNP4R2" \
"SqPoWeNlFePeSimA42R100N8WdNP4R3" \
"SqPoWeNlFePeStrA42R100N16WdNP1R0" \
"SqPoWeNlFePeStrA42R100N16WdNP2R0" \
"SqPoWeNlFePeStrA42R100N16WdNP2R1" \
"SqPoWeNlFePeStrA42R100N16WdNP3R0" \
"SqPoWeNlFePeStrA42R100N16WdNP3R1" \
"SqPoWeNlFePeStrA42R100N16WdNP3R2" \
"SqPoWeNlFePeStrA42R100N16WdNP4R0" \
"SqPoWeNlFePeStrA42R100N16WdNP4R1" \
"SqPoWeNlFePeStrA42R100N16WdNP4R2" \
"SqPoWeNlFePeStrA42R100N16WdNP4R3" \
"SqPoWeNlFePeStrA42R100N32WdNP1R0" \
"SqPoWeNlFePeStrA42R100N32WdNP2R0" \
"SqPoWeNlFePeStrA42R100N32WdNP2R1" \
"SqPoWeNlFePeStrA42R100N32WdNP3R0" \
"SqPoWeNlFePeStrA42R100N32WdNP3R1" \
"SqPoWeNlFePeStrA42R100N32WdNP3R2" \
"SqPoWeNlFePeStrA42R100N32WdNP4R0" \
"SqPoWeNlFePeStrA42R100N32WdNP4R1" \
"SqPoWeNlFePeStrA42R100N32WdNP4R2" \
"SqPoWeNlFePeStrA42R100N32WdNP4R3" \
"SqPoWeNlFePeStrA42R100N8WdNP1R0" \
"SqPoWeNlFePeStrA42R100N8WdNP2R0" \
"SqPoWeNlFePeStrA42R100N8WdNP2R1" \
"SqPoWeNlFePeStrA42R100N8WdNP3R0" \
"SqPoWeNlFePeStrA42R100N8WdNP3R1" \
"SqPoWeNlFePeStrA42R100N8WdNP3R2" \
"SqPoWeNlFePeStrA42R100N8WdNP4R0" \
"SqPoWeNlFePeStrA42R100N8WdNP4R1" \
"SqPoWeNlFePeStrA42R100N8WdNP4R2" \
"SqPoWeNlFePeStrA42R100N8WdNP4R3" )

touch validation.log

for i in "${files[@]}"
do
	grep "RAYITS" $TEMPVALIDATADIR/$i > RAYITS_new
  grep "RAYITS" ./../$VALIDATADIR/$i > RAYITS_old
  
  DIFF=$(diff RAYITS_new RAYITS_old)
  if [ "$DIFF" != "" ]
  then
    echo "File not the same: $i" >> validation.log
  fi
done

cd $CURRENTDIR

#make $PROGRAM && mpirun -np 2 ./$PROGRAM $RUNARGUMENTS

###############################################################################

#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 1 ./$PROGRAM --ns_solver 1 --visc 0 --ang 0 --rey 0 --noel 64 --itstimedir ray_temp

# I should get:
# RAYITS: 0    19 34    26.5(2)



