#!/bin/bash

#egrep "RAYDOING|RAYAV" loop_ns1.dat

tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.00625 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.003125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100
tetgen -a0.0015625 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69 --rey 100


