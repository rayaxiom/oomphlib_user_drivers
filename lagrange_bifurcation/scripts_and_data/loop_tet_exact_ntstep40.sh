#!/bin/bash

#egrep "RAYDOING|RAYAV" loop_tet_exact.dat

tetgen -a0.4 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 0 --mesh_type 0 --ntsteps 40

tetgen -a0.4 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 100 --mesh_type 0 --ntsteps 40

tetgen -a0.4 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 0 --rey 200 --mesh_type 0 --ntsteps 40

###############################################################################


tetgen -a0.4 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 0 --mesh_type 0 --ntsteps 40
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 0 --mesh_type 0 --ntsteps 40

tetgen -a0.4 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 100 --mesh_type 0 --ntsteps 40
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 100 --mesh_type 0 --ntsteps 40

tetgen -a0.4 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.2 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.1 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.05 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.025 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 200 --mesh_type 0 --ntsteps 40
tetgen -a0.0125 fsi_bifurcation_fluid.poly
./unstructured_three_d_fluid --w_solver 0 --ns_solver 0 --visc 1 --rey 200 --mesh_type 0 --ntsteps 40


