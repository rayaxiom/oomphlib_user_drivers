//LIC//====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.85. June 9, 2008.
//LIC//
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================

#include <sstream>
#include <iomanip>
#include <ios>

// Oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"

// My own header
#include "./../rayheader.h"

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"

//using namespace std;
using namespace oomph;

namespace SquareLagrange
{
  // Domain dimenensions.
  double Y_min = 0.0;
  double Y_max = 1.0;
  double X_min = 0.0;
  double X_max = 1.0;
  double Lx = 1.0;
  double Ly = 1.0;

  // CL - set directly from the commandline.
  // To set from CL - a CL value is set, this is changed depending on that
  // value.

  // Set the defaults.
  unsigned NS_solver = 1; //CL, 0 = SuperLU, 1 - LSC
  unsigned F_solver = 0; //CL, 0 - SuperLU, 1 - AMG
  unsigned P_solver = 0; //CL, 0 - SuperLU, 1 - AMG
  unsigned Vis = 0; //CL, 0 - Simple, 1 - Stress divergence
  double Rey = 100.0; //CL, Reynolds number
  double Ang = 0.0; //CL, Reynolds number
  unsigned Noel = 4; //CL, Number of elements in 1D

  std::string Prob_str = "PROBSTR"; //Set from CL, a unique identifier.
  std::string NS_str = "Nl"; //Set from CL, e - Exact, l - LSC
  std::string F_str = "Fe"; //Set from CL, e - Exact, a - AMG
  std::string P_str = "Pe"; //Set from CL, e - Exact, a - AMG
  std::string Vis_str = "Sim"; //Set from CL, Sim - Simple, Str = Stress Diver.
  std::string Rey_str = "R100"; //Set from CL, Reynolds number
  std::string Noel_str = "N4"; //Set from CL, Number of elements in 1D
                              // null.
  bool Doc_prec = false; // To set from CL
  bool Doc_soln = false; // To set from CL
  bool Print_hypre = true;

  std::string Label = ""; // To be set as the label for this problem. Contains
                          // all the information for this run.
  std::string Soln_dir = ""; // Where to put the solution.
  std::string Doc_prec_dir = ""; // Where to put the solution.

  std::string Itstime_dir = ""; //Set from CL, directory to output the 
                                // iteration counts and timing results.

  // Used to determine if we are using the TrilinosAztecOOSolver solver or not.
  // This cannot be determined by the OOMPH_HAS_TRILINOS ifdef since we may be
  // using OOMPH-LIB's GMRES even if we have Trilinos. This should be set in
  // the problem constuctor as soon as we set the linear_solver_pt() for the
  // problem.
  bool Using_trilinos_solver = false;

  double Rey_start = 0.0;
  double Rey_incre = 50.0;
  double Rey_end = 500.0;

  // Object to store the linear solver iterations and times.
  DocLinearSolverInfo* Doc_linear_solver_info_pt;

  double f_amg_strength = -1.0;
  double f_amg_damping = -1.0;
  int f_amg_coarsening = -1;
  int f_amg_smoother = -1;
  int f_amg_iterations = -1;
  int f_amg_smoother_iterations = -1;
  
  double p_amg_strength = -1.0;
  double p_amg_damping = -1.0;
  int p_amg_coarsening = -1;
  int p_amg_smoother = -1;
  int p_amg_iterations = -1;
  int p_amg_smoother_iterations = -1;
}

namespace oomph
{
//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi
//========================================================================
 template<class ELEMENT>
 class SlopingQuadMesh : public RectangularQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingQuadMesh(const unsigned& nx, const unsigned& ny,
                  const double& lx,  const double& ly, const double& phi ) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
   {
    // Find out how many nodes there are
    unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      double x=nod_pt->x(0);
      double y=nod_pt->x(1);

      // Set new nodal coordinates
      nod_pt->x(0)=x*cos(phi)-y*sin(phi);
      nod_pt->x(1)=x*sin(phi)+y*cos(phi);
     }
   }
 };
} // end of namespace oomph


//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class TiltedCavityProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 TiltedCavityProblem();

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
   // Initialise counters for each newton solve.
   Doc_linear_solver_info_pt->setup_new_time_step();
 }

 /// \short Update after solve is empty
 void actions_after_newton_solve()
 {
 }

 void actions_after_newton_step()
 {
   unsigned iters = 0;
   double preconditioner_setup_time = 0.0;
   double solver_time = 0.0;

   // Get the iteration counts.
#ifdef PARANOID
   IterativeLinearSolver* iterative_solver_pt
     = dynamic_cast<IterativeLinearSolver*>
       (this->linear_solver_pt());
   if(iterative_solver_pt == 0)
   {
     std::ostringstream error_message;
     error_message << "Cannot cast the solver pointer." << std::endl;

     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
   else
   {
     iters = iterative_solver_pt->iterations();
     preconditioner_setup_time 
       = iterative_solver_pt->preconditioner_pt()->setup_time();
   }
#else
   iters = static_cast<IterativeLinearSolver*>
             (this->linear_solver_pt())->iterations();
   preconditioner_setup_time = static_cast<IterativeLinearSolver*>
             (this->linear_solver_pt())->preconditioner_pt()->setup_time();

#endif

   // Get the preconditioner setup time.
   


   // Set the solver time.
   if(SquareLagrange::Using_trilinos_solver)
   {
     TrilinosAztecOOSolver* trilinos_solver_pt 
       = dynamic_cast<TrilinosAztecOOSolver*>(this->linear_solver_pt());
     solver_time = trilinos_solver_pt->linear_solver_solution_time();
   }
   else
   {
     solver_time = linear_solver_pt()->linear_solver_solution_time();
   }

   Doc_linear_solver_info_pt->add_iteration_and_time
     (iters,preconditioner_setup_time,solver_time);
 }

 /// Doc the solution
 void doc_solution();

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

 void create_impenetrable_lagrange_elements(const unsigned &b,
                                            Mesh* const &bulk_mesh_pt,
                                            Mesh* const &surface_mesh_pt);

 void set_inflow_BC(const unsigned &b,
                    Mesh* const &bulk_mesh_pt);
 void set_nonslip_BC(const unsigned &b,
                     Mesh* const &bulk_mesh_pt);
 void set_slip_BC(const unsigned &b, const unsigned &free_direction,
                     Mesh* const &bulk_mesh_pt);

private:

 /// Pointer to the "bulk" mesh
 //SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 //Mesh* Surface_mesh_T_pt;
 //Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;
};



//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem()
{
 // Alias the namespace for convenience
 namespace SL = SquareLagrange;
 
 Doc_linear_solver_info_pt = SL::Doc_linear_solver_info_pt;

 // Assign the boundaries:
 //             2 slip
 //         ----------
 //         |        |
 // 3 Inflow|        |1 P.O.
 //         |        |
 //         ----------
 //             0 non slip
 unsigned if_b=3; // inflow
 unsigned po_b=1; // Parallel Outflow.

 /// Setup the mesh
 // # of elements in x-direction
 unsigned nx=SL::Noel;

 // # of elements in y-direction
 unsigned ny=SL::Noel;

 // Domain length in x-direction
 double lx=SL::Lx;

 // Domain length in y-direction
 double ly=SL::Ly;

 //Bulk_mesh_pt =
 // new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,SL::Ang);
 mesh_pt() = new SimpleRectangularQuadMesh<QTaylorHoodElement<2> >
             (nx,ny,lx,ly); 

 // Create a "surface mesh" that will contain only
 // ImposeParallelOutflowElements in boundary 1
 // The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 //Surface_mesh_P_pt = new Mesh;

 // Create ImposeParallelOutflowElement from all elements that are
 // adjacent to the Neumann boundary.
 //create_parall_outflow_lagrange_elements(po_b,
 //                                        Bulk_mesh_pt,Surface_mesh_P_pt);
 //create_impenetrable_lagrange_elements(po_b,
 //                                        Bulk_mesh_pt,Surface_mesh_P_pt);

 unsigned num_bound = mesh_pt()->nboundary();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
   {
     unsigned num_nod=mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
     {
       // Get node
       Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

     }
   }
 }

 // Which boundary are we dealing with?
 unsigned current_bound = -1;
 
 // The number of nodes on a boundary.
 unsigned num_nod = -1;

 // Inflow is at boundary 3
 current_bound = 3;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
   double x1=nod_pt->x(1);

   double u=(x1-0.0)*(2.0-x1);


   nod_pt->set_value(0,u);
   nod_pt->set_value(1,0);
 }

 // Top boundary is slip.
 current_bound = 2;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   if(!nod_pt->is_on_boundary(3))
   {
     nod_pt->unpin(0);
     nod_pt->pin(1);

     nod_pt->set_value(1,0.0);
   }
 }


 // Outflow is at the right
 current_bound = 1;
 // unpin the x direction, to enforce parallel outflow
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(po_b,inod);
   if(!(nod_pt->is_on_boundary(0)) && !(nod_pt->is_on_boundary(2)) )
     {
      nod_pt->unpin(0);
     }
  }

 // Add the two sub meshes to the problem
 //add_sub_mesh(Bulk_mesh_pt);
 //add_sub_mesh(Surface_mesh_P_pt);
 //add_sub_mesh(Surface_mesh_T_pt);
 
 // Combine all submeshes into a single Mesh
 //build_global_mesh();
 //
 // Inflow
// set_inflow_BC(3,mesh_pt());

 // Bottom boundary
// set_nonslip_BC(0,mesh_pt());
// set_nonslip_BC(2,mesh_pt());

 // Right boundary (Outflow)
// set_slip_BC(1,1,mesh_pt());

 // Top boundary
// set_slip_BC(2,1,mesh_pt());

 
 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = mesh_pt()->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &SL::Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assgn equation numbers
 std::cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;

 // The preconditioner for the fluid block:
 if(SL::NS_solver == 0) // Exact solve.
 {
   ExactBlockPreconditioner<CRDoubleMatrix>* exact_block_prec_pt
     = new ExactBlockPreconditioner<CRDoubleMatrix>;
   exact_block_prec_pt->set_nmesh(1);
   exact_block_prec_pt->set_mesh(0,mesh_pt());
   Prec_pt = exact_block_prec_pt;
 }
 else if(SL::NS_solver == 1) // LSC
 {
   NavierStokesSchurComplementPreconditioner* lsc_prec_pt =
     new NavierStokesSchurComplementPreconditioner(this);

   lsc_prec_pt->set_navier_stokes_mesh(mesh_pt());

   Prec_pt = lsc_prec_pt;

   // F block solve
   // Preconditioner for the F block:
   Preconditioner* f_preconditioner_pt = 0;
   // SL::F_solver == 0 is default, so do nothing.
   if(SL::F_solver == 1)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = new HyprePreconditioner;

     // Cast it so we can fiddle with hypre settings
     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(f_preconditioner_pt);

     Hypre_default_settings::
     set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);
#else
     std::ostringstream error_message;
     error_message << "No Hypre detected...\n"
                   << "Cannot set a hypre preconiditoner for f solve";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
#endif
   }
   // Set the preconditioner in the LSC preconditioner.
   lsc_prec_pt->set_f_preconditioner(f_preconditioner_pt);

   // P block solve
   //SL::P_solver == 0 is default, so do nothing.
   if(SL::P_solver == 1)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;

     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(p_preconditioner_pt);

     Hypre_default_settings::
     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

     lsc_prec_pt->set_p_preconditioner(p_preconditioner_pt);
#else
     std::ostringstream error_message;
     error_message << "No Hypre detected...\n"
                   << "Cannot set a hypre preconiditoner for p solve.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
#endif
   }
 } // if for using LSC as NS prec.
 else
 {
   pause("There is no solver for NS.");
 }

// Build solve and preconditioner
//#ifdef OOMPH_HAS_TRILINOS
// TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
// trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
// Solver_pt = trilinos_solver_pt;
// SL::Using_trilinos_solver = true;
//#else
 Solver_pt = new GMRES<CRDoubleMatrix>;
 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.
 static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
 SL::Using_trilinos_solver = false;
//#endif

 Solver_pt->tolerance() = 1.0e-6;
 this->newton_solver_tolerance() = 1.0e-6;

 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution()
{

  namespace SL = SquareLagrange;
  
  std::ofstream some_file;
  std::stringstream filename;
  filename << SL::Soln_dir<<"/"<<SL::Label<<".dat";

  // Number of plot points
  unsigned npts=5;

  // Output solution
  some_file.open(filename.str().c_str());
  mesh_pt()->output(some_file,npts);
  some_file.close();
}

//============RAYRAY===========
/// RAYRAY - nonslip = pin all directions tp zero
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
set_nonslip_BC(const unsigned &b,
               Mesh* const &bulk_mesh_pt)
{
  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
   
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);
    
    // pin all velocity components and set the value to zero.
    for (unsigned velo_i = 0; velo_i < dim; velo_i++) 
    {
      nod_pt->pin(velo_i);
      nod_pt->set_value(velo_i,0);
    }
   }
}

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
set_slip_BC(const unsigned &b, const unsigned &free_direction,
            Mesh* const &bulk_mesh_pt)
{
  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
   
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);
    
    // pin all velocity components and set the value to zero.
    // Except for the free_direction
    for (unsigned velo_i = 0; velo_i < dim; velo_i++) 
    {
      if(velo_i != free_direction)
       {
        nod_pt->pin(velo_i);
        nod_pt->set_value(velo_i,0);
       }
    }
   }
}

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
set_inflow_BC(const unsigned &b,
              Mesh* const &bulk_mesh_pt)
{

 // Alias the namespace for convenience
 namespace SL = SquareLagrange;

 // Check that the dimension is correct.
#ifdef PARANOID
  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
  if(dim != 2)
   {
     std::ostringstream err_msg;
     err_msg << "Inflow implemented for dim = 2 only." << std::endl;

     throw OomphLibError(err_msg.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
#endif

  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
   
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);

    // Pin both velocity components
    nod_pt->pin(0);
    nod_pt->pin(1);

    // Get the x and y cartesian coordinates.
    double x0 = nod_pt->x(0);
    double x1 = nod_pt->x(1);

    // Tilt x1 back the coordinate so we get the original coordinate.
    double x1_old = x1;
    
    // Now calculate the parabolic inflow at this point
    //double u0_old = (x1_old - SL::Y_min)*(SL::Y_max - x1_old);
    double u0_old = (x1_old - SL::Y_min)*(1.0 - x1_old);

    // Now apply the rotation to u0_old, using rotation matrices.
    // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
    // velocity in the x direction only. There is no velocity
    // in the y direction.
    double u0=u0_old;
    double u1=0;
    
    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,u1); 
   }
}

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_parall_outflow_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   // What is the index of the face of element e along boundary b?
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                          face_index,0);


   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
}

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_impenetrable_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeImpenetrabilityElement<ELEMENT>* flux_element_pt = new
    ImposeImpenetrabilityElement<ELEMENT>(bulk_elem_pt,
                                          face_index);
//   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
//    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
//                                          face_index);

   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
}



int str2int(const std::string &str)
{
  std::stringstream ss(str);
  int n;
  ss >> n;
  return n;
}

unsigned str2unsigned(const std::string &str)
{
  std::stringstream ss(str);
  unsigned n;
  ss >> n;
  return n;
}

double str2double(const std::string &str)
{
  std::stringstream ss(str);
  double n;
  ss >> n;
  return n;
}

//===start_of_main======================================================
/// Driver code
//======================================================================
int main(int argc, char* argv[])
{
 #ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);
 #endif

 // Alias the namespace for convenience.
 namespace SL = SquareLagrange;

 // Set up doc info
 DocLinearSolverInfo doc_linear_solver_info;
 
 SL::Doc_linear_solver_info_pt = &doc_linear_solver_info;

 SL::Soln_dir = "RESLT";

 SL::Doc_prec_dir = "rawdata";

 // Store commandline arguments
 CommandLineArgs::setup(argc,argv);


 // Flag to output the solution.
 CommandLineArgs::specify_command_line_flag("--doc_soln", &SL::Soln_dir);
 // Flag to output the preconditioner, used for debugging.

 CommandLineArgs::specify_command_line_flag("--ns_solver", &SL::NS_solver);
 CommandLineArgs::specify_command_line_flag("--p_solver", &SL::P_solver);
 CommandLineArgs::specify_command_line_flag("--f_solver", &SL::F_solver);
 CommandLineArgs::specify_command_line_flag("--visc", &SL::Vis);
 CommandLineArgs::specify_command_line_flag("--rey", &SL::Rey);
 CommandLineArgs::specify_command_line_flag("--rey_start", &SL::Rey_start);
 CommandLineArgs::specify_command_line_flag("--rey_incre", &SL::Rey_incre);
 CommandLineArgs::specify_command_line_flag("--rey_end", &SL::Rey_end);
 CommandLineArgs::specify_command_line_flag("--noel", &SL::Noel);
 
 // Iteration count and times directory.
 CommandLineArgs::specify_command_line_flag("--itstimedir", &SL::Itstime_dir);

 // NS_F block AMG parameters
 CommandLineArgs::specify_command_line_flag("--f_amg_str", &SL::f_amg_strength);
 CommandLineArgs::specify_command_line_flag("--f_amg_damp", &SL::f_amg_damping);
 CommandLineArgs::specify_command_line_flag("--f_amg_coarse", &SL::f_amg_coarsening);
 CommandLineArgs::specify_command_line_flag("--f_amg_smoo", &SL::f_amg_smoother);
 CommandLineArgs::specify_command_line_flag("--f_amg_iter", &SL::f_amg_iterations);
 CommandLineArgs::specify_command_line_flag("--f_amg_smiter", &SL::f_amg_smoother_iterations);

 // NS_P block AMG parameters
 CommandLineArgs::specify_command_line_flag("--p_amg_str", &SL::p_amg_strength);
 CommandLineArgs::specify_command_line_flag("--p_amg_damp", &SL::p_amg_damping);
 CommandLineArgs::specify_command_line_flag("--p_amg_coarse", &SL::p_amg_coarsening);
 CommandLineArgs::specify_command_line_flag("--p_amg_smoo", &SL::p_amg_smoother);
 CommandLineArgs::specify_command_line_flag("--p_amg_iter", &SL::p_amg_iterations);
 CommandLineArgs::specify_command_line_flag("--p_amg_smiter", &SL::p_amg_smoother_iterations);


 // Parse the above flags.
 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags();

 ////////////////////////////////////////////////////
 // Now set up the flags/parameters for the problem//
 ////////////////////////////////////////////////////
 
 // Document the solution? Default is false.
 if(CommandLineArgs::command_line_flag_has_been_set("--doc_soln"))
 {
   SL::Doc_soln = true;
 }

 // Set a string to identify the problem. This is unique to each problem,
 // so we hard code this. 2DStrPo = 2 dimension, straight parallel outflow.
 // straight describes the velocity flow field. Po = Parallel outflow
 // describes the boundary type.
 SL::Prob_str = "SqVan_semi";

 // Set the strings to identify the preconditioning,
 // This is used purely for book keeping purposes.
 
 // Default: NS_solver = 1, NS_str = Nl
 if(CommandLineArgs::command_line_flag_has_been_set("--ns_solver"))
 {
  switch(SL::NS_solver)
  {
    case 0:
      SL::NS_str = "Ne";
      SL::P_str = "";
      SL::F_str = "";
      break;
    case 1:
      SL::NS_str = "Nl";
      break;
    default:
      std::cout << "Do not recognise NS: " << SL::NS_solver << "\n"
                << "Exact solve = 0\n"
                << "LSC = 1\n"
                << "Using default: LSC for NS block (NS_solver = 1)"<<std::endl;
  }  // switch
 } // if

 // Default: This can only be set if NS_solver != 0 i.e. we are using LSC for 
 // the NS block
 // Default: P_solver = 0, P_str = Pe
 if(CommandLineArgs::command_line_flag_has_been_set("--p_solver"))
 {
  if(SL::NS_solver == 0)
  {
    pause("NS solve is exact. There cannot be a P solver.");
  }

  switch(SL::P_solver)
  {
    case 0:
      SL::P_str = "Pe";
      break;
    case 1:
      SL::P_str = "Pa";
      break;
    default:
      std::cout << "Do not recognise P: " << SL::P_solver << "\n"
                << "Exact preconditioning = 0\n"
                << "AMG = 1\n"
                << "Using default: Exact P solve (P_solver = 0)"<< std::endl;
  }  // switch
 } // if

 // This can only be set if we're using LSC for the NS block.
 // Default: 0, Fe
 if(CommandLineArgs::command_line_flag_has_been_set("--f_solver"))
 {
  if(SL::NS_solver == 0)
  {
    pause("NS solve is exact. There cannot be an F solver.");
  }

  switch(SL::F_solver)
  {
    case 0:
      SL::F_str = "Fe";
      break;
    case 1:
      SL::F_str = "Fa";
      break;
    default:
      std::cout << "Do not recognise F: " << SL::F_solver << "\n"
                << "Exact preconditioning = 0\n"
                << "AMG = xxx Look in the code...\n"
                << "Using default: Exact F solve (F_solver = 0)"<< std::endl;
  }  // switch
 } // if
 
 // Set the viscuous term.
 // Default: 0, Sim
 if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
 {
   if (SL::Vis == 0)
   {
     SL::Vis_str = "Sim";
     NavierStokesEquations<2>::Gamma[0]=0.0;
     NavierStokesEquations<2>::Gamma[1]=0.0;

   }
   else if (SL::Vis == 1)
   {
     SL::Vis_str = "Str";
     NavierStokesEquations<2>::Gamma[0]=1.0;
     NavierStokesEquations<2>::Gamma[1]=1.0;
   } // else - setting viscuous term.
   else
   {
     std::cout << "There is no such Viscous term, using 0 = simple." 
               << std::endl; 
   }
 }

 // Set Noel_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
 {
   std::ostringstream strs;
   strs << "N" << SL::Noel;
   SL::Noel_str = strs.str();
 }

 // Solve with Taylor-Hood element, set up problem
 // Set Rey_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--rey"))
 {
   if(SL::Rey >= 0)
   {
     std::ostringstream strs;
     strs << "R" << SL::Rey;
     SL::Rey_str = strs.str();
   }
   else
   {
     std::cout << "Looping Reynolds" << std::endl; 
     std::cout << "Rey_start: " << SL::Rey_start << std::endl; 
     std::cout << "Rey_incre: " << SL::Rey_incre << std::endl; 
     std::cout << "Rey_end: " << SL::Rey_end << std::endl; 
   }
 }

 TiltedCavityProblem< QTaylorHoodElement<2> > problem;


///////////////////////////////////////////////////////////////////////////////

 if(SL::Rey < 0)
 {
   unsigned rey_increment = 0;

   for (SL::Rey = SL::Rey_start; 
        SL::Rey <= SL::Rey_end; SL::Rey += SL::Rey_incre)
   {
     std::ostringstream strs;
     strs << "R" << SL::Rey;
     SL::Rey_str = strs.str();

     // Setup the label. Used for doc solution and preconditioner.
     SL::Label = SL::Prob_str
                 + SL::NS_str + SL::F_str + SL::P_str
                 + SL::Vis_str + SL::Rey_str
                 + SL::Noel_str;

     time_t rawtime;
     time(&rawtime);

     std::cout << "RAYDOING: "
       << SL::Label
       << " on " << ctime(&rawtime) << std::endl;


     // Solve the problem
     problem.newton_solve();

     //Output solution
     if(SL::Doc_soln){problem.doc_solution();}
  
     //////////////////////////////////////////////////////////////////////////
     ////////////// Outputting results ////////////////////////////////////////
     //////////////////////////////////////////////////////////////////////////
  
  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // my rank and number of processors. This is used later for putting the data.
  unsigned my_rank = comm_pt->my_rank();
  
  // Output the iteration counts and times if my_rank == 0
  if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir") 
     && (my_rank == 0))
  {
   
   // Create the File...
   std::ostringstream filename_stream;
   filename_stream << SL::Itstime_dir<<"/"<<SL::Label;
   std::ofstream outfile;
   outfile.open(filename_stream.str().c_str());

   // We now output the iteration and time.
   Vector<Vector<Vector<double> > > iters_times
     = SL::Doc_linear_solver_info_pt->iterations_and_times();

   // Below outputs the iteration counts and time.
   // Output the number of iterations
   // Since this is a steady state problem, there is only
   // one "time step".
   //*
   // Loop over the time steps and output the iterations, prec setup time and
   // linear solver time.
   //unsigned ntimestep = iters_times.size();
   //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
   {
     // New timestep:
     outfile << "RAYITS:\t" << rey_increment << "\t";
     std::cout << "RAYITS:\t" << rey_increment << "\t";
     
     // Loop through the Newtom Steps
     unsigned nnewtonstep = iters_times[rey_increment].size();
     unsigned sum_of_newtonstep_iters = 0;
     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
         innewtonstep++)
     {
       sum_of_newtonstep_iters += iters_times[rey_increment][innewtonstep][0];
       outfile << iters_times[rey_increment][innewtonstep][0] << " ";
       std::cout << iters_times[rey_increment][innewtonstep][0] << " ";
     }
     double average_its = ((double)sum_of_newtonstep_iters)
       / ((double)nnewtonstep);

     // Print to one decimal place if the average is not an exact
     // integer. Otherwise we print normally.
     std::streamsize cout_precision = std::cout.precision();
     ((unsigned(average_its*10))%10)?
       outfile << "\t"<< std::fixed << std::setprecision(1)
       << average_its << "(" << nnewtonstep << ")" << std::endl:
       outfile << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
     outfile << std::setprecision(cout_precision);

     ((unsigned(average_its*10))%10)?
       std::cout << "\t"<< std::fixed << std::setprecision(1)
       << average_its << "(" << nnewtonstep << ")" << std::endl:
       std::cout << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
     std::cout << std::setprecision(cout_precision);
   }

   // Now doing the preconditioner setup time.
   //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
   {
     // New timestep:
     outfile << "RAYPRECSETUP:\t" << rey_increment << "\t";
     std::cout << "RAYPRECSETUP:\t" << rey_increment << "\t";
     // Loop through the Newtom Steps
     unsigned nnewtonstep = iters_times[rey_increment].size();
     double sum_of_newtonstep_times = 0;
     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
         innewtonstep++)
     {
       sum_of_newtonstep_times += iters_times[rey_increment][innewtonstep][1];
       outfile << iters_times[rey_increment][innewtonstep][1] << " ";
       std::cout << iters_times[rey_increment][innewtonstep][1] << " ";
     }
     double average_time = ((double)sum_of_newtonstep_times)
       / ((double)nnewtonstep);

     // Print to one decimal place if the average is not an exact
     // integer. Otherwise we print normally.
     outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
     std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
   }

   // Now doing the linear solver time.
   //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
   {
     // New timestep:
     outfile << "RAYLINSOLVER:\t" << rey_increment << "\t";
     std::cout << "RAYLINSOLVER:\t" << rey_increment << "\t";
     // Loop through the Newtom Steps
     unsigned nnewtonstep = iters_times[rey_increment].size();
     double sum_of_newtonstep_times = 0;
     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
         innewtonstep++)
     {
       sum_of_newtonstep_times += iters_times[rey_increment][innewtonstep][2];
       outfile << iters_times[rey_increment][innewtonstep][2] << " ";
       std::cout << iters_times[rey_increment][innewtonstep][2] << " ";
     }
     double average_time = ((double)sum_of_newtonstep_times)
       / ((double)nnewtonstep);

     // Print to one decimal place if the average is not an exact
     // integer. Otherwise we print normally.
     outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
     std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
   }
   outfile.close();
  } // if my_rank == 0

  rey_increment++;





   }
 }
 else
 {
   // Setup the label. Used for doc solution and preconditioner.
   SL::Label = SL::Prob_str
               + SL::NS_str + SL::F_str + SL::P_str
               + SL::Vis_str + SL::Rey_str
               + SL::Noel_str;

   time_t rawtime;
   time(&rawtime);

   std::cout << "RAYDOING: "
     << SL::Label
     << " on " << ctime(&rawtime) << std::endl;

   // Solve the problem
   problem.newton_solve();

   //Output solution
   if(SL::Doc_soln){problem.doc_solution();}
   
  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // my rank and number of processors. This is used later for putting the data.
  unsigned my_rank = comm_pt->my_rank();
  
  // Output the iteration counts and times if my_rank == 0
  if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir") 
     && (my_rank == 0))
  {
   
   // Create the File...
   std::ostringstream filename_stream;
   filename_stream << SL::Itstime_dir<<"/"<<SL::Label;
   std::ofstream outfile;
   outfile.open(filename_stream.str().c_str());

   // We now output the iteration and time.
   Vector<Vector<Vector<double> > > iters_times
     = SL::Doc_linear_solver_info_pt->iterations_and_times();

   // Below outputs the iteration counts and time.
   // Output the number of iterations
   // Since this is a steady state problem, there is only
   // one "time step".
   //*
   // Loop over the time steps and output the iterations, prec setup time and
   // linear solver time.
   unsigned ntimestep = iters_times.size();
   for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
   {
     // New timestep:
     outfile << "RAYITS:\t" << intimestep << "\t";
     std::cout << "RAYITS:\t" << intimestep << "\t";
     
     // Loop through the Newtom Steps
     unsigned nnewtonstep = iters_times[intimestep].size();
     unsigned sum_of_newtonstep_iters = 0;
     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
         innewtonstep++)
     {
       sum_of_newtonstep_iters += iters_times[intimestep][innewtonstep][0];
       outfile << iters_times[intimestep][innewtonstep][0] << " ";
       std::cout << iters_times[intimestep][innewtonstep][0] << " ";
     }
     double average_its = ((double)sum_of_newtonstep_iters)
       / ((double)nnewtonstep);

     // Print to one decimal place if the average is not an exact
     // integer. Otherwise we print normally.
     std::streamsize cout_precision = std::cout.precision();
     ((unsigned(average_its*10))%10)?
       outfile << "\t"<< std::fixed << std::setprecision(1)
       << average_its << "(" << nnewtonstep << ")" << std::endl:
       outfile << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
     outfile << std::setprecision(cout_precision);

     ((unsigned(average_its*10))%10)?
       std::cout << "\t"<< std::fixed << std::setprecision(1)
       << average_its << "(" << nnewtonstep << ")" << std::endl:
       std::cout << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
     std::cout << std::setprecision(cout_precision);
   }

   // Now doing the preconditioner setup time.
   for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
   {
     // New timestep:
     outfile << "RAYPRECSETUP:\t" << intimestep << "\t";
     std::cout << "RAYPRECSETUP:\t" << intimestep << "\t";
     // Loop through the Newtom Steps
     unsigned nnewtonstep = iters_times[intimestep].size();
     double sum_of_newtonstep_times = 0;
     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
         innewtonstep++)
     {
       sum_of_newtonstep_times += iters_times[intimestep][innewtonstep][1];
       outfile << iters_times[intimestep][innewtonstep][1] << " ";
       std::cout << iters_times[intimestep][innewtonstep][1] << " ";
     }
     double average_time = ((double)sum_of_newtonstep_times)
       / ((double)nnewtonstep);

     // Print to one decimal place if the average is not an exact
     // integer. Otherwise we print normally.
     outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
     std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
   }

   // Now doing the linear solver time.
   for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
   {
     // New timestep:
     outfile << "RAYLINSOLVER:\t" << intimestep << "\t";
     std::cout << "RAYLINSOLVER:\t" << intimestep << "\t";
     // Loop through the Newtom Steps
     unsigned nnewtonstep = iters_times[intimestep].size();
     double sum_of_newtonstep_times = 0;
     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
         innewtonstep++)
     {
       sum_of_newtonstep_times += iters_times[intimestep][innewtonstep][2];
       outfile << iters_times[intimestep][innewtonstep][2] << " ";
       std::cout << iters_times[intimestep][innewtonstep][2] << " ";
     }
     double average_time = ((double)sum_of_newtonstep_times)
       / ((double)nnewtonstep);

     // Print to one decimal place if the average is not an exact
     // integer. Otherwise we print normally.
     outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
     std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
   }
   outfile.close();
  } // if my_rank == 0
 } // else do not loop reynolds


#ifdef OOMPH_HAS_MPI
// finalize MPI
MPI_Helpers::finalize();
#endif
 return(EXIT_SUCCESS);
} // end_of_main
