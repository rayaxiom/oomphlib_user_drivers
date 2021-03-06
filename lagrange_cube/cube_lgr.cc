//LIC// ====================================================================
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

#include <fenv.h>

#include <sstream>

// Oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
#include "../rayheader.h"

// The 3D mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"

using namespace oomph;

//===start_of_namespace=================================================
/// Namepspace for global parameters
//======================================================================
//namespace GParam
//{
// /// Current Reynolds number
// double re = 0.0;
// 
// double ang_x = 0.0;
// double ang_y = 0.0;
// double ang_z = 0.0;
// 
// /// Lagrange Multiplier ID 
// const unsigned lagrange_multiplier_id = 42;
// 
// /// Constant for pi
// const double pi = 4.0*atan(1.0);
//
//
// /// Traction at the outflow boundary
// void prescribed_traction(const double& t,
//                          const Vector<double>& x,
//                          Vector<double>& traction)
// {
//  traction.resize(3);
//  traction[0]=1.0;
//  traction[1]=0.0;
//  traction[2]=0.0;
// } 
//} // end of namespace


namespace oomph {


////=========================================================================
///// Wrapper class 
////========================================================================
// template <class ELEMENT>
// class BulkFluidSubjectToLagrangeMultiplierElement : public virtual ELEMENT
// {
// 
// public:
// 
//  /// Default constructor
//  BulkFluidSubjectToLagrangeMultiplierElement() : ELEMENT() {}
// 
//  /// \short Returns the number of DOF types associated with this element:
//  /// Twice the number of its spatial dimension plus one (for the pressure)
//  unsigned ndof_types()
//   {
//    return 2*ELEMENT::dim()+1;
//   }
// 
//  /// \short Create a list of pairs for all unknowns in this element,
//  /// so that the first entry in each pair contains the global equation
//  /// number of the unknown, while the second one contains the number
//  /// of the "DOF" that this unknown is associated with.
//  /// (Function can obviously only be called if the equation numbering
//  /// scheme has been set up.)\n
//  /// E.g. in a 2D problem there are 5 types of DOF:\n
//  /// 0 - x velocity (without lagr mult )\n
//  /// 1 - y velocity (without lagr mult )\n
//  /// 2 - pressure
//  /// 3 - x velocity (with lagr mult )\n
//  /// 4 - y velocity (with lagr mult )\n
//  void get_dof_numbers_for_unknowns(
//   std::list<std::pair<unsigned long,unsigned> >& block_lookup_list)
//   {
//    // temporary pair (used to store block lookup prior to being added to list
//    std::pair<unsigned,unsigned> block_lookup;
//   
//    // number of nodes
//    const unsigned n_node = this->nnode();
//   
//    //Get the dimension of the node
//    const unsigned nodal_dim = ELEMENT::nodal_dimension();
//    //cout << "nodal_dim: " << nodal_dim << endl;
//   
//    //Integer storage for local unknown
//    int local_unknown=0;
//   
//    //Loop over the nodes
//    for(unsigned n=0;n<n_node;n++)
//     {
//      // Check if node has been resized, i.e. if Lagrange multipliers
//      // have been attached; in that case the veloc dofs get
//      // a different dof type.
//      unsigned offset = 0;
//      if (this->node_pt(n)->nvalue() != this->required_nvalue(n))
//       {
//        offset = ELEMENT::dim()+1;
//       }
//     
//      //Loop over dimension for velocity components
//      for(unsigned i=0;i<nodal_dim;i++)
//       {
//        //If the variable is free
//        local_unknown = ELEMENT::nodal_local_eqn(n,i);
//       
//        // ignore pinned values
//        if (local_unknown >= 0)
//         {
//          // store block lookup in temporary pair: First entry in pair
//          // is global equation number; second entry is block type
//          block_lookup.first = this->eqn_number(local_unknown);
//          block_lookup.second = offset+i;
//         
//          // add to list
//          block_lookup_list.push_front(block_lookup);  
//         
//         }
//       
//        // Pressure (Taylor Hood only!)
//        if (this->required_nvalue(n)==(ELEMENT::dim()+1))
//         {         
//          //If the variable is free
//          local_unknown = ELEMENT::nodal_local_eqn(n,ELEMENT::dim());
//         
//          // ignore pinned values
//          if (local_unknown >= 0)
//           {
//            // store block lookup in temporary pair: First entry in pair
//            // is global equation number; second entry is block type
//            block_lookup.first = this->eqn_number(local_unknown);
//            block_lookup.second = ELEMENT::dim();
//           
//            // add to list
//            block_lookup_list.push_front(block_lookup);
//           
//           }
//         }
//       }
//     }
//   }
// };
   
////===========start_face_geometry=============================================
///// FaceGeometry of wrapped element is the same as the underlying element
////===========================================================================
// template<class ELEMENT>
// class FaceGeometry<BulkFluidSubjectToLagrangeMultiplierElement<ELEMENT> > :
//  public virtual FaceGeometry<ELEMENT>
// {
// };



//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from SimpleCubicMesh:
/// the same mesh rotated with an angle phi
//========================================================================
 template<class ELEMENT> 
 class SlopingCubicMesh : public SimpleCubicMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                   const double& lx,  const double& ly, const double& lz,
                   const double& ang_x, 
                   const double& ang_y, 
                   const double& ang_z) :
   SimpleCubicMesh<ELEMENT>(nx,ny,nz,lx,ly,lz)
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
      double z=nod_pt->x(2);

      // Set new nodal coordinates by premultiplying by R_xyz.
      nod_pt->x(0)=x*cos(ang_y)*cos(ang_z)
       -y*cos(ang_y)*sin(ang_z)
       +z*sin(ang_y);
      nod_pt->x(1)=x*(cos(ang_x)*sin(ang_z) + cos(ang_z)*sin(ang_x)*sin(ang_y))
       +y*(cos(ang_x)*cos(ang_z) - sin(ang_x)*sin(ang_y)*sin(ang_z))
       -z*(cos(ang_y)*sin(ang_x));
      nod_pt->x(2)=x*(sin(ang_x)*sin(ang_z) - cos(ang_x)*cos(ang_z)*sin(ang_y))
       +y*(cos(ang_z)*sin(ang_x) + cos(ang_x)*sin(ang_y)*sin(ang_z))
       +z*(cos(ang_x)*cos(ang_y));
     }
   }
 };
}



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class DrivenCavityThreeDProblem : public Problem
{

public:

 /// Constructor
 DrivenCavityThreeDProblem();


 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end of fix_pressure

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by 
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b, 
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);
                                  
 /// Create traction elements on outflow boundary
 void create_traction_elements();

 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

private:

 void set_prec_and_solver();


 // Number the boundaries.
 unsigned Front_bound;
 unsigned Back_bound;
 unsigned Left_bound;
 unsigned Right_bound;
 unsigned Top_bound;
 unsigned Bottom_bound;

 /// ID of imposed flow boundary
 unsigned Imposed_flow_boundary;

 /// ID of neumann boundary
 unsigned Neumann_boundary;

 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh for Lagrange multiplier elements
 Mesh* Surface_mesh_pt1;
 
 /// Pointer to the "surface" mesh for Lagrange multiplier elements
 Mesh* Surface_mesh_pt2;

 // Preconditioner
 Preconditioner* Prec_pt;

 // Solver
 IterativeLinearSolver* Solver_pt;

 // Object to hold the iteration counts.
 DocLinearSolverInfo* Doc_linear_solver_info_pt;
};



//==start_of_constructor==================================================
/// Constructor for DrivenCavity problem 
// 
//                                               z
//                                     4 Back    |
//      ____________                             |
//     /|          /|                            |_______y
//    / |         / |          5 Top            /
//   /  |        /  |                          /
//  /   |       /   |                        x/
// /____|______/    |     1 Left    3 Right
// |    |______|____|      
// |   /       |   /
// |  /        |  /            0 Bottom
// | /         | /
// |/__________|/        2 Front
//  
//
// Front_bound = 2;
// Back_bound = 4;
// Left_bound = 1;
// Right_bound = 3;
// Top_bound = 5;
// Bottom_bound = 0;
//
// Imposed_flow_boundary=4;
// Neumann_boundary=2;
// 
//========================================================================
template<class ELEMENT> 
DrivenCavityThreeDProblem<ELEMENT>::DrivenCavityThreeDProblem()
{ 
 namespace CL = CubeLagrange;

 // Set the boundaries so we do not have to use numbers!
 Front_bound = 2;
 Back_bound = 4;
 Left_bound = 1;
 Right_bound = 3;
 Top_bound = 5;
 Bottom_bound = 0;

 // Set the boundaries, you see, it's much easier to read with the variable
 // names from above.
 Imposed_flow_boundary = Back_bound;
 Neumann_boundary = Front_bound;

 // If the pointer to the documentation object is not null (i.e. it has been
 // set), then we set it.
 if(CL::Doc_linear_solver_info_pt != 0)
 {
   Doc_linear_solver_info_pt = CL::Doc_linear_solver_info_pt;
 }

 // Build and assign mesh
 Bulk_mesh_pt = 
  new SlopingCubicMesh<ELEMENT >(CL::Noelx, CL::Noely, CL::Noelz,
                                 CL::Lx, CL::Ly, CL::Lz,
                                 CL::Angx,CL::Angy,CL::Angz);

 // Create "surface mesh" that will contain only the Lagrange multiplier 
 // elements.
 Surface_mesh_pt1 = new Mesh;

 create_parall_outflow_lagrange_elements(Neumann_boundary,
                                         Bulk_mesh_pt,
                                         Surface_mesh_pt1);
 
 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt1);
// add_sub_mesh(Surface_mesh_pt2);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default. 
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u, v and w velocities)
     unsigned spatial_dim = 3;
     for (unsigned i=0;i<spatial_dim;i++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  }

 
 // Impose inflow and negative inflow (outflow).
  unsigned num_nod= Bulk_mesh_pt->nboundary_node(Imposed_flow_boundary);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Imposed_flow_boundary,inod);
    double x=nod_pt->x(0);
    double y=nod_pt->x(1);
    double z=nod_pt->x(2);
      
    // Reverse the tilting
    // We apply R_zyx
    double tilt_back_y = x*(cos(-CL::Angy)*sin(-CL::Angz))
     + y*(cos(-CL::Angx)*cos(-CL::Angz) 
          + sin(-CL::Angx)*sin(-CL::Angy)*sin(-CL::Angz))
     + z*(cos(-CL::Angx)*sin(-CL::Angy)*sin(-CL::Angz) 
          - cos(-CL::Angz)*sin(-CL::Angx));
    double tilt_back_z = -x*sin(-CL::Angy)
     +y*cos(-CL::Angy)*sin(-CL::Angx)
     +z*cos(-CL::Angx)*cos(-CL::Angy);
      
    // The imposed velocity along the x boundary.
    double u = 0.0;
    u=(tilt_back_y)*(1.0-tilt_back_y)*(tilt_back_z)*(1.0-tilt_back_z);

//    if (tilt_back_y>0.5)
//     {
//      u=(tilt_back_y-0.5)*(1.0-tilt_back_y)*(tilt_back_z)*(1.0-tilt_back_z);
//     }
//    else
//     {
//      u=-(tilt_back_y)*(0.5-tilt_back_y)*(tilt_back_z)*(1.0-tilt_back_z);
//     }
      
    // Now apply Rxyz to u, using rotation matrices.
    // We have velocity in the x direction only.
    // Thus the vector to rotate is [u,0,0] since the imposed
    // velocity in the y and z direction is 0.
    double u_x =u*cos(CL::Angy)*cos(CL::Angz);
    double u_y =u*(cos(CL::Angx)*sin(CL::Angz) 
                   + cos(CL::Angz)*sin(CL::Angx)*sin(CL::Angy));
    double u_z =u*(sin(CL::Angx)*sin(CL::Angz) 
                   - cos(CL::Angx)*cos(CL::Angz)*sin(CL::Angy));
      
    nod_pt->set_value(0,u_x);
    nod_pt->set_value(1,u_y);
    nod_pt->set_value(2,u_z);
   }
   
 // Unpin the nodes constrained by the Lagrange multiplier. 
  num_nod= Bulk_mesh_pt->nboundary_node(Neumann_boundary);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Neumann_boundary,inod);
    // Only free if node is on a single boundary
    std::set<unsigned>* bnd_pt=0;
    nod_pt->get_boundaries_pt(bnd_pt);
    if (bnd_pt!=0)
     {
      if (bnd_pt->size()<2)
       {
        nod_pt->unpin(0);
        nod_pt->unpin(1); // Why was only 0 unpinned?
        nod_pt->unpin(2); // Why was only 0 unpinned?
       }
     }
   }

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = 
    dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &CL::Rey;
  }
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

 set_prec_and_solver();
 
}

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void DrivenCavityThreeDProblem<ELEMENT>::set_prec_and_solver()
{
  // Alias the namespace for convenience
  namespace CL = CubeLagrange;

  ////// Build the preconditioner
  LagrangeEnforcedflowPreconditioner* prec_pt
    = new LagrangeEnforcedflowPreconditioner;

  Prec_pt = prec_pt;

  Vector<Mesh*> mesh_pt;
  mesh_pt.resize(2);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_pt1;
  //meshes_pt[2] = Surface_mesh_T_pt;
  
  prec_pt->set_meshes(mesh_pt);

  if(!CL::Use_axnorm)
  {
    prec_pt->scaling_sigma() = CL::Scaling_sigma;
  }

 // W solver. SuperLU is default so we do nothing here.
 if(CL::W_solver == 0)
 {
 }
 else
 {
   std::cout << "Other W solvers not complemented yet. Using default SuperLU"
             << std::endl;
 }

 // The preconditioner for the fluid block:
 if(CL::NS_solver == 0) // Exact solve.
 {}
 else if(CL::NS_solver == 1) // LSC
 {
   NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
     new NavierStokesSchurComplementPreconditioner(this);

   prec_pt->set_navier_stokes_lsc_preconditioner(ns_preconditioner_pt);
   ns_preconditioner_pt->set_navier_stokes_mesh(Bulk_mesh_pt);

   // F block solve
   // Preconditioner for the F block:
   Preconditioner* f_preconditioner_pt = 0;
   // CL::F_solver == 0 is default, so do nothing.
   if(CL::F_solver == 11)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_2D_poison_problem();
#endif
   }
   else if(CL::F_solver == 12)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_navier_stokes_momentum_block();
#endif
   }
   else if(CL::F_solver == 13)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPGSStrn075();
#endif
   }
   else if(CL::F_solver == 14)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSGSStrn075();
#endif
   }
   else if(CL::F_solver == 15)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPPilutStrn075();
#endif
   }
   else if(CL::F_solver == 16)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSPilutStrn075();
#endif
   }
   else if(CL::F_solver == 17)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_augmented_momentum_block();
#endif
   }
   else if(CL::F_solver == 81)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPGSStrn0668();
#endif
   }
   else if(CL::F_solver == 82)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPJStrn0668();
#endif
   }
   else if(CL::F_solver == 83)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPPilutStrn0668();
#endif
   }
   else if(CL::F_solver == 84)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSGSStrn0668();
#endif
   }
   else if(CL::F_solver == 85)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSJStrn0668();
#endif
   }
   else if(CL::F_solver == 86)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSPilutStrn0668();
#endif
   }
   else if(CL::F_solver == 2)
   {
//     f_preconditioner_pt = new RayBlockDiagonalPreconditioner<CRDoubleMatrix>;
     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
   }
   else if(CL::F_solver == 3)
   {
     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
#ifdef OOMPH_HAS_HYPRE
     dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
       (f_preconditioner_pt)->set_subsidiary_preconditioner_function
       (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_for_2D_poison_problem);
#endif
   }
   else if (CL::F_solver == 69)
   {
#ifdef OOMPH_HAS_HYPRE
     // AMG coarsening: Ruge-Stuben
     RayGlobalAMGParam::amg_coarsening = 1;
     
     // AMG smoother: Gauss-Seidel
     RayGlobalAMGParam::amg_smoother=0;
     
     // There is no damping with GS, otherwise we set the parameter:
     // RayGlobalAMGParam::amg_damping

     // Different amg strength for simple/stress divergence for viscuous term.
     if(CL::Vis == 0)
     {
       // Simple form
       RayGlobalAMGParam::amg_strength = 0.25;
     }
     else
     {
       // Stress divergence form
       RayGlobalAMGParam::amg_strength = 0.668;
     }
     
     // Setup the preconditioner.
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_using_2D_poisson_base();
#endif
   }
   else if (CL::F_solver == 96)
   {
#ifdef OOMPH_HAS_HYPRE
     // AMG coarsening:
     // Set: RayGlobalAMGParam::amg_coarsening = 
     // 0 - CLJP
     // 1 - RS
     
     // AMG smoother:
     // Set: RayGlobalAMGParam::amg_smoother = 
     // 0 - Jacobi (Need to set damping as well)
     // 1 - Gauss-Seidel
     // 2 - Pilut
     
     // There is no damping with GS, otherwise we set the parameter:
     // RayGlobalAMGParam::amg_damping

     RayGlobalAMGParam::amg_strength = CL::f_amg_strength;
     RayGlobalAMGParam::amg_damping = CL::f_amg_damping;
     RayGlobalAMGParam::amg_coarsening = CL::f_amg_coarsening;
     RayGlobalAMGParam::amg_smoother = CL::f_amg_smoother;
     RayGlobalAMGParam::amg_iterations = CL::f_amg_iterations;
     RayGlobalAMGParam::amg_smoother_iterations = CL::f_amg_smoother_iterations;
     RayGlobalAMGParam::print_hypre = CL::Print_hypre;


     // Different amg strength for simple/stress divergence for viscuous term.
     if(RayGlobalAMGParam::amg_strength < 0.0)
     {
       if(CL::Vis == 0)
       {
         // Simple form
         RayGlobalAMGParam::amg_strength = 0.25;
       }
       else
       {
         // Stress divergence form
         RayGlobalAMGParam::amg_strength = 0.668;
       }
     }
     
     // Setup the preconditioner.
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_ray();
#endif
   }

   // Set the preconditioner in the LSC preconditioner.
   ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
   
   // P block solve
   //CL::P_solver == 0 is default, so do nothing.
   if(CL::P_solver == 1)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;

     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(p_preconditioner_pt);

     Hypre_default_settings::
     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
   }
   else if(CL::P_solver == 96)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = 0;
     
//* 
     RayGlobalAMGParam::amg_iterations = CL::p_amg_iterations;
     RayGlobalAMGParam::amg_smoother_iterations = CL::p_amg_smoother_iterations;
     RayGlobalAMGParam::amg_smoother = CL::p_amg_smoother;
     RayGlobalAMGParam::amg_strength = CL::p_amg_strength;
   //RayGlobalAMGParam::amg_damping = CL::p_amg_damping;
     RayGlobalAMGParam::amg_coarsening = CL::p_amg_coarsening;
     RayGlobalAMGParam::print_hypre = CL::Print_hypre;

//     std::cout << "p_amg_iterations:" << CL::p_amg_iterations << std::endl; 
//     std::cout << "p_amg_smoother_iterations" << CL::p_amg_smoother_iterations << std::endl; 
//     std::cout << "p_amg_strength" << CL::p_amg_strength << std::endl;
//     std::cout << "p_amg_coarsening" << CL::p_amg_coarsening << std::endl; 
// */

     p_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::set_hypre_ray();

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
   }
   else if(CL::P_solver == 2)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;

     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(p_preconditioner_pt);

     hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

     // Setup v-cycles
     hypre_preconditioner_pt->set_amg_iterations(2);
     hypre_preconditioner_pt->amg_smoother_iterations() = 2;
     
     // Setup smoother
     // simple: 0 - DJ, 1 - GS
     // compelx: Pilut - 7
     hypre_preconditioner_pt->amg_using_simple_smoothing();
     hypre_preconditioner_pt->amg_simple_smoother() = 0;
     // only applicable for DJ
     hypre_preconditioner_pt->amg_damping() = 0.8;

     // Setup coarsening
     // 0 - CLJP
     // 1 - RS
     hypre_preconditioner_pt->amg_coarsening() = 1;

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif

   }
 } // if for using LSC as NS prec.
 else
 {
   pause("There is no solver for NS.");
 }


 // Set the doc info for book keeping purposes.
 prec_pt->set_doc_linear_solver_info_pt(CL::Doc_linear_solver_info_pt);

 if(CL::Use_block_diagonal_w)
 {
   prec_pt->use_block_diagonal_w_block();
 }
 else
 {
   prec_pt->use_diagonal_w_block();
 }

 if(CL::Doc_prec)
 {
   prec_pt->enable_doc_prec();
 }

 // Set the label, use to output information from the preconditioner, such
 // as the block matrices and the rhs vector
 prec_pt->set_label_pt(&CL::Label);
 prec_pt->set_doc_prec_directory_pt(&CL::Doc_prec_dir);

 // Build solve and preconditioner
#ifdef OOMPH_HAS_TRILINOS
 TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
 trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
 Solver_pt = trilinos_solver_pt;
 CL::Using_trilinos_solver = true;
#else
 Solver_pt = new GMRES<CRDoubleMatrix>;
 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.
 static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
 CL::Using_trilinos_solver = false;
#endif

 Solver_pt->tolerance() = 1.0e-6;
 this->newton_solver_tolerance() = 1.0e-6;

 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;
}



////============start_of_create_traction_elements==========================
///// Create Navier-Stokes traction elements on outflow boundary
////=======================================================================
//template<class ELEMENT>
//void DrivenCavityThreeDProblem<ELEMENT>::create_traction_elements()
//{
//
// unsigned b=Neumann_boundary;
//
// // How many bulk elements are adjacent to boundary b?
// unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
//
// // Loop over the bulk elements adjacent to boundary b?
// for(unsigned e=0;e<n_element;e++)
//  {
//   // Get pointer to the bulk element that is adjacent to boundary b
//   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
//    Bulk_mesh_pt->boundary_element_pt(b,e));
//   
//   //What is the index of the face of element e along boundary b
//   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
//   
//   // Build the corresponding prescribed-flux element
//   NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
//    NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
//
//   //Add the prescribed-flux element to the surface mesh
//   Surface_mesh_pt1->add_element_pt(flux_element_pt);
//   
//   // Set the pointer to the prescribed traction function
//   flux_element_pt->traction_fct_pt() = &GParam::prescribed_traction;
//   
//  }
//
//}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void DrivenCavityThreeDProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 some_file.open("pressure_dofs.dat");
 unsigned nnod=Bulk_mesh_pt->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   if (Bulk_mesh_pt->node_pt(j)->nvalue()==3)
    {
     some_file << Bulk_mesh_pt->node_pt(j)->x(0) << " " 
               << Bulk_mesh_pt->node_pt(j)->x(1) << " " 
               << Bulk_mesh_pt->node_pt(j)->eqn_number(2) << " " 
               << std::endl;
    }
  }
 some_file.close();
}



//============start_of_create_parall_outflow_lagrange_elements===========
/// Create my_lagrange_element  on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the 
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void DrivenCavityThreeDProblem<ELEMENT>::
create_parall_outflow_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 namespace CL = CubeLagrange;
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
   
   // Build the corresponding my_lagrange_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new 
    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                          face_index,
                                          CL::Lagrange_multiplier_id);
                                          
   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

 
   // Loop over the nodes 
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);
     
     // Is the node also on boundary 0, 1, 3 or 5?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(1))
         || (nod_pt->is_on_boundary(3))||(nod_pt->is_on_boundary(5)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       // Cast to a boundary node
       BoundaryNodeBase *bnod_pt =
        dynamic_cast<BoundaryNodeBase*>(nod_pt);
       
       // Get the index of the first Lagrange multiplier
       unsigned first_lmi=bnod_pt->
        index_of_first_value_assigned_by_face_element(
         CL::Lagrange_multiplier_id);
       
       // There is only one lagrange multiplier.
       for (unsigned j=0;j<2;j++)
        {
         nod_pt->pin(first_lmi+j);
        }
      }
    }
  }
} // FoFunc create_parall_outflow_lagrange_elements(...)



//===start_of_main======================================================
/// Driver code 
//======================================================================
int main(int argc, char* argv[]) 
{
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif


 namespace CL = CubeLagrange;
 
 DocLinearSolverInfo doc_linear_solver_info;

 CL::Doc_linear_solver_info_pt = &doc_linear_solver_info;


 // Store commandline arguments
 CommandLineArgs::setup(argc,argv);

 CommandLineArgs::specify_command_line_flag("--dist_prob");

 // Flag to output the solution.
 CommandLineArgs::specify_command_line_flag("--doc_soln", &CL::Soln_dir);
 // Flag to output the preconditioner, used for debugging.
 CommandLineArgs::specify_command_line_flag("--doc_prec", &CL::Doc_prec_dir);
 
 // A problem ID, there are eight different types of problems.
 // Check the header file.
 CommandLineArgs::specify_command_line_flag("--prob_id",&CL::Prob_id);

 CommandLineArgs::specify_command_line_flag("--w_solver", &CL::W_solver);
 CommandLineArgs::specify_command_line_flag("--ns_solver", &CL::NS_solver);
 CommandLineArgs::specify_command_line_flag("--p_solver", &CL::P_solver);
 CommandLineArgs::specify_command_line_flag("--f_solver", &CL::F_solver);
 CommandLineArgs::specify_command_line_flag("--visc", &CL::Vis);
 CommandLineArgs::specify_command_line_flag("--ang", &CL::Ang_deg);
 CommandLineArgs::specify_command_line_flag("--angx", &CL::Angx_deg);
 CommandLineArgs::specify_command_line_flag("--angy", &CL::Angy_deg);
 CommandLineArgs::specify_command_line_flag("--angz", &CL::Angz_deg);
 CommandLineArgs::specify_command_line_flag("--rey", &CL::Rey);
 CommandLineArgs::specify_command_line_flag("--rey_start", &CL::Rey_start);
 CommandLineArgs::specify_command_line_flag("--rey_incre", &CL::Rey_incre);
 CommandLineArgs::specify_command_line_flag("--rey_end", &CL::Rey_end);
 CommandLineArgs::specify_command_line_flag("--noel", &CL::Noel);
 CommandLineArgs::specify_command_line_flag("--noelx", &CL::Noelx);
 CommandLineArgs::specify_command_line_flag("--noely", &CL::Noely);
 CommandLineArgs::specify_command_line_flag("--noelz", &CL::Noelz);
 CommandLineArgs::specify_command_line_flag("--sigma",
                                            &CL::Scaling_sigma);
 CommandLineArgs::specify_command_line_flag("--bdw");
 
 // Iteration count and times directory.
 CommandLineArgs::specify_command_line_flag("--itstimedir", &CL::Itstime_dir);

 // NS_F block AMG parameters
 CommandLineArgs::specify_command_line_flag("--f_amg_str", &CL::f_amg_strength);
 CommandLineArgs::specify_command_line_flag("--f_amg_damp", &CL::f_amg_damping);
 CommandLineArgs::specify_command_line_flag("--f_amg_coarse", &CL::f_amg_coarsening);
 CommandLineArgs::specify_command_line_flag("--f_amg_smoo", &CL::f_amg_smoother);
 CommandLineArgs::specify_command_line_flag("--f_amg_iter", &CL::f_amg_iterations);
 CommandLineArgs::specify_command_line_flag("--f_amg_smiter", &CL::f_amg_smoother_iterations);

 // NS_P block AMG parameters
 CommandLineArgs::specify_command_line_flag("--p_amg_str", &CL::p_amg_strength);
 CommandLineArgs::specify_command_line_flag("--p_amg_damp", &CL::p_amg_damping);
 CommandLineArgs::specify_command_line_flag("--p_amg_coarse", &CL::p_amg_coarsening);
 CommandLineArgs::specify_command_line_flag("--p_amg_smoo", &CL::p_amg_smoother);
 CommandLineArgs::specify_command_line_flag("--p_amg_iter", &CL::p_amg_iterations);
 CommandLineArgs::specify_command_line_flag("--p_amg_smiter", &CL::p_amg_smoother_iterations);


 // Parse the above flags.
 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags(); 


 if(CommandLineArgs::command_line_flag_has_been_set("--dist_prob"))
 {
   CL::Distribute_problem = true;
 }
 else
 {
   CL::Distribute_problem = false;
 }


 // Document the solution? Default is false.
 if(CommandLineArgs::command_line_flag_has_been_set("--doc_soln"))
 {
   // The argument immediately after --doc_soln is put into CL::Soln_dir.
   // If this begins with "--", then no solution directory has been provided.
   std::size_t found = CL::Soln_dir.find("--");
   
   // Check if they have set the solution directory.
   if(found != std::string::npos)
   {
     std::ostringstream err_msg;
     err_msg << "Please provide the doc_soln directory "
             << "after the argument --doc_soln.\n" 
             << "This must not start with \"--\"." << std::endl;

     throw OomphLibError(err_msg.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
   else
   {
     CL::Doc_soln = true;
   }
 }

 // Document the preconditioner? Default is false.
 if(CommandLineArgs::command_line_flag_has_been_set("--doc_prec"))
 {
   // The argument immediately after --doc_prec is put into CL::Doc_prec_dir.
   // If this begins with "--", then no prec directory has been provided.
   std::size_t found = CL::Doc_prec_dir.find("--");

   // Check if they have set the doc_prec directory.
   if(found != std::string::npos)
   {
     std::ostringstream err_msg;
     err_msg << "Please provide the doc_prec directory "
             << "after the argument --doc_prec.\n" 
             << "This must not start with \"--\"." << std::endl;

     throw OomphLibError(err_msg.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
   else
   {
     CL::Doc_prec = true;
   }
 }


 // Set a problem id to identify the problem.
 // This is used for book keeping purposes.
 if(CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
 {
   // The argument immediately after --prob_id is put into CL::Prob_id.
   // If this begins with "--", then no problem id has been provided.

   // Maybe I should check if CL::Prob_id is a number or a string...

   // We only accept problem IDs as defined below.
   // Creating a set of acceptable IDs
//   int prob_id_array[]= {10,11,12,13,
//                         20,21,22,23};
   int prob_id_array[]= {10,11};

   bool inset = check_if_in_set<int>(prob_id_array,2,CL::Prob_id);

   // Check if they have provided an acceptable ID.
   // If a new element has been inserted, it means the user has provided an
   // ID not in the set.
   if(inset == false)
   {
     std::ostringstream err_msg;
     err_msg << "Please provide a problem id to identify the problem after "
             << "after the argument --prob_id.\n" 
             << "Acceptable IDs are:\n"
             << "10 = (SqTmp) Square, custom stuff...\n"
             << "11 = (SqPo) Square, Parallel outflow (para inflow)\n"
             << "12 = (SqTf) Square, Tangential flow (Semi para inflow)\n"
             << "13 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)\n"
             << "\n"
             << "20 = (AwTmp) Annulus wedge, custom stuff...\n"
             << "21 = (AwPo) Annulus wedge, Parallel outflow (para inflow)\n"
             << "22 = (AwTf) Annulus wedge, Tangential flow (semi para inflow)\n"
             << "23 = (AwTfPo) Annulus wedge, Tan. flow, Para. outflow (semi para inflow)\n"
             << std::endl;

     throw OomphLibError(err_msg.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
 }

 // Set the viscuous term.
 // Default: 0, Sim
 if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
 {
   if (CL::Vis == 0)
   {
     NavierStokesEquations<3>::Gamma[0]=0.0;
     NavierStokesEquations<3>::Gamma[1]=0.0;
     NavierStokesEquations<3>::Gamma[2]=0.0;
   }
   else if (CL::Vis == 1)
   {
     NavierStokesEquations<3>::Gamma[0]=1.0;
     NavierStokesEquations<3>::Gamma[1]=1.0;
     NavierStokesEquations<3>::Gamma[2]=1.0;
   } // else - setting viscuous term.
   else
   {
     std::ostringstream err_msg;
     err_msg << "Do not recognise viscuous term: " << CL::Vis << ".\n"
             << "Vis = 0 for simple form\n"
             << "Vis = 1 for stress divergence form\n"
             << std::endl;
     throw OomphLibError(err_msg.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
 }

 if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
 {
   CL::Angx = CL::Ang_deg * (MathematicalConstants::Pi / 180.0);
   CL::Angy = CL::Angx;
   CL::Angz = CL::Angx;
 }
 else if(CommandLineArgs::command_line_flag_has_been_set("--angx")
         && CommandLineArgs::command_line_flag_has_been_set("--angy")
         && CommandLineArgs::command_line_flag_has_been_set("--angz"))
 {
   CL::Angx = CL::Angx_deg * (MathematicalConstants::Pi / 180.0);
   CL::Angy = CL::Angy_deg * (MathematicalConstants::Pi / 180.0);
   CL::Angz = CL::Angz_deg * (MathematicalConstants::Pi / 180.0);
 }
 else
 {
   std::ostringstream err_msg;
   err_msg << "Please supply either --ang or all three --angx --angy --angz\n"
           << std::endl;
   throw OomphLibError(err_msg.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
 }


 // Check if the Reynolds numbers have been set.
 if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
    &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
    &&CommandLineArgs::command_line_flag_has_been_set("--rey_end")
    &&CommandLineArgs::command_line_flag_has_been_set("--rey"))
 {
   std::ostringstream err_msg;
   err_msg << "You have set all --rey* argument, please choose carefully!\n"
           << std::endl;
   throw OomphLibError(err_msg.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
 }
 else if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
    &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
    &&CommandLineArgs::command_line_flag_has_been_set("--rey_end"))
 {
   std::cout << "Looping Reynolds: \n"
             << "Rey_start = " << CL::Rey_start << std::endl; 
   std::cout << "Rey_incre = " << CL::Rey_incre << std::endl; 
   std::cout << "Rey_end = " << CL::Rey_end << std::endl; 
 }
 else if(!CommandLineArgs::command_line_flag_has_been_set("--rey"))
 {
   std::ostringstream err_msg;
   err_msg << "No Reynolds numbers have been set.\n"
           << "For a single Reynolds number, use --rey.\n"
           << "For looping through Reynolds numbers, use:\n"
           << "--rey_start --rey_incre --rey_end.\n"
           << std::endl;
   throw OomphLibError(err_msg.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
 }

 if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
 {
   // This will use the copy constructor.
   CL::Noelx = CL::Noel;
   CL::Noely = CL::Noel;
   CL::Noelz = CL::Noel;
 }

 CL::Label = CL::create_label();

 time_t rawtime;
 time(&rawtime);

 std::cout << "RAYDOING: "
           << CL::Label
           << " on " << ctime(&rawtime) << std::endl;
    
 //problem.distribute();

 DrivenCavityThreeDProblem<QTaylorHoodElement<3> > problem;

 // Set up doc info
 DocInfo doc_info;
 doc_info.number()=0;
 doc_info.set_directory("RESLT");
 
 // Solve the problem 
 problem.newton_solve();
           
 // Doc solution
 problem.doc_solution(doc_info);
 doc_info.number()++;

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
} // end of main
