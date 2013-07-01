//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
//Driver for 2D rectangular driven cavity

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 double Phi_min=0.0;
 double Phi_max=90.0;

 /// Reynolds number
 double Re=200.0;

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain
//====================================================================
template<class ELEMENT>
class RectangularDrivenCavityProblem : public Problem
{

public:


 /// Constructor
 RectangularDrivenCavityProblem();

 /// Destructor (empty)
 ~RectangularDrivenCavityProblem(){}

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Update the after solve (empty)
 void actions_after_newton_solve(){}


 /// \short Update the problem specs before solve. 
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
 {
  // Setup inflow
  unsigned ibound=3; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    
    // Pointer to node:
    Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
    
    // Get the x/y coordinates
    double x=nod_pt->x(0);
    double y=nod_pt->x(1);
    
    double phi=atan2(y,x);
    const double pi = 4.0*atan(1.0);
    double u=1;

//    double u=-(phi-Global_Physical_Variables::Phi_min*pi/180.0)*
//     (phi-Global_Physical_Variables::Phi_max*pi/180.0);
    nod_pt->set_value(0,u*cos(phi));
    nod_pt->set_value(1,u*sin(phi));
   }
  
 } // end_of_actions_before_newton_solve

 // Access function for the specific mesh
 SimpleRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<SimpleRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }


 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RectangularDrivenCavity problem 
//========================================================================
template<class ELEMENT>
RectangularDrivenCavityProblem<ELEMENT>::RectangularDrivenCavityProblem()
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=30;

 // # of elements in y-direction
 unsigned n_y=30;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);


 // Calculate the value of pi
 const double pi = 4.0*atan(1.0);
 
 const double r_min=1.0;
 const double r_max=3.0;

 // Find out how many nodes there are
 unsigned n_node=mesh_pt()->nnode();
 
 // Loop over all nodes
 for (unsigned n=0;n<n_node;n++)
  {
   // Pointer to node:
   Node* nod_pt=mesh_pt()->node_pt(n);
   
   // Get the x/y coordinates
   double x_old=nod_pt->x(0);
   double y_old=nod_pt->x(1);
   
   // Map from the old x/y to the new r/phi:
   double r=r_min+(r_max-r_min)*x_old;
   double phi=(Global_Physical_Variables::Phi_min+
               (Global_Physical_Variables::Phi_max-
                Global_Physical_Variables::Phi_min)*y_old)*pi/180.0;
   
   // Set new nodal coordinates
   nod_pt->x(0)=r*cos(phi);
   nod_pt->x(1)=r*sin(phi);
  }

 // hierher
 bool add_parallel_outflow=true;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     if (ibound == 3) // inflow boundary
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
      }
     else if(ibound == 0) // bottom boundary
      {
       // pin y, leave x free
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     else if(ibound == 2)
      {
       // pin x, leave y free
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
     else
      {
       if (!add_parallel_outflow)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
        }
      }
    }
  } // end loop over boundaries

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements



 if (add_parallel_outflow)
  {   
   unsigned b=1;
   
   // How many bulk elements are adjacent to boundary b?
   unsigned n_element = mesh_pt()->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      mesh_pt()->boundary_element_pt(b,e));
     
     //What is the index of the face of the element e along boundary b
     int face_index = mesh_pt()->face_index_at_boundary(b,e);
     
     // Build the corresponding lagrange element
     ImposeParallelOutflowElement<ELEMENT>* el_pt = new 
      ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,face_index);
     
     // Add it to the mesh
     mesh_pt()->add_element_pt(el_pt);

     // Loop over the nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
     {
       Node* nod_pt = el_pt->node_pt(j);

       // Is the node also on boundary 0 or 2?
       if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
       {
         // How many nodal values were used by the "bulk" element
         // that originally created this node?
         unsigned n_bulk_value=el_pt->nbulk_value(j);

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

 // Now set the first pressure value in element 0 to 0.0
 //fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RectangularDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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
 mesh_pt()->output(some_file,npts);
 some_file.close();
} // end_of_doc_solution





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for RectangularDrivenCavity test problem -- test drive
/// with two different types of element.
//=====================================================================
int main()
{

 // Set up doc info
 // ---------------

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;

 // ---------------
 // end of Set up doc info

 // Doing QTaylorHoodElements
 {
  
  // Build the problem with QTaylorHoodElements
  RectangularDrivenCavityProblem<QTaylorHoodElement<2> > problem;
  cout << "Doing QTaylorHoodElement<2>" << std::endl;
  
  // Solve the problem
  problem.newton_solve();
  
  // Outpt the solution
  problem.doc_solution(doc_info);

  // Step number
  doc_info.number()++;

 } // end of QTaylorHoodElements


} // end_of_main










