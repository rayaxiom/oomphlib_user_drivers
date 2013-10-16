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
//
// This is a header for shared data I keep.
//  - Ray
//
#ifndef RAY_LAGRANGE_SQUARE_HEADER
#define RAY_LAGRANGE_SQUARE_HEADER

// Oomph-lib includes
#include "generic.h"

//using namespace std;
using namespace oomph;

namespace RayGlobalAMGParam
{
  double amg_strength = -1.0;
  double amg_damping = -1.0;
  int amg_coarsening = -1;
  int amg_smoother = -1;

  int amg_iterations = -1;
  int amg_smoother_iterations = -1;
}

#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{                 
  void print_hypre_settings(Preconditioner* preconditioner_pt)
  {
    HyprePreconditioner* h_prec_pt = 
      static_cast<HyprePreconditioner*>(preconditioner_pt);
    
    // Print AMG iterations:
    std::cout << "HyprePreconditioner settings are: " << std::endl;
    std::cout << "Max_iter: " << h_prec_pt->amg_iterations() << std::endl;
    std::cout << "smoother iter: " << h_prec_pt->amg_smoother_iterations() << std::endl;
    std::cout << "Hypre_method: " << h_prec_pt->hypre_method() << std::endl;
    std::cout << "internal_preconditioner: " << h_prec_pt->internal_preconditioner() << std::endl;
    std::cout << "AMG_using_simple_smoothing: " << h_prec_pt->amg_using_simple_smoothing_flag() << std::endl; 
    std::cout << "AMG_simple_smoother: " << h_prec_pt->amg_simple_smoother() << std::endl; 
    std::cout << "AMG_complex_smoother: " << h_prec_pt->amg_complex_smoother() << std::endl; 
    std::cout << "AMG_coarsening: " << h_prec_pt->amg_coarsening() << std::endl;
    std::cout << "AMG_max_levels: " << h_prec_pt->amg_max_levels() << std::endl;
    std::cout << "AMG_damping: " << h_prec_pt->amg_damping() << std::endl;
    std::cout << "AMG_strength: " << h_prec_pt->amg_strength() << std::endl;;
    std::cout << "AMG_max_row_sum: " << h_prec_pt->amg_max_row_sum() << std::endl;
    std::cout << "AMG_truncation: " << h_prec_pt->amg_truncation() << std::endl;
    std::cout << "\n" << std::endl;
  }

  void reset_amg_param()
   {
     RayGlobalAMGParam::amg_strength = -1.0;
     RayGlobalAMGParam::amg_damping = -1.0;
     RayGlobalAMGParam::amg_coarsening = -1;
     RayGlobalAMGParam::amg_smoother = -1;
     RayGlobalAMGParam::amg_iterations = -1;
     RayGlobalAMGParam::amg_smoother_iterations = -1;
   }

  Preconditioner* set_hypre_for_2D_poison_problem()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    if(!(RayGlobalAMGParam::amg_strength < 0))
    {
      hypre_preconditioner_pt->amg_strength() = RayGlobalAMGParam::amg_strength;
      std::cout << "RAYAMGSTR: " << hypre_preconditioner_pt->amg_strength() << std::endl; 
    }

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_navier_stokes_momentum_block()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);

    if(!(RayGlobalAMGParam::amg_strength < 0))
    {
      hypre_preconditioner_pt->amg_strength() = RayGlobalAMGParam::amg_strength;
      std::cout << "RAYAMGSTR: " << hypre_preconditioner_pt->amg_strength() << std::endl; 
    }

    if(!(RayGlobalAMGParam::amg_damping < 0))
    {
      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;
      std::cout << "RAYDAMP: " << hypre_preconditioner_pt->amg_damping() << std::endl;
    }
 
    return another_preconditioner_pt;
  } 
  ///////////////////////////////////////////////
  ///////////////////////////////////////////////

  Preconditioner* set_hypre_ray(const bool& print_hypre=false)
  {
    // Everything can be set! But defaults are used...

    // Create the preconditioner and cast it.
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

   ////////// SETTING STUFF //////////

   // AMG preconditioner
   hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

    // Set the amg_iterations.
    if(RayGlobalAMGParam::amg_iterations == -1)
     {
       RayGlobalAMGParam::amg_iterations = 1;
     }

    hypre_preconditioner_pt
      ->set_amg_iterations(RayGlobalAMGParam::amg_iterations);
    
    // Set the number of amg smoother iterations.
    if(RayGlobalAMGParam::amg_smoother_iterations == -1)
     {
       RayGlobalAMGParam::amg_smoother_iterations = 2;
     }

    hypre_preconditioner_pt
      ->amg_smoother_iterations() = RayGlobalAMGParam::amg_smoother_iterations;

    // Store this information.
    std::stringstream cyclestream;
    cyclestream << "its"<<RayGlobalAMGParam::amg_iterations
                << "smits" << RayGlobalAMGParam::amg_smoother_iterations;
    
    string cycle_str = cyclestream.str();
    string smoother_str = "";
    string damping_str = "";
    string coarsening_str = "";
    string strength_str = "";

    ///////////////////////////////////////////////////////////////////////////
    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    // Setting the smoother:
    if(RayGlobalAMGParam::amg_smoother == 1)
    {
      smoother_str = "GS";
      // Setting up Gauss-Seidel
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 1;
    }
    else if(RayGlobalAMGParam::amg_smoother == 0)
    {
      smoother_str = "J";
      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;
      
      // Setting up Jacobi with damping.
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 0;
      if(!(RayGlobalAMGParam::amg_damping < 0))
      {
        std::ostringstream strs;
        strs << "Dmp" << hypre_preconditioner_pt->amg_damping();
        damping_str = strs.str(); 
      }
      else
      {
        std::cout << "Please set your damping using --amg_damping" << std::endl; 
        pause("Please do not continue."); 
      }
    }
    else if(RayGlobalAMGParam::amg_smoother == 2)
    {
      smoother_str = "Pilut";
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother() = 7;
    }
    else
    {
      std::cout << "You supplied smoother: " << RayGlobalAMGParam::amg_smoother << std::endl;
      std::cout << "No such smoother. 0 is Jacobi, 1 is GS, 2 is Pilut" << std::endl;
      std::cout << "Please set your smoother using --smoother" << std::endl;
      pause("Please do not continue.");
    }


    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    if(RayGlobalAMGParam::amg_coarsening == 0)
    {
      coarsening_str = "CLJP";
      hypre_preconditioner_pt->amg_coarsening() = 0;
    }
    else if(RayGlobalAMGParam::amg_coarsening == 1)
    {
      coarsening_str = "RS";
      hypre_preconditioner_pt->amg_coarsening() = 1;
    }
    else
    {
      std::cout << "There is no such coarsening: " << RayGlobalAMGParam::amg_coarsening << std::endl;
      std::cout << "0 - CLJP, 1 - RS, use --amg_coarsening" << std::endl;
      pause("Do not continue"); 
      
    }
    
    if(!(RayGlobalAMGParam::amg_strength < 0))
    {
      hypre_preconditioner_pt->amg_strength() = RayGlobalAMGParam::amg_strength;
      std::ostringstream strs;
      strs << "Strn" << hypre_preconditioner_pt->amg_strength();
      strength_str = strs.str(); 
    }
    else
    {
      std::cout << "Please set the amg_strengh using --amg_strength" << std::endl;
      pause("Do not continue");
    }

    if(print_hypre)
     {
      print_hypre_settings(another_preconditioner_pt);
     }
    
    std::cout << "RAYHYPRE: " << cycle_str
                              << coarsening_str 
                              << smoother_str
                              << damping_str
                              << strength_str
                              << std::endl;
    reset_amg_param();
    
    return another_preconditioner_pt;
  }

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  Preconditioner* set_hypre_using_2D_poisson_base()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    string smoother_str = "";
    string damping_str = "";
    string coarsening_str = "";
    string strength_str = "";
    ///////////////////////////////////////////////////////////////////////////
    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    // Setting the smoother:
    if(RayGlobalAMGParam::amg_smoother == 0)
    {
      smoother_str = "GS";
      // Setting up Gauss-Seidel
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 1;
    }
    else if(RayGlobalAMGParam::amg_smoother == 1)
    {
      smoother_str = "J";
      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;
      
      // Setting up Jacobi with damping.
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 0;
      if(!(RayGlobalAMGParam::amg_damping < 0))
      {
        std::ostringstream strs;
        strs << "Dmp" << hypre_preconditioner_pt->amg_damping();
        damping_str = strs.str(); 
      }
      else
      {
        std::cout << "Please set your damping using --amg_damping" << std::endl; 
        pause("Please do not continue."); 
      }
    }
    else if(RayGlobalAMGParam::amg_smoother == 2)
    {
      smoother_str = "Pilut";
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother() = 7;
    }
    else
    {
      std::cout << "You supplied smoother: " << RayGlobalAMGParam::amg_smoother << std::endl;
      std::cout << "No such smoother. 0 is GS, 1 is Jacobi, 2 is Pilut" << std::endl;
      std::cout << "Please set your smoother using --smoother" << std::endl;
      pause("Please do not continue.");
    }


    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    if(RayGlobalAMGParam::amg_coarsening == 0)
    {
      coarsening_str = "CLJP";
      hypre_preconditioner_pt->amg_coarsening() = 0;
    }
    else if(RayGlobalAMGParam::amg_coarsening == 1)
    {
      coarsening_str = "RS";
      hypre_preconditioner_pt->amg_coarsening() = 1;
    }
    else
    {
      std::cout << "There is no such coarsening: " << RayGlobalAMGParam::amg_coarsening << std::endl;
      std::cout << "0 - CLJP, 1 - RS, use --amg_coarsening" << std::endl;
      pause("Do not continue"); 
      
    }
    
    if(!(RayGlobalAMGParam::amg_strength < 0))
    {
      hypre_preconditioner_pt->amg_strength() = RayGlobalAMGParam::amg_strength;
      std::ostringstream strs;
      strs << "Strn" << hypre_preconditioner_pt->amg_strength();
      strength_str = strs.str(); 
    }
    else
    {
      std::cout << "Please set the amg_strengh using --amg_strength" << std::endl;
      pause("Do not continue");
    }
    
    std::cout << "RAYHYPRE: " << coarsening_str 
                              << smoother_str
                              << damping_str
                              << strength_str
                              << std::endl;
    
    return another_preconditioner_pt;
  }
  /////////////////////////////////////////////////////////////////////////////

  Preconditioner* set_hypre_for_augmented_momentum_block()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);

    hypre_preconditioner_pt->amg_strength() = 0.668;

    hypre_preconditioner_pt->amg_damping() = 1.0;
 
    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_CLJPGSStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_RSGSStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_CLJPPilutStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_RSPilutStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  }
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
  Preconditioner* set_hypre_for_CLJPGSStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;
    //hypre_preconditioner_pt->amg_damping() = 1.0;
   
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_CLJPJStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 0;
    hypre_preconditioner_pt->amg_damping() = 1.0;
   
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_CLJPPilutStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 0;
    //hypre_preconditioner_pt->amg_damping() = 1.0;
   
    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
  Preconditioner* set_hypre_for_RSGSStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;
    //hypre_preconditioner_pt->amg_damping() = 1.0;
   
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_RSJStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 0;
    hypre_preconditioner_pt->amg_damping() = 1.0;
   
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_RSPilutStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 0;
    //hypre_preconditioner_pt->amg_damping() = 1.0;
   
    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
} 
#endif            







#endif

