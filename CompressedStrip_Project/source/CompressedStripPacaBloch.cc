/* ---------------------------------------------------------------------
 *
 *
 *
 *
 *
 * ---------------------------------------------------------------------
 *
 *
 * Author: Andrew Akerson
 */

#ifndef COMPRESSEDSTRIPPACABLOCH_CC_
#define COMPRESSEDSTRIPPACABLOCH_CC_
#include "CompressedStripPacaBloch.h"

#include <fstream>
#include <iostream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "mapping_q_eulerian_hp.h"


#define MU_VALUE 1.0
#define NU_VALUE 0.33

#define DIM 2



namespace compressed_strip
{
  using namespace dealii;


  /****************************************************************
                       Function Definitions
  ****************************************************************/

  double NuFunction::value (const Point<DIM>  &p,
                                 const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0, ExcNotImplemented());
    Assert (p.dimension == 2, ExcNotImplemented())

    // Put your function for nu. p(0) is x1 value, p(1) is the x2 value

    return nuValue;
  }

  void NuFunction::value_list(const std::vector< Point< DIM > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = NuFunction::value(points[i], component);

  }

  double MuFunction::value (const Point<DIM>  &p,
                                 const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0, ExcNotImplemented());

    // function for mu. p(0) is x1 value, p(1) is the x2 value
    double muValue;
    if (pieceConstFlag == true)
    {
      muValue = (p(1) >= L1 ? kappa : mu0 );
    }
    else
      muValue = mu0*exp(kappa*p(1));

    return muValue;
  }

  void MuFunction::value_list(const std::vector< Point< DIM > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = MuFunction::value(points[i], component);
  }

  double FE_solution::value (const Point<DIM>  &p,
                             const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0 || component == 1, ExcNotImplemented());
    Assert (p.dimension == 2, ExcNotImplemented())

    // Put your function for nu. p(0) is x1 value, p(1) is the x2 value
    Vector<double> value(2);
    VectorTools::point_value(*dof_handler_ptr, *present_solution, p, value);

    return value(component);
  }



  // computes right hand side values if we were to have body forces. But it just
  // always returns zeros because we don't.
  void right_hand_side (const std::vector<Point<DIM> > &points,
                        std::vector<Tensor<1, DIM> >   &values)
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));
    Assert (DIM >= 2, ExcNotImplemented());


    // not imposing body forces or tractions
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      values[point_n][0] = 0.0;
      values[point_n][1] = 0.0;
    }
  }


  // Functions that get the tensors used in the stiffness and rhs calculations.

  inline
  Tensor<2,DIM>
  ElasticProblem::get_deformation_gradient(std::vector<Tensor<1,DIM> > old_solution_gradient)
  {
    Tensor<2,DIM> tmp = F0;

    for (unsigned int i = 0; i < DIM; i ++)
      for(unsigned int j = 0; j < DIM; j++)
      {
        tmp[i][j] += old_solution_gradient[i][j];
      }

    return tmp;
  }


  ElasticProblem::ElasticProblem ()
    :
    dof_handler (triangulation)
  {

  }




  ElasticProblem::~ElasticProblem ()
  {
    dof_handler.clear ();
    delete mu;
    if(FE_solution_function != NULL)
      delete FE_solution_function;
  }


  void ElasticProblem::create_mesh()
  {

    // creates our strip.
    Point<DIM> corner1, corner2;
    corner1(0) = -domain_dimensions[0]/2.0;
    corner1(1) =  0.0;
    corner2(0) =  domain_dimensions[0]/2.0;
    corner2(1) =  domain_dimensions[1];
    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);

    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell1 = dof_handler.begin_active(),
    endc1 = dof_handler.end();
    for (; cell1!=endc1; ++cell1)
    {
      cell1->set_active_fe_index(section_FE_id[0]);
    }


    // Now we will refine this mesh
    const int numSections = number_sections - 1;
    for (int i = 0 ; i < numSections; i ++)
    {
      double section_x2 = (0.5/numSections)*i + 0.5;
      typename hp::DoFHandler<DIM>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (; cell!=endc; ++cell)
      {
        for (unsigned int v=0;
             v < GeometryInfo<2>::vertices_per_cell;
             ++v)
          {
            const double x2_pos = (cell->vertex(v))(1);
            if ( x2_pos > section_x2)
              {
//                  cell->set_active_fe_index(cell->active_fe_index() + 1);
                cell->set_refine_flag ();
                break;
              }
          }
      }
      triangulation.execute_coarsening_and_refinement ();

      typename hp::DoFHandler<DIM>::active_cell_iterator
      cell2 = dof_handler.begin_active(),
      endc2 = dof_handler.end();
      for (; cell2!=endc2; ++cell2)
      {
        for (unsigned int v=0;
             v < GeometryInfo<2>::vertices_per_cell;
             ++v)
          {
            const double x2_pos = (cell2->vertex(v))(1);
            if ( x2_pos > section_x2)
              {
                cell2->set_active_fe_index(section_FE_id[i + 1]);
//                  cell->set_refine_flag ();
                break;
              }
          }
      }
    }


    // Make sure to renumber the boundaries
    renumber_boundary_ids();

  }

  void ElasticProblem::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.

    // only reinit the present solution if it wasn't already loaded in
    // Also, only distribute the dofs if it hasn't been done so already.
    if (fileLoadFlag == false)
    {

      dof_handler.distribute_dofs (fe_collection);
      present_solution.reinit (dof_handler.n_dofs());
    }

    setup_system_constraints();

    unsigned int number_dofs = dof_handler.n_dofs();

    newton_update.reinit(number_dofs);
    system_rhs.reinit (number_dofs);
    drhs_dlambda.reinit (number_dofs);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

    system_matrix_petsc.reinit (sparsity_pattern);

    setup_bloch();

    evaluation_point = present_solution;

  }

  void ElasticProblem::setup_system_constraints ()
  {

    constraints.clear ();

    make_periodicity_constraints();
    make_ave_x1_constraints();
    make_symmetry_constraints();

    constraints.close ();

    // now do hanging nodes. Because some of the constraints might refer to the same dof
    // for both the symmetry constraint and the hanging node constraint, we will make them
    // separate, then merge them, giving precedence to the hanging node constraints;
    ConstraintMatrix hanging_node_constraints;
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close();

    constraints.merge(hanging_node_constraints, ConstraintMatrix::MergeConflictBehavior::right_object_wins);
  }

  void ElasticProblem::make_periodicity_constraints()
  {
    const FiniteElement<DIM, DIM> &fe = dof_handler.begin_active()->get_fe();
    (void)fe;


    const unsigned int   number_dofs = dof_handler.n_dofs();

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> is_x1_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    std::vector<bool> x1_x2_components = {true, true};
    ComponentMask x1_x2_mask(x1_x2_components);
    std::vector<bool> boundary_1_dof (number_dofs, false);
    std::vector<bool> boundary_2_dof (number_dofs, false);

    std::set< types::boundary_id > boundary_id_1;
    boundary_id_1.insert(1);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);

    DoFTools::extract_boundary_dofs(dof_handler,
                                    x1_x2_mask,
                                    boundary_1_dof,
                                    boundary_id_1);

    DoFTools::extract_boundary_dofs (dof_handler,
                                     x1_x2_mask,
                                     boundary_2_dof,
                                     boundary_id_2);

    // get the coords of the dofs
    std::vector<Point<DIM>> support_points(number_dofs);
    map_dofs_to_support_points(dof_handler, support_points);

    for(unsigned int i = 0; i < number_dofs; i++)
    {
      if(boundary_1_dof[i])
      {
        bool foundFlag = 0;
        for (unsigned int j = 0 ; j < number_dofs; j++)
        {
          if (boundary_2_dof[j] && (is_x1_comp[i] == is_x1_comp[j]) &&
              fabs(support_points[i](1) - support_points[j](1))< 1e-4 && i != j)
            // is a match
          {
            constraints.add_line(i);
            constraints.add_entry(i, j, 1.0);
            foundFlag = true;
          }
        }
        if(foundFlag == false)
        {
          std::cout << "could not find match for boundary 0 dof " << i << " Exiting." << std::endl;
          exit(-1);
        }
      }
    }
  }

  void ElasticProblem::make_ave_x1_constraints()
  {
    unsigned int number_dofs = dof_handler.n_dofs();

    // now do constraints that the average x1 displacement is zero
    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> boundary_dof_x1 (number_dofs, false);
    std::vector<bool> boundary_12_dof_x1 (number_dofs, false);

    std::set< types::boundary_id > boundary_id_12;
    boundary_id_12.insert(1);
    boundary_id_12.insert(2);

    DoFTools::extract_boundary_dofs(dof_handler,
                                    x1_mask,
                                    boundary_dof_x1);

    DoFTools::extract_boundary_dofs (dof_handler,
                                     x1_mask,
                                     boundary_12_dof_x1,
                                     boundary_id_12);



    unsigned int first_boundary_dof = 0;
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if ((boundary_dof_x1[i] == true) && (boundary_12_dof_x1[i] == false))
      {
        first_boundary_dof = i;
        break;
      }
    }
    constraints.add_line (first_boundary_dof);
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if (i == first_boundary_dof)
        continue;

      if(boundary_dof_x1[i] == true)
        constraints.add_entry (first_boundary_dof, i, -1);

    }
  }

  void ElasticProblem::make_symmetry_constraints()
  {
    unsigned int number_dofs = dof_handler.n_dofs();

    // get the vector that tells us it the dof is an x2 component
    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    std::vector<bool> is_x2_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);

    // get the vector that tells us if it is a boundary 3 x2 dof
    // (these are already constrained to zero)

    std::vector<bool> boundary_3_dof_x2 (number_dofs, false);

    std::set< types::boundary_id > boundary_id_3;
    boundary_id_3.insert(3);

    DoFTools::extract_boundary_dofs(dof_handler,
                                      x2_mask,
                                      boundary_3_dof_x2,
                                      boundary_id_3);

    std::vector<bool> x1_x2_components = {true, true};
    ComponentMask x1_x2_mask(x1_x2_components);
    std::vector<bool> boundary_1_dof (number_dofs, false);
    std::vector<bool> boundary_2_dof (number_dofs, false);

    std::set< types::boundary_id > boundary_id_1;
    boundary_id_1.insert(1);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);

    DoFTools::extract_boundary_dofs(dof_handler,
                                   x1_x2_mask,
                                   boundary_1_dof,
                                   boundary_id_1);

    DoFTools::extract_boundary_dofs (dof_handler,
                                    x1_x2_mask,
                                    boundary_2_dof,
                                    boundary_id_2);

    // get the coords of the dofs
    std::vector<Point<DIM>> support_points(number_dofs);
    map_dofs_to_support_points(dof_handler, support_points);

    matched_dofs.resize(number_dofs, -1);
    for(unsigned int i = 0; i < number_dofs; i++)
    {
     if (boundary_1_dof[i] || boundary_2_dof[i])
       continue;
     if (is_x2_comp[i])
     {
       if ( boundary_3_dof_x2[i] || matched_dofs[i] != -1 || (support_points[i](0) == 0.0))
       {
         // this dof has already been constrained or don't need to be
         continue;
       }
       else
       {
         double x1_coord = support_points[i](0);
         double x2_coord = support_points[i](1);
         for (unsigned int j = 0; j < number_dofs; j++)
         {
           if (is_x2_comp[j] && (fabs(x1_coord + support_points[j](0)) < 1e-12 ) &&
                                 (fabs(x2_coord - support_points[j](1)) < 1e-12) && i != j)
           {
             constraints.add_line (i);
             constraints.add_entry (i, j, 1);
             matched_dofs[i] = j;
             matched_dofs[j] = i;
             break;
           }
         }
       }
     }
     else
     {
       if ( matched_dofs[i] != -1 || (support_points[i](0) == 0.0))
       {
         // this dof has already been constrained or don't need to be
         continue;
       }
       else
       {
         double x1_coord = support_points[i](0);
         double x2_coord = support_points[i](1);
         for (unsigned int j = 0; j < number_dofs; j++)
         {
           if ((!is_x2_comp[j]) && (fabs(x1_coord + support_points[j](0)) < 1e-12 ) &&
                                 (fabs(x2_coord - support_points[j](1)) < 1e-12) && i != j)
           {
             constraints.add_line (i);
             constraints.add_entry (i, j, -1);
             matched_dofs[i] = j;
             matched_dofs[j] = i;
             break;
           }
         }
       }
     }
    }

  }

  void ElasticProblem::setup_bloch()
  {

    const unsigned int number_dofs = dof_handler.n_dofs();

    // first, setup the hanging nodes for the bloch constraints.
    // Have to manually do them for the bloch indicies
    bloch_hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler, bloch_hanging_node_constraints);

    ConstraintMatrix hanging_nodes_shifted;
    hanging_nodes_shifted.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, hanging_nodes_shifted);
    hanging_nodes_shifted.shift(number_dofs);
    hanging_nodes_shifted.close();

    bloch_hanging_node_constraints.merge(hanging_nodes_shifted);

    // now to do the sparsity stuff, have to make the bloch constraint matrix. Just use
    // some random wave number, say 0.1 (will give entries for all the ones that would
    // generally be there.
    update_bloch_wave_constraints(0.1);

    // manually do the sparsity stuff for the double sized matrix
    // this follows the dealii make_sparsity_pattern routine closely.
    DynamicSparsityPattern dsp_bloch(2*number_dofs, 2*number_dofs);

    std::vector<types::global_dof_index> local_dof_indices;
    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);
      constraints_bloch.add_entries_local_to_global (local_dof_indices,
                                                            dsp_bloch,
                                                           /*keep_constrained_dofs*/ true);
      // increment them by the number of dofs
      for(unsigned int i = 0; i < dofs_per_cell; i++)
        local_dof_indices[i] += number_dofs;


      constraints_bloch.add_entries_local_to_global (local_dof_indices,
                                                            dsp_bloch,
                                                           /*keep_constrained_dofs*/ true);
    }

    sparsity_pattern_bloch.copy_from(dsp_bloch);

    bloch_matrix.reinit(sparsity_pattern_bloch);


//    stiffness_matrix.reinit (2*number_dofs, 2*number_dofs,
//                             dof_handler.max_couplings_between_dofs());

    stiffness_matrix.reinit(sparsity_pattern_bloch);

    // in parallel.
    IndexSet local_dofs_index_set = dof_handler.locally_owned_dofs ();
    IndexSet eigenfunction_index_set(2*number_dofs);
    for(unsigned int i = 0 ; i < number_dofs; i ++)
    {
      if(local_dofs_index_set.is_element(i))
      {
        eigenfunction_index_set.add_index(i);
        eigenfunction_index_set.add_index(i + number_dofs);
      }
    }

    eigenfunctions_.resize (5);
    for (unsigned int i=0; i<eigenfunctions_.size (); ++i)
      eigenfunctions_[i].reinit (eigenfunction_index_set, MPI_COMM_WORLD);

    eigenvalues_.resize (eigenfunctions_.size ());
  }

  void ElasticProblem::setup_bloch_matched_dofs()
  {
    bloch_boundry_1_dof_index.clear();
    bloch_boundry_2_dof_index.clear();


    // get the boundary 1 and boundary 2 dofs
    const unsigned int   number_dofs = dof_handler.n_dofs();

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> boundary_1_dof (number_dofs, false);
    std::vector<bool> boundary_2_dof (number_dofs, false);

    std::set< types::boundary_id > boundary_id_1;
    boundary_id_1.insert(1);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);

    std::vector<bool> x1_x2_components = {true, true};
    ComponentMask x1_x2_mask(x1_x2_components);

    DoFTools::extract_boundary_dofs(dof_handler,
                                   x1_x2_mask,
                                   boundary_1_dof,
                                   boundary_id_1);

    DoFTools::extract_boundary_dofs (dof_handler,
                                    x1_x2_mask,
                                    boundary_2_dof,
                                    boundary_id_2);

    std::vector<bool> is_x1_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    // get the coords of the dofs
    std::vector<Point<DIM>> support_points(number_dofs);
    map_dofs_to_support_points(dof_handler, support_points);

    for(unsigned int i = 0; i < number_dofs; i ++)
    {
     if (boundary_1_dof[i])
     {
       double x2_coord = support_points[i](1);
       for (unsigned int j = 0; j < number_dofs; j++)
       {
         if (boundary_2_dof[j] &&
             (fabs(x2_coord - support_points[j](1)) < 1e-12 ) &&
             (is_x1_comp[i] == is_x1_comp[j]) )
         {

           bloch_boundry_1_dof_index.push_back(i);
           bloch_boundry_2_dof_index.push_back(j);

           break;
         }
       }
     }
    }

    bloch_matched_flag = true;

  }

  void ElasticProblem::update_bloch_wave_constraints(double wave_ratio)
  {

    if(!bloch_matched_flag)
    {
      setup_bloch_matched_dofs();
    }
    constraints_bloch.clear();

    double k_x1 = 2.0*(4.0*atan(1.0))*wave_ratio;

    // get the boundary 1 and boundary 2 dofs
    const unsigned int   number_dofs = dof_handler.n_dofs();

    for (unsigned int k = 0; k < bloch_boundry_1_dof_index.size(); k++)
    {
      unsigned int i = bloch_boundry_1_dof_index[k];
      unsigned int j = bloch_boundry_2_dof_index[k];

      constraints_bloch.add_line (i);
      constraints_bloch.add_entry (i, j, cos(k_x1)); //-cos(k_x1));
      constraints_bloch.add_entry (i, (j + number_dofs) , sin(k_x1)); //-cos(k_x1));

      constraints_bloch.add_line (i + number_dofs);
      constraints_bloch.add_entry ((i + number_dofs), j, -sin(k_x1));
      constraints_bloch.add_entry ((i + number_dofs), (j + number_dofs), cos(k_x1));
    }

    constraints_bloch.close();

    // Because some of the constraints might refer to the same dof
    // for both the bloch constraint and the hanging node constraint, we will make them
    // separate, then merge them, giving precedence to the hanging node constraints.

    constraints_bloch.merge(bloch_hanging_node_constraints, ConstraintMatrix::MergeConflictBehavior::right_object_wins);

  }

  void ElasticProblem::set_boundary_values()
  {
    // this sets the boundary values of the solution vector so that the Newton step
    // can use homogeneous direchlet conditions on the set boundaries. It also makes sure
    // that the periodic faces' DoFs start with the same values (zero).

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::map<types::global_dof_index,double> boundary_values;

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 1,
                                                 ZeroFunction<DIM>(DIM),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 2,
                                                 ZeroFunction<DIM>(DIM),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 3,
                                                 ZeroFunction<DIM>(DIM),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 4,
                                                 ZeroFunction<DIM>(DIM),
                                                 boundary_values,
                                                 x1_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        present_solution,
                                        system_rhs);
    for (std::map<types::global_dof_index, double>::const_iterator
         p = boundary_values.begin();
         p != boundary_values.end(); ++p)
      present_solution(p->first) = p->second;

  }

  void ElasticProblem::update_F0(const double lambda)
  {
    // Update the uniform diagonal F0 tensor.

    Assert ((lambda >= 0 && lambda <= 1.0), ExcInternalError());

    double lambda1 = 1 - lambda;
    double lambda2;

    // solve quadratic for lambda2
    double a = (1.0 + (2.0*NU_VALUE/(1.0 - NU_VALUE))*lambda1*lambda1);
    double b = (-2.0*NU_VALUE/(1.0 - NU_VALUE))*lambda1;
    double c = -1;


    // make sure we get the root that is greater than 1.0
    lambda2 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);

    F0[0][0] = lambda1;
    F0[1][1] = lambda2;
    F0[0][1] = 0.0;
    F0[1][0] = 0.0;
  }

  void ElasticProblem::add_small_pertubations(double amplitude, bool firstTime)
  {
    // This is just to make the solution vector non-zero so it doesn't start at the right answer.

    if(firstTime)
      std::srand(5); //some seed

    for(unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      if (matched_dofs[i] != -1)
      {
        double random_val =  0.5*(2.0*(std::rand()/double(RAND_MAX)) - 1.0)*amplitude;
        present_solution[i] += random_val;
        present_solution[matched_dofs[i]] += random_val;
      }
      else
        present_solution[i] += (2.0*(std::rand()/double(RAND_MAX)) - 1.0)*amplitude;
    }

    evaluation_point = present_solution;

  }

  void ElasticProblem::assemble_system_matrix()
  {
    // Assembling the system matrix. I chose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_matrix = 0.0;

//    QGauss<1> quad_x(n_qpoints_x);
//    QGauss<1> quad_y(n_qpoints_y);
//
//    QAnisotropic<DIM> quadrature_formula(quad_x, quad_y);

//    QGauss<DIM>  quadrature_formula(1);

    hp::FEValues<DIM> hp_fe_values (fe_collection, quadrature_collection,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

//    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix;


    std::vector<types::global_dof_index> local_dof_indices;


    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_matrix = 0.0;

      hp_fe_values.reinit (cell);
      const FEValues<DIM> &fe_values = hp_fe_values.get_present_fe_values();
      unsigned int n_q_points = fe_values.n_quadrature_points;

      std::vector<double>     nu_values(n_q_points);
      std::vector<double>     mu_values(n_q_points);

      std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,DIM>>(DIM));

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<4,DIM> d2W_dFdF = nh.get_incremental_moduli_tensor(nu_values[q_point],
                                           mu_values[q_point], F_inv, II_F);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = cell->get_fe().system_to_component_index(n).first;

          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {
            const unsigned int component_m = cell->get_fe().system_to_component_index(m).first;

            for (unsigned int j=0; j<DIM; ++j)
              for (unsigned int l=0; l<DIM; ++l)
              {
                cell_matrix(n,m) += d2W_dFdF[component_n][j][component_m][l]*
                        fe_values.shape_grad(n, q_point)[j]*fe_values.shape_grad(m, q_point)[l]*fe_values.JxW(q_point);

              }
          }
        }
      }
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        for (unsigned int m=0; m<dofs_per_cell; ++m)
        {
          system_matrix.add (local_dof_indices[n],
                             local_dof_indices[m],
                             cell_matrix(n,m));
        }
    }
  }

  void ElasticProblem::assemble_system_rhs()
  {
    // Assembling the system rhs. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_rhs = 0.0;

    hp::FEValues<DIM> hp_fe_values (fe_collection, quadrature_collection,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);


    Vector<double>       cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;


    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      cell_rhs.reinit(dofs_per_cell);
      cell_rhs = 0.0;

      hp_fe_values.reinit (cell);
      const FEValues<DIM> &fe_values = hp_fe_values.get_present_fe_values();
      unsigned int n_q_points = fe_values.n_quadrature_points;
      std::vector<Tensor<1, DIM> > rhs_values (n_q_points);

      std::vector<double>     nu_values (n_q_points);
      std::vector<double>     mu_values (n_q_points);

      std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,DIM>>(DIM));

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);
      right_hand_side (fe_values.get_quadrature_points(), rhs_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<2,DIM> dW_dF = nh.get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point],
                                                        F, F_inv, II_F);
        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = cell->get_fe().system_to_component_index(n).first;

          for(unsigned int j = 0; j<DIM; ++j)
          {
            cell_rhs(n) -= dW_dF[component_n][j]*fe_values.shape_grad(n, q_point)[j]*fe_values.JxW(q_point);
          }

          cell_rhs(n) += fe_values.shape_value(n, q_point)*rhs_values[q_point][component_n]*fe_values.JxW(q_point);
        }
      }
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        system_rhs(local_dof_indices[n]) += cell_rhs(n);

    }

    constraints.condense (system_rhs);

  }

  void ElasticProblem::assemble_system_energy_and_congugate_lambda()
  {
    double lambda_eval = present_lambda;

    system_energy = 0.0;
    congugate_lambda = 0.0;

    hp::FEValues<DIM> hp_fe_values (fe_collection, quadrature_collection,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    std::vector<types::global_dof_index> local_dof_indices;

    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      hp_fe_values.reinit (cell);
      const FEValues<DIM> &fe_values = hp_fe_values.get_present_fe_values();
      unsigned int n_q_points = fe_values.n_quadrature_points;

      std::vector<Tensor<1, DIM> > rhs_values (n_q_points);

      std::vector<double>     nu_values (n_q_points);
      std::vector<double>     mu_values (n_q_points);

      std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points,
                                                  std::vector<Tensor<1,DIM>>(DIM));

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<2,DIM> dW_dF = nh.get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point],
                                                        F, F_inv, II_F);

        // get the dF_dlambda
        double dlambda1_dlambda, dlambda2_dlambda, a_nu;

        dlambda1_dlambda = -1.0;

        a_nu = (2.0*NU_VALUE/(1.0 - NU_VALUE));
        double lambda1 = 1.0 - lambda_eval;
        double term1 = 2.0*a_nu*(1.0 + a_nu*lambda1*lambda1);
        double term2 = -4.0*a_nu*a_nu*lambda1*lambda1;
        double term3 = (2.0*a_nu*a_nu*lambda1 + 8.0*a_nu*lambda1)*(1.0 + a_nu * lambda1*lambda1)/
                              sqrt(a_nu*a_nu*lambda1*lambda1 + 4.0 *a_nu*lambda1*lambda1 + 4.0);
        double term4 = -sqrt(a_nu*a_nu*lambda1*lambda1 + 4.0 *a_nu*lambda1*lambda1 + 4.0)*4.0*a_nu*lambda1;

        dlambda2_dlambda = -(term1 + term2 + term3 + term4)/(4.0*(1 + a_nu*lambda1*lambda1)*(1 + a_nu*lambda1*lambda1));

        congugate_lambda += (dlambda1_dlambda*dW_dF[0][0] + dlambda2_dlambda*dW_dF[1][1])
                             *fe_values.JxW(q_point);

        double W = nh.get_energy(nu_values[q_point], mu_values[q_point], F, II_F);
        system_energy += W*fe_values.JxW(q_point);
      }
    }

  }

  void ElasticProblem::assemble_drhs_dlambda()
  {
    double lambda_eval = present_lambda;
    update_F0(lambda_eval);

    drhs_dlambda = 0.0;

    hp::FEValues<DIM> hp_fe_values (fe_collection, quadrature_collection,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);


    Vector<double>       cell_drhs_dlambda;

    std::vector<types::global_dof_index> local_dof_indices;




    // construct the dF_dlambda
    double dlambda1_dlambda, dlambda2_dlambda, a_nu;
    dlambda1_dlambda = -1.0;
    a_nu = (2.0*NU_VALUE/(1.0 - NU_VALUE));
    double lambda1 = 1.0 - lambda_eval;
    double term1 = 2.0*a_nu*(1.0 + a_nu*lambda1*lambda1);
    double term2 = -4.0*a_nu*a_nu*lambda1*lambda1;
    double term3 = (2.0*a_nu*a_nu*lambda1 + 8.0*a_nu*lambda1)*(1.0 + a_nu * lambda1*lambda1)/
                          sqrt(a_nu*a_nu*lambda1*lambda1 + 4.0 *a_nu*lambda1*lambda1 + 4.0);
    double term4 = -sqrt(a_nu*a_nu*lambda1*lambda1 + 4.0 *a_nu*lambda1*lambda1 + 4.0)*4.0*a_nu*lambda1;
    dlambda2_dlambda = -(term1 + term2 + term3 + term4)/(4.0*(1 + a_nu*lambda1*lambda1)*(1 + a_nu*lambda1*lambda1));

    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      cell_drhs_dlambda.reinit(dofs_per_cell);
      cell_drhs_dlambda = 0.0;

      hp_fe_values.reinit (cell);
      const FEValues<DIM> &fe_values = hp_fe_values.get_present_fe_values();
      unsigned int n_q_points = fe_values.n_quadrature_points;

      std::vector<double>     nu_values (n_q_points);
      std::vector<double>     mu_values (n_q_points);

       std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,DIM>>(DIM));

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<4,DIM> d2W_dFdF = nh.get_incremental_moduli_tensor(nu_values[q_point],
                                           mu_values[q_point], F_inv, II_F);

        Tensor<2,DIM> dW_dF_dlambda;
        for (unsigned int i = 0 ; i < DIM; i ++)
          for (unsigned int j = 0 ; j < DIM; j ++)
              dW_dF_dlambda[i][j] = d2W_dFdF[i][j][0][0]*dlambda1_dlambda  +  d2W_dFdF[i][j][1][1]*dlambda2_dlambda;



        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = cell->get_fe().system_to_component_index(n).first;

          for(unsigned int j = 0; j<DIM; ++j)
          {
            cell_drhs_dlambda(n) += dW_dF_dlambda[component_n][j]*fe_values.shape_grad(n, q_point)[j]*fe_values.JxW(q_point);
          }
        }
      }
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        drhs_dlambda(local_dof_indices[n]) += cell_drhs_dlambda(n);
    }

    constraints.condense (drhs_dlambda);
  }

  void ElasticProblem::assemble_bloch_matrix()
  {
    bloch_matrix = 0.0;

    // just copy the system matrix into the upper left and lower right blocks of
    // the bloch_matrix
    unsigned int number_dofs = dof_handler.n_dofs();
    for (unsigned int row = 0; row < system_matrix.m(); ++row)
    {
      const typename SparseMatrix<double>::const_iterator end_row = system_matrix.end(row);
      for (typename SparseMatrix<double>::const_iterator entry = system_matrix.begin(row);
                             entry != end_row; ++entry)
       {
         bloch_matrix.set(row, entry->column(),entry->value());
         bloch_matrix.set(row + number_dofs, entry->column() + number_dofs,entry->value());
       }
    }

  }

  void ElasticProblem::compute_first_bif_amplitude()
  {

    u1u1 = 0.0;
    excee = 0.0;

    hp::FEValues<DIM> hp_fe_values (fe_collection, quadrature_collection,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);


    typename hp::DoFHandler<DIM>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      hp_fe_values.reinit (cell, last_q_index);
      const FEValues<DIM> &fe_values = hp_fe_values.get_present_fe_values();
      unsigned int n_q_points = fe_values.n_quadrature_points;

      std::vector<Vector<double>> solution_vals(n_q_points,
                                                Vector<double>(DIM));

      std::vector<Vector<double>> u1(n_q_points, Vector<double>(DIM));

      fe_values.get_function_values(present_solution, solution_vals);

      u1_value_list(fe_values.get_quadrature_points(), u1);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        excee += solution_vals[q_point]*(u1[q_point])*fe_values.JxW(q_point);
        u1u1 += (u1[q_point])*(u1[q_point])*fe_values.JxW(q_point);
      }

    }


    u1u1 = u1u1/domain_dimensions[0];
    excee = excee/(domain_dimensions[0]);
  }

  void ElasticProblem::u1_value_list(const std::vector< Point< DIM > > &  points, std::vector<Vector<double>> & value_list)
  {
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
     {
       Vector<double> tmp(2);
       Vector<std::complex<double>> tmp_complex(2);
       double x1 = points[point_n](0);
       double x2 = points[point_n](1);

       unsigned int t = 0;
       if(x2 > L1)
         t = 4;

       for(unsigned int i = 0; i < 4; i++)
       {
         std::complex<double> Ae = A[i + t]*std::exp(alphas[i]*x2);
         std::complex<double> AeB = Ae*B[i];

         tmp_complex[0] += Ae;
         tmp_complex[1] += AeB;

       }

       tmp[0] = -sin(w_c*x1)*(tmp_complex[0]).real();
       tmp[1] = cos(w_c*x1)*(tmp_complex[1]).real();

       value_list[point_n] = tmp;
     }
  }

  void ElasticProblem::apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs)
  {
    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      if ((*homogenous_dirichlet_dofs)[i] == true)
        (*rhs)[i] = 0.0;
    }
  }

  void ElasticProblem::apply_boundaries_and_constraints_system_matrix()
  {
    constraints.condense (system_matrix);

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side3_components = {false, true};
    ComponentMask side3_mask(side3_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              3,
                                              ZeroFunction<DIM, double>(DIM),
                                              boundary_values,
                                              side3_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);

  }

  void ElasticProblem::apply_boundaries_and_constraints_bloch_matrix()
  {
    constraints_bloch.condense(bloch_matrix);

    const unsigned int number_dofs = dof_handler.n_dofs();

    std::vector<bool> side2_components = {false, true};
    ComponentMask x2_mask(side2_components);
    std::vector<bool> boundary_3_dof_x2 (number_dofs, false);

    std::set< types::boundary_id > boundary_id_3;
    boundary_id_3.insert(3);

    DoFTools::extract_boundary_dofs(dof_handler,
                                      x2_mask,
                                      boundary_3_dof_x2,
                                      boundary_id_3);

    unsigned int m = bloch_matrix.m();
    // set values on the diagonal to the first diagonal element,
    // or 1 if it is nonexistent
    // This all follows the dealii built in apply_boundaries closely
    double first_nonzero_diagonal_entry = 1.0;
    for (unsigned int i=0; i<m; ++i)
    {
      if (bloch_matrix.diag_element(i) != 0.0)
      {
        first_nonzero_diagonal_entry = fabs(bloch_matrix.diag_element(i));
        break;
      }
    }
    // now march through matrix, zeroing out rows and columns.
    // If there is a current value on the diagonal of the constrained
    // boundary dof, don't touch it. If there is not one, then we can
    // just set it equal to the first nonzero entry we just found
    for (unsigned int row = 0; row < m; ++row)
    {

      const typename SparseMatrix<double>::iterator end_row = bloch_matrix.end(row);
      for (typename SparseMatrix<double>::iterator entry = bloch_matrix.begin(row);
                    entry != end_row; ++entry)
      {
        if((boundary_3_dof_x2[row%number_dofs] || boundary_3_dof_x2[entry->column()%number_dofs])
            && (row != entry->column()))
        {
          entry->value() = 0.0;
        }
        else if(boundary_3_dof_x2[row%number_dofs] && (row == entry->column()) &&
                (bloch_matrix.diag_element(row%number_dofs) == 0.0))
        {
          entry->value() = first_nonzero_diagonal_entry;
        }

      }

    }

   }

  void ElasticProblem::newton_iterate()
  {
    /* This function executes the newton iterations until it converges to a solution
     * or it exceeds the maximum number of iterations.
     */

    double current_residual;
    unsigned int iteration = 0;

    // get the dofs that we will apply dirichlet condition to
    std::vector<bool> boundary_3_dof_x2 (dof_handler.n_dofs(), false);

    std::set< types::boundary_id > boundary_id_3;
    boundary_id_3.insert(3);


    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    DoFTools::extract_boundary_dofs(dof_handler,
                                       x2_mask,
                                       boundary_3_dof_x2,
                                       boundary_id_3);

    // Starts by getting the residual in the current configuration.
    evaluation_point = present_solution;

    assemble_system_rhs();
    apply_boundaries_to_rhs(&system_rhs, &boundary_3_dof_x2);

    current_residual = system_rhs.l2_norm();

    // Loops until coverge or go over max iterations
    while((current_residual > tol) &&
             (iteration < maxIter))
    {
      // Assemble the stiffness matrix
      assemble_system_matrix();
      apply_boundaries_and_constraints_system_matrix();

      // solve for the newton step
      solve();

      // Find the step length and add it to the current solution.
      // This function also calls assemble_system_rhs() so we don't need to
      // do another rhs call.
      line_search_and_add_step_length(current_residual, &boundary_3_dof_x2);

      evaluation_point = present_solution;
      current_residual = system_rhs.l2_norm();
      iteration ++;
    }

    // output iterations for convergance.
    std::cout << "    Converging Iterations : " << iteration << "\n";

  }

  void ElasticProblem::path_follow_PACA_iterate(Vector<double> *solVectorDir,
                                                          double lambdaDir, double ds)
  {
    // The setup stuff, may need to take  in more inputs. But for now we will do with this

    Vector<double> previousSolution = present_solution;
    double previousLambda = present_lambda;
    double current_residual = 0.0;
    unsigned int iteration = 0;

    // Scale the input tangent to unity just incase it wasnt already
    double scalingVal = solVectorDir->norm_sqr() + lambdaDir*lambdaDir;
    scalingVal = 1.0/sqrt(scalingVal);
    lambdaDir *= scalingVal;
    (*solVectorDir) *= scalingVal;

    // get the dofs that we will apply dirichlet condition to
    std::vector<bool> boundary_3_dof_x2 (dof_handler.n_dofs(), false);

    std::set< types::boundary_id > boundary_id_3;
    boundary_id_3.insert(3);


    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    DoFTools::extract_boundary_dofs(dof_handler,
                                       x2_mask,
                                       boundary_3_dof_x2,
                                       boundary_id_3);

    // Starts by getting the residual in the initial guess.
    present_lambda += lambdaDir*ds;
    present_solution.add(ds, *solVectorDir);
    evaluation_point = present_solution;

    update_F0(present_lambda);

    assemble_system_rhs();
    apply_boundaries_to_rhs(&system_rhs, &boundary_3_dof_x2);

    lambda_diff = present_lambda - previousLambda;
    solution_diff = evaluation_point;
    solution_diff -= previousSolution;

    current_residual = system_rhs.norm_sqr() +
                            0.25*( solution_diff.norm_sqr() + lambda_diff*lambda_diff - ds*ds)*
                            ( solution_diff.norm_sqr() + lambda_diff*lambda_diff - ds*ds);
    current_residual = sqrt(current_residual);

    while ((current_residual > tol) &&
        (iteration < maxIter))
    {

      evaluation_point = present_solution;

      update_F0(present_lambda);

      // Assemble rhs
      assemble_system_rhs();
      apply_boundaries_to_rhs(&system_rhs, &boundary_3_dof_x2);

      // Assemble the stiffness matrix
      assemble_system_matrix();
      apply_boundaries_and_constraints_system_matrix();

      // Assemble the drhs_dlambda
      assemble_drhs_dlambda();
      apply_boundaries_to_rhs(&drhs_dlambda, &boundary_3_dof_x2);

      // Get the solution diff (bottom boarder of system mat)
      solution_diff = present_solution;
      solution_diff -= previousSolution;
      apply_boundaries_to_rhs(&solution_diff, &boundary_3_dof_x2);

      // get the bottom corner of system mat
      lambda_diff = present_lambda - previousLambda;

      // get bottom element of the rhs
      rhs_bottom = -0.5*( solution_diff.norm_sqr() + lambda_diff*lambda_diff - ds*ds);

      // Solve it with the boarder matrix solver
      solve_boarder_matrix_system();

      // add the update using the line search
      line_search_and_add_step_length_PACA(current_residual, &boundary_3_dof_x2, &previousSolution, previousLambda, ds);

      solution_diff = present_solution;
      solution_diff -= previousSolution;

      current_residual = system_rhs.norm_sqr() +
                              0.25*( solution_diff.norm_sqr() + lambda_diff*lambda_diff - ds*ds)*
                              ( solution_diff.norm_sqr() + lambda_diff*lambda_diff - ds*ds);
      current_residual = sqrt(current_residual);

      iteration ++;

    }
    // output iterations for convergance.
    std::cout << "\n    Converging Iterations : " << iteration << "         Residual : " << current_residual << std::endl;
  }

  void ElasticProblem::line_search_and_add_step_length_PACA(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs,
                                                                 Vector<double> *previousSolution, double previousLambda, double ds)
  {
   /* this function makes sure that the step sizes we are taking with
    * the newton iteration are making the residual get smaller.
    * Something very similar is used in the dealii step 57?. for now
    * it doesn't do much but we might need to give it more capabilites in the future.
    */

    double current_residual = 0.0;
    double lambda_eval = 0.0;
    double lambdaDiff = 0.0;
    Vector<double> solutionDiff;

    for(double alpha = 1.0; alpha > 1e-5; alpha *=0.5)
    {
      evaluation_point = present_solution;
      evaluation_point.add(alpha, newton_update);

      lambda_eval = present_lambda + alpha*lambda_update;

      // want to just skip if the lambda eval is out of range...
      if (lambda_eval > 1.0 || lambda_eval < 0.0)
        continue;

      update_F0(lambda_eval);
      assemble_system_rhs();
      apply_boundaries_to_rhs(&system_rhs, homogenous_dirichlet_dofs);


      solutionDiff = present_solution;
      solutionDiff -= *previousSolution;
      lambdaDiff = lambda_eval - previousLambda;

      current_residual = system_rhs.norm_sqr() +
                              0.25*( solutionDiff.norm_sqr() + lambdaDiff*lambdaDiff - ds*ds)*
                              ( solutionDiff.norm_sqr() + lambdaDiff*lambdaDiff - ds*ds);
      current_residual = sqrt(current_residual);
      if(current_residual < last_residual)
        break;

    }
    present_solution = evaluation_point;
    present_lambda = lambda_eval;
  }

  void ElasticProblem::line_search_and_add_step_length(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs)
  {
   /* this function makes sure that the step sizes we are taking with
    * the newton iteration are making the residual get smaller.
    * Something very similar is used in the dealii step 57?. for now
    * it doesn't do much but we might need to give it more capabilites in the future.
    */

    double current_residual;
    for(double alpha = 1.0; alpha > 1e-5; alpha *=0.5)
    {
      evaluation_point = present_solution;
      evaluation_point.add(alpha, newton_update);

      assemble_system_rhs();
      apply_boundaries_to_rhs(&system_rhs, homogenous_dirichlet_dofs);

      current_residual = system_rhs.l2_norm();

      if(current_residual < last_residual)
        break;

    }
    present_solution = evaluation_point;
  }

  void ElasticProblem::solve ()
  {


    // direct solver for the system

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);
    constraints.distribute (newton_update);

  }

  void ElasticProblem::solve_boarder_matrix_system()
  {

    // steps to solving board system as described in
    // Govaerts's Numerical methods for bifurcation of dynamical equillibria 3.6.2

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);

    Vector<double> omega;
    Vector<double> v;
    Vector<double> excee;
    Vector<double> f1;
    double delta_star, delta, y1, g1, y2;

    // step 1
    A_direct.vmult(omega, solution_diff);

    // step 2
    delta_star = lambda_diff - omega*drhs_dlambda;

    // step 3
    A_direct.vmult(v, drhs_dlambda);

    // step 4
    delta = lambda_diff - solution_diff*v;

    // step 5
    y1 =(rhs_bottom - omega*system_rhs)/delta_star;

    // step 6
    f1 = system_rhs;
    f1.add(-y1, drhs_dlambda);

    // step 7
    g1 = rhs_bottom - lambda_diff*y1;

    // step 8
    A_direct.vmult(excee, f1);

    // step 9
    y2 = (g1 - solution_diff*excee)/delta;

    // step 10
    newton_update = excee;
    newton_update.add(-y2, v);

    // step 11
    lambda_update = y1 + y2;

    constraints.distribute (newton_update);
  }

  unsigned int ElasticProblem::get_system_eigenvalues(double lambda_eval, const int cycle)
  {
    // get the system's eigenvalues. I really don't need to have it take in lambda_eval
    // and could just have it use the present_lambda, but its fine. the cycle is the
    // the output file will have appended on its name. -1 for no output. Outputs the number
    // of negative eigenvalues in the range -10 < lambda < 0.3

    update_F0(lambda_eval);
    evaluation_point = present_solution;
    assemble_system_matrix();
    apply_boundaries_and_constraints_system_matrix();

//    // copy current sparse system matrix to a full matrix.
//    LAPACKFullMatrix<double> system_matrix_full;
//    system_matrix_full.copy_from(system_matrix);
//
//    Vector<double> eigenvalues;
//    FullMatrix<double> eigenvectors;
//
//    system_matrix_full.compute_eigenvalues_symmetric(-10, 0.3, 1e-5, eigenvalues, eigenvectors);

//    unsigned int numEigs = 20;
//
//    std::vector<Vector<double> >        eigenfunctions(numEigs);
//    for(unsigned int i = 0; i < numEigs; i++)
//      eigenfunctions[i].reinit(2*dof_handler.n_dofs());
//
//    std::vector<std::complex<double>>   eigenvalues(numEigs);
//    SolverControl solver_control (2*dof_handler.n_dofs(), 1e-9);
//    SparseDirectUMFPACK inverse;
//    inverse.initialize (system_matrix);
//    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;
//    ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors, ArpackSolver::WhichEigenvalues::largest_magnitude, true);   //, ArpackSolver::WhichEigenvalues::smallest_magnitude, false);
//    ArpackSolver eigensolver (solver_control, additional_data);
//    eigensolver.solve (system_matrix,
//                       IdentityMatrix(2*dof_handler.n_dofs()),
//                       inverse,
//                       eigenvalues,
//                       eigenfunctions,
//                       eigenvalues.size(),
//                       0);


    std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
    std::vector<double>                     eigenvalues;
    IndexSet eigenfunction_index_set(dof_handler.n_dofs());
    eigenfunction_index_set.add_range(0, dof_handler.n_dofs());
    eigenfunctions.resize (1);
    for (unsigned int i=0; i<eigenfunctions.size (); ++i)
      eigenfunctions[i].reinit (eigenfunction_index_set, MPI_COMM_WORLD);

    eigenvalues.resize (eigenfunctions.size());

    system_matrix_petsc = 0.0;
    for (unsigned int row = 0; row < system_matrix.m(); ++row)
    {
      const typename SparseMatrix<double>::const_iterator end_row = system_matrix.end(row);
      for (typename SparseMatrix<double>::const_iterator entry = system_matrix.begin(row);
                             entry != end_row; ++entry)
       {
         if(fabs(entry->value()) > 1e-10 )
           system_matrix_petsc.set(row, entry->column(),entry->value());
       }
    }

    system_matrix_petsc.compress(VectorOperation::insert);
    SolverControl solver_control (dof_handler.n_dofs(), 1e-6);
    SLEPcWrappers::SolverKrylovSchur eigensolver (solver_control);
//    SLEPcWrappers::SolverLAPACK eigensolver (solver_control);

    // Before we actually solve for the eigenfunctions and -values, we have to
    // also select which set of eigenvalues to solve for. Lets select those
    // eigenvalues and corresponding eigenfunctions with the smallest real
    // part (in fact, the problem we solve here is symmetric and so the
    // eigenvalues are purely real). After that, we can actually let SLEPc do
    // its work:
//    eigensolver.set_target_eigenvalue(-10.0);
    eigensolver.set_which_eigenpairs (EPS_SMALLEST_REAL);

    eigensolver.set_problem_type (EPS_HEP);

    eigensolver.solve (system_matrix_petsc,
                       eigenvalues, eigenfunctions,
                       eigenfunctions.size());


    unsigned int num_neg_eigs = 0;
    if (cycle != -1)
    {
      std::string filename(output_directory);
          filename += "/eigenvalues";

      // see if the directory exists
      struct stat st;
      if (stat(filename.c_str(), &st) == -1)
        mkdir(filename.c_str(), 0700);

      filename += "/eigenvalues-";
      filename += std::to_string(cycle);

      std::ofstream outputFile;
      outputFile.open(filename.c_str());

      //outputFile << "# eigenvalues of the system matrix" << std::endl;

      for (unsigned int i = 0 ; i < eigenvalues.size(); i ++)
      {
        double nextEigenVal = eigenvalues[i];

        outputFile << std::setprecision(15) << nextEigenVal << std::endl;

        if (nextEigenVal < 0.0)
        {
          num_neg_eigs ++;
        }

      }

      // outputFile << "\nIs positive definite : " << num_neg_eigs << std::endl;
      outputFile.close();
    }
    else
    {
      for (unsigned int i = 0; i < eigenvalues.size(); i++)
      {
        if (eigenvalues[i] < 0.0)
        {
          num_neg_eigs ++;
          break;
        }
      }
    }

    return num_neg_eigs;
  }

  double ElasticProblem::get_bloch_eigenvalues(const int cycle, const int step, double wave_ratio, unsigned int indx, bool first)
  {
    // The cycle is the number appended onto the end of the file. Step is the current step number,
    // and wave_ratio is the wave ratio value between 0 and 0.5

    if(first)
    {
      evaluation_point = present_solution;
      update_F0(present_lambda);
      assemble_system_matrix();
    }
    assemble_bloch_matrix();

    update_bloch_wave_constraints(wave_ratio);

    apply_boundaries_and_constraints_bloch_matrix();

//    stiffness_matrix = 0;
//    for (unsigned int row = 0; row < bloch_matrix.m(); ++row)
//    {
//      const typename SparseMatrix<double>::const_iterator end_row = bloch_matrix.end(row);
//      for (typename SparseMatrix<double>::const_iterator entry = bloch_matrix.begin(row);
//                             entry != end_row; ++entry)
//       {
//         if(fabs(entry->value()) > 1e-8)
//           stiffness_matrix.set(row, entry->column(),entry->value());
//       }
//    }
//
//    stiffness_matrix.compress(VectorOperation::insert);
//    SolverControl solver_control (2*dof_handler.n_dofs(), 1e-5);
//    SLEPcWrappers::SolverKrylovSchur eigensolver (solver_control);
////    SLEPcWrappers::SolverLAPACK eigensolver (solver_control);
//
//    // Before we actually solve for the eigenfunctions and -values, we have to
//    // also select which set of eigenvalues to solve for. Lets select those
//    // eigenvalues and corresponding eigenfunctions with the smallest real
//    // part (in fact, the problem we solve here is symmetric and so the
//    // eigenvalues are purely real). After that, we can actually let SLEPc do
//    // its work:
////    eigensolver.set_target_eigenvalue(-10.0);
//    eigensolver.set_which_eigenpairs (EPS_SMALLEST_REAL);
//
//    eigensolver.set_problem_type (EPS_HEP);
//
//    eigensolver.solve (stiffness_matrix,
//                       eigenvalues_, eigenfunctions_,
//                       eigenfunctions_.size());

    // copy current sparse system matrix to a full matrix.

//    LAPACKFullMatrix<double> bloch_matrix_full;
//    bloch_matrix_full.copy_from(bloch_matrix);
//
//
//    Vector<double> eigenvalues;
//    FullMatrix<double> eigenvectors;
//
//    bloch_matrix_full.compute_eigenvalues_symmetric(-10, 0.1, 1e-5, eigenvalues, eigenvectors);


    unsigned int numEigs = 20;

    std::vector<Vector<double> >        eigenfunctions(numEigs+1);
    for(unsigned int i = 0; i < numEigs; i++)
      eigenfunctions[i].reinit(2*dof_handler.n_dofs());

    std::vector<std::complex<double>>   eigenvalues(numEigs);
    SolverControl solver_control (2*dof_handler.n_dofs(), 1e-9);
    SparseDirectUMFPACK inverse;
    inverse.initialize (bloch_matrix);
    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;
    ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors, ArpackSolver::WhichEigenvalues::largest_magnitude, true);   //, ArpackSolver::WhichEigenvalues::smallest_magnitude, false);
    ArpackSolver eigensolver (solver_control, additional_data);
    eigensolver.solve (bloch_matrix,
                       IdentityMatrix(2*dof_handler.n_dofs()),
                       inverse,
                       eigenvalues,
                       eigenfunctions,
                       eigenvalues.size(),
                       0);

    std::string filename(output_directory);
    filename += "/bloch_eigenvalues";

    // see if the directory exists
    struct stat st;
    if (stat(filename.c_str(), &st) == -1)
      mkdir(filename.c_str(), 0700);

    // open the file
    filename += "/bloch_eigenvalues";
    filename += std::to_string(indx);
    filename += "-";
    filename += std::to_string(step);
    filename += "-";
    filename += std::to_string(cycle);

    std::ofstream outputFile;
    outputFile.open(filename.c_str());

    // output the "step" and the wave_ratio
    outputFile << step << std::endl;
    outputFile << wave_ratio << std::endl;

    // output the eigenvalues to the file
    for (unsigned int i = 0 ; i < eigenvalues.size(); i ++)
    {
      outputFile << std::setprecision(15) <<  eigenvalues[i].real() << std::endl;
    }

    outputFile.close();

    return eigenvalues[0].real();

  }

  void ElasticProblem::set_unstable_eigenvector_as_initial_tangent(unsigned int indx)
  {
    unsigned int number_dofs = dof_handler.n_dofs();

    evaluation_point = present_solution;
    assemble_system_matrix();
    apply_boundaries_and_constraints_system_matrix();

//    // copy current sparse system matrix to a full matrix.
//    LAPACKFullMatrix<double> system_matrix_full;
//    system_matrix_full.copy_from(system_matrix);
//
//    Vector<double> eigenvalues;
//    FullMatrix<double> eigenvectors;
//
//    // get the eigenvalues
//    system_matrix_full.compute_eigenvalues_symmetric(-10, 0.3, 1e-6, eigenvalues, eigenvectors);

    std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
    std::vector<double>                     eigenvalues;
    IndexSet eigenfunction_index_set(dof_handler.n_dofs());
    eigenfunction_index_set.add_range(0, dof_handler.n_dofs());
    eigenfunctions.resize (8);
    for (unsigned int i=0; i<eigenfunctions.size (); ++i)
      eigenfunctions[i].reinit (eigenfunction_index_set, MPI_COMM_WORLD);

    eigenvalues.resize (eigenfunctions.size());

    system_matrix_petsc = 0.0;
    for (unsigned int row = 0; row < system_matrix.m(); ++row)
    {
      const typename SparseMatrix<double>::const_iterator end_row = system_matrix.end(row);
      for (typename SparseMatrix<double>::const_iterator entry = system_matrix.begin(row);
                             entry != end_row; ++entry)
       {
         if(fabs(entry->value()) > 1e-8)
           system_matrix_petsc.set(row, entry->column(),entry->value());
       }
    }

    system_matrix_petsc.compress(VectorOperation::insert);
    SolverControl solver_control (dof_handler.n_dofs(), 1e-6);
    SLEPcWrappers::SolverKrylovSchur eigensolver (solver_control);
//    SLEPcWrappers::SolverLAPACK eigensolver (solver_control);

    // Before we actually solve for the eigenfunctions and -values, we have to
    // also select which set of eigenvalues to solve for. Lets select those
    // eigenvalues and corresponding eigenfunctions with the smallest real
    // part (in fact, the problem we solve here is symmetric and so the
    // eigenvalues are purely real). After that, we can actually let SLEPc do
    // its work:
//    eigensolver.set_target_eigenvalue(-10.0);
    eigensolver.set_which_eigenpairs (EPS_SMALLEST_REAL);

    eigensolver.set_problem_type (EPS_HEP);

    eigensolver.solve (system_matrix_petsc,
                       eigenvalues, eigenfunctions,
                       eigenfunctions.size());

    bool eigenvector_flag = false;
    unsigned int numNegVals = 0;
    for (unsigned int i = 0 ; i < eigenvalues.size(); i ++)
    {
      double nextEigenVal = eigenvalues[i];
      if (nextEigenVal < 0.0)
      {
        numNegVals ++;

        // set the eigenvector corresponnding the negative eigenvalue specified by index
        if (numNegVals == indx)
        {
          eigenvector_flag = true;

          initial_solution_tangent.reinit(dof_handler.n_dofs());
          for (unsigned int j = 0; j < dof_handler.n_dofs(); j++)
            initial_solution_tangent[j] = eigenfunctions[i][j];

          break;
        }
      }
    }

    if (eigenvector_flag == false)
    {
      std::cout << "Unable to set unstable eigenvetor because there was not " << indx << " negative eigenvalues. Exiting" << std::endl;
      exit(-1);
    }
    // So the eigenvector had bad entries for the constrained dofs, so need to distribute the constraints
    constraints.distribute(initial_solution_tangent);

    // also, normalize the tangent vector
    initial_solution_tangent *= (1.0/initial_solution_tangent.l2_norm());


    // Want to make sure that localization happens in the middle of the structure.
    // This can be accomplished by makeing sure the x_2 component of the initial tangent
    // is negatice at x_1 = 0, x_2 = 1
    std::vector<Point<DIM>> support_points(number_dofs);
    map_dofs_to_support_points(dof_handler, support_points);

    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    std::vector<bool> is_x2_comp(number_dofs);

    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);

    bool foundFlag = false;
    for(unsigned int i = 0; i < number_dofs; i++)
    {
      if ( (fabs(support_points[i](0)) < 1e-5) && (fabs(support_points[i](1) - 1.0) < 1e-5)
            && is_x2_comp[i] == true)
      {
        foundFlag = true;
        if (initial_solution_tangent[i] > 0.0)
          initial_solution_tangent *= -1.0;

        break;
      }
    }

    if (foundFlag == false)
    {
      std::cout << "Did not find a node at x_1 = 0, x_2 = 1. Exiting.\n";
      exit(-1);
    }


    // and finally, just set the lambda tangent to 0.0.
    initial_lambda_tangent = 0.0;
  }

  double ElasticProblem::bisect_find_lambda_critical(double lowerBound, double upperBound,
                                                          double tol, unsigned int maxIter)
  {
    if (lowerBound < 0.0)
      lowerBound = 0.0;

    unsigned int N = 1;
    double middleVal = 0.0;
    while ( N < maxIter)
    {
      middleVal = (upperBound + lowerBound)/2.0;

      if ((upperBound - lowerBound)/2 < tol)
        return middleVal;

      N += 1;

      if (get_system_eigenvalues(middleVal, -1) == get_system_eigenvalues(lowerBound, -1))
        lowerBound = middleVal;
      else
        upperBound = middleVal;
    }

    return middleVal;
  }

  void ElasticProblem::output_results (const unsigned int cycle)
  {

    std::vector<std::string> solution_names;
    switch (DIM)
      {
      case 1:
        solution_names.push_back ("displacement");
        break;
      case 2:
        solution_names.push_back ("x1_displacement");
        solution_names.push_back ("x2_displacement");
        break;
      case 3:
        solution_names.push_back ("x1_displacement");
        solution_names.push_back ("x2_displacement");
        solution_names.push_back ("x3_displacement");
        break;
      default:
        Assert (false, ExcNotImplemented());
        break;
      }

    // output the total displacements. this requires adding in the uniform solution on top of the displacements

    std::string filename0(output_directory);
    filename0 += "/total_displacement";

    // see if the directory exists...
    struct stat st;
    if (stat(filename0.c_str(), &st) == -1)
      mkdir(filename0.c_str(), 0700);

    filename0 += "/total_displacement-";
    filename0 += std::to_string(cycle);

    filename0 += ".vtk";
    std::ofstream output_totalDisp (filename0.c_str());

    DataOut<DIM, hp::DoFHandler<DIM>> data_out_totalDisp;

    data_out_totalDisp.attach_dof_handler (dof_handler);


    // Get the total displacement of each of the points.
    // Get the points of the dofs so we can do some shifting...
    std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
    map_dofs_to_support_points(dof_handler, support_points);

    const unsigned int   number_dofs = dof_handler.n_dofs();

    Vector<double> shifted_solution(number_dofs);

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> is_x1_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    for (unsigned int i = 0; i < number_dofs; i++)
    {
      if(is_x1_comp[i])
      {
        // it is an x1 component
        shifted_solution[i] = present_solution[i] + support_points[i](0)*(F0[0][0] - 1.0);
      }
      else
      {
        // it is an x2 component
        shifted_solution[i] = present_solution[i] + support_points[i](1)*(F0[1][1] - 1.0);
      }
    }

    data_out_totalDisp.add_data_vector (shifted_solution, solution_names);
    data_out_totalDisp.build_patches (max_degree);
    data_out_totalDisp.write_vtk (output_totalDisp);


    // Now output the displacements from uniform solution

    std::string filename1(output_directory);
    filename1 += "/displacement_from_uniform";

    // see if the directory exists...
    if (stat(filename1.c_str(), &st) == -1)
      mkdir(filename1.c_str(), 0700);

    filename1 += "/displacement_from_uniform-";
    filename1 += std::to_string(cycle);
    filename1 += ".vtk";
    std::ofstream output_disp_from_uniform (filename1.c_str());

    DataOut<DIM, hp::DoFHandler<DIM>> data_out_disp_from_uniform;

    data_out_disp_from_uniform.attach_dof_handler (dof_handler);

    data_out_disp_from_uniform.add_data_vector (present_solution, solution_names);
    data_out_disp_from_uniform.build_patches (max_degree);
    data_out_disp_from_uniform.write_vtk (output_disp_from_uniform);

    // Now output the deformed mesh

    // just need to shift the corrdinates of the verticies by the shifted solution vector

    DataOut<DIM, hp::DoFHandler<DIM>> deformed_data_out;

    deformed_data_out.attach_dof_handler(dof_handler);
    deformed_data_out.add_data_vector(shifted_solution, solution_names);
    hp::MappingCollection<DIM, DIM> mapping;
    MappingQEulerian_hp<DIM> q_mapping(1,  dof_handler, shifted_solution);
    deformed_data_out.build_patches(q_mapping, 1);

    std::string filename2(output_directory);
    filename2 += "/deformed_mesh";
    // see if the directory exists...
    if (stat(filename2.c_str(), &st) == -1)
      mkdir(filename2.c_str(), 0700);
    filename2 += "/deformed_mesh-";
    filename2 += std::to_string(cycle);
    filename2 += ".vtk";
    std::ofstream output_deformed_mesh(filename2.c_str());
    deformed_data_out.write_vtk(output_deformed_mesh);


  }

  void ElasticProblem::output_load_info(std::vector<double> lambda_values,
                                             std::vector<double> energy_values,
                                             std::vector<double> congugate_lambda_values,
                                             std::vector<double> displacement_magnitude,
                                             const unsigned int cycle,
                                             std::vector<bool> stable) const
  {

    // output the lambda value, system energy, and the congugate lambda vales for each step

    std::string filename(output_directory);
    filename += "/load_info";
    // see if the directory exists...
    struct stat st;
    if (stat(filename.c_str(), &st) == -1)
      mkdir(filename.c_str(), 0700);

    filename += "/load_info";
    filename += std::to_string(cycle);
    filename += ".txt";
    std::ofstream load_data_output;
    load_data_output.open(filename.c_str());
    load_data_output << "# lambda";
    load_data_output << std::setw(25) << "energy" ;
    load_data_output << std::setw(25) << "congugate_lambda" << std::endl;
    load_data_output << std::endl;
    for(unsigned int i = 0; i < lambda_values.size(); i ++)
    {
      load_data_output << std::setprecision(15) << std::setw(8) << lambda_values[i];
      load_data_output << std::setprecision(15) << std::setw(25) << energy_values[i];
      load_data_output << std::setprecision(15) << std::setw(25) << congugate_lambda_values[i];
      load_data_output << std::setprecision(15) << std::setw(25) << displacement_magnitude[i];
      if(stable.size() > 1)
        load_data_output << std::setw(25) << stable[i] << std::endl;
      else
        load_data_output << std::endl;
    }

    load_data_output.close();
  }

  void ElasticProblem::read_input_file(char* filename)
  {
    FILE* fid;
    int endOfFileFlag;
    char nextLine[MAXLINE];

    int valuesWritten;
    bool fileReadErrorFlag = false;

    grid_dimensions.resize(2);
    domain_dimensions.resize(2);

    fid = std::fopen(filename, "r");
    if (fid == NULL)
    {
      std::cout << "Unable to open file \"" << filename  << "\"" <<  std::endl;
      fileReadErrorFlag = true;
    }
    else
    {
      std::vector<unsigned int>::iterator last;

      // Read in the output name
      char directory_name[MAXLINE];
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%s", directory_name);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in the number of sections
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u", &number_sections);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in the polynomial degree of each section
      section_polynomial_degrees.clear();
      unsigned int next_polynomial_degree;
      unsigned int values_count; values_count = 0;
      char* nextValueChar;
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      nextValueChar = strtok(nextLine, ", \t");
      for(unsigned int i = 0; i < number_sections; i++)
      {
        valuesWritten = sscanf(nextValueChar, "%u", &next_polynomial_degree);
        if (valuesWritten != 1)
        {
          fileReadErrorFlag = true;
          goto fileClose;
        }

        section_polynomial_degrees.push_back(next_polynomial_degree);
        values_count++;

        nextValueChar = strtok(NULL, ", \t");
        if(nextValueChar == NULL)
          break;
      }
      if(values_count != number_sections)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      // now find unique values and stuff
      FE_id_polynomial_degree = section_polynomial_degrees;
      last = std::unique(FE_id_polynomial_degree.begin(), FE_id_polynomial_degree.end());
      FE_id_polynomial_degree.erase(last, FE_id_polynomial_degree.end());
      section_FE_id.clear();
      section_FE_id.resize(section_polynomial_degrees.size());
      for(unsigned int i = 0; i < section_polynomial_degrees.size(); i ++)
      {
        for(unsigned int j = 0; j < FE_id_polynomial_degree.size(); j++)
        {
          if(section_polynomial_degrees[i] == FE_id_polynomial_degree[j])
          {
            section_FE_id[i] = j;
            break;
          }
        }
      }
      max_degree = *std::max_element(FE_id_polynomial_degree.begin(), FE_id_polynomial_degree.end());

      // Read in the number of unit cells
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u", &number_unit_cells);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // append the number of unit cells used to the directory name
      char num_cells_str[MAXLINE];
      sprintf(num_cells_str, "%u", number_unit_cells);
      strcat(directory_name, "_");
      strcat(directory_name, num_cells_str);

      sprintf(output_directory, "output/");
      strcat(output_directory, directory_name);

      // Read in the grid dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u", &grid_dimensions[0], &grid_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      grid_dimensions[0] *= number_unit_cells;

      // read in the absolute tolerance of newton iteration
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg  %u", &tol, &maxIter);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
      }

      // read in exponential growth parameter and possibly the l1 value

      double val1;
      double val2;
      double val3;
      double nuVal;
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg", &val1, &val2, &val3);
      if(valuesWritten < 1 || (valuesWritten == 1 && val1 < 1.0))
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      else if ((valuesWritten == 1) && val1 > 1.0)
      {
        // treat as just a kappa value
        kappa = val1;
        L1 = 1.1;
        nuVal = NU_VALUE;

        pieceConstFlag = false;

        mu = new MuFunction(kappa);

        std::cout << "Using Exponential mu with: k = " << kappa << ", nu = " << nuVal << std::endl;
      }
      else if((valuesWritten == 2) && val1 < 1.0)
      {
        // treat as nu and a kappa
        nuVal = val1;
        kappa = val2;
        L1 = 1.1;

        pieceConstFlag = false;

        mu = new MuFunction(kappa);
        nu.set_nu_value(nuVal);

        std::cout << "Using Exponential mu with: k = " << kappa << ", nu = " << nuVal << std::endl;

      }
      else if((valuesWritten == 2) && val1 >= 1.0)
      {
        // treat as a kappa and a l1
        kappa = val1;
        L1 = val2;
        nuVal = NU_VALUE;

        pieceConstFlag = true;

        mu = new MuFunction(kappa, L1);

        std::cout << "Using Piecewise Constant mu with: k = " << kappa << ", L1 = " << L1 << ", nu = " << nuVal << std::endl;

      }
      else if(valuesWritten == 3)
      {
        // treat as a nu, kappa, and l1
        nuVal = val1;
        kappa = val2;
        L1 = val3;

        pieceConstFlag = true;

        mu = new MuFunction(kappa, L1);
        nu.set_nu_value(nuVal);

        std::cout << "Using Piecewise Constant mu with: k = " << kappa << ", L1 = " << L1 << ", nu = " << nuVal << std::endl;
      }


      // read in critical lambda
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg", &critical_lambda_analytical);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the critical frequency
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg", &critical_frequency);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the ds and load steps
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %u", &ds, &load_steps);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the output every
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u", &output_every);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      fileClose:
      {
        fclose(fid);
      }
    }

    if (fileReadErrorFlag)
    {
      // default parameter values
      std::cout << "Error reading input file, Exiting.\n" << std::endl;
      exit(1);
    }
    else
      std::cout << "Input file successfully read" << std::endl;

    // set the domain dimensions
    domain_dimensions[0] = number_unit_cells*2.0*(4.0*atan(1.0))/critical_frequency;
    domain_dimensions[1] = 1.0;

    // make the output directory
    struct stat st;
    if (stat("./output", &st) == -1)
       mkdir("./output", 0700);

    if (stat(output_directory, &st) == -1)
      mkdir(output_directory, 0700);


    init_fe_and_quad_collection();

  }

  void ElasticProblem::init_fe_and_quad_collection()
  {
    for(unsigned int i = 0; i < FE_id_polynomial_degree.size(); i++)
    {
      unsigned int next_degree = FE_id_polynomial_degree[i];
      fe_collection.push_back(FESystem<DIM>(FE_Q<DIM>(next_degree), DIM));
      if(next_degree == 1)
      {
        QGauss<1> quad_x(1);
        QGauss<1> quad_y(2);
        QAnisotropic<DIM> quadrature_formula_reduced(quad_x, quad_y);
        quadrature_collection.push_back(quadrature_formula_reduced);
      }
      else
        quadrature_collection.push_back(QGauss<DIM>(next_degree + 1));

    }
    quadrature_collection.push_back(QGauss<DIM>(10));

    last_q_index = FE_id_polynomial_degree.size();
  }

  void ElasticProblem::read_asymptotic_input_file(char* filename)
  {
    FILE* fid;
    int endOfFileFlag;
    char nextLine[MAXLINE];

    int valuesWritten;
    bool fileReadErrorFlag = false;
    unsigned int dummy;
    unsigned int dummy1;

    grid_dimensions.resize(2);
    domain_dimensions.resize(2);

    fid = std::fopen(filename, "r");
    if (fid == NULL)
    {
     std::cout << "Unable to open file \"" << filename  << "\"" <<  std::endl;
     fileReadErrorFlag = true;
    }
    else
    {
      double nuVal = 0.0;
      double kappa_new = 0.0;
      double L1_new = 0.0;


      // Read in the output name
      char directory_name[MAXLINE];
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%s", directory_name);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in the number of unit cells
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u", &dummy);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // Read in the grid dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u", &dummy, &dummy1);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }


      grid_dimensions[0] *= number_unit_cells;

      // read in exponential growth parameter and possibly the L1 value and possibly the nu value. Compare to what we already loaded too

      double val1;
      double val2;
      double val3;
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg", &val1, &val2, &val3);
      if(valuesWritten < 1 || (valuesWritten == 1 && val1 < 1.0))
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      else if ((valuesWritten == 1) && val1 >= 1.0)
      {
        // treat as just a kappa value
        kappa_new = val1;
        L1_new = 1.1;
        nuVal = NU_VALUE;

        std::cout << "Using Exponential mu with: k = " << kappa_new << ", nu = " << nuVal << std::endl;
      }
      else if((valuesWritten == 2) && val1 < 1.0)
      {
        // treat as nu and a kappa
        nuVal = val1;
        kappa_new = val2;
        L1_new = 1.1;

        std::cout << "Using Exponential mu with: k = " << kappa_new << ", nu = " << nuVal << std::endl;

      }
      else if((valuesWritten == 2) && val1 >= 1.0)
      {
        // treat as a kappa and a l1
        kappa_new = val1;
        L1_new = val2;
        nuVal = NU_VALUE;

        std::cout << "Using Piecewise Constant mu with: k = " << kappa_new << ", L1 = " << L1_new << ", nu = " << nuVal << std::endl;

      }
      else if(valuesWritten == 3)
      {
        // treat as a nu, kappa, and l1
        nuVal = val1;
        kappa_new = val2;
        L1_new = val3;

        std::cout << "Using Piecewise Constant mu with: k = " << kappa_new << ", L1 = " << L1_new << ", nu = " << nuVal << std::endl;
      }

      // check that everything matches
      if(!(kappa_new == kappa && L1_new == L1 && nuVal == nu.get_nu_value()))
      {
        std::cout << "Parameters for asymptotic input file do not match the original input file\n";
        fileReadErrorFlag = true;
        goto fileClose;
      }


      // read in critical lambda
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg", &critical_lambda_analytical);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the critical frequency
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg", &critical_frequency);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      w_c = critical_frequency;


      double re1, re2, im1, im2, re3, re4, im3, im4;

      // read in the roots to the characteristic, alphas
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
         &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      alphas.resize(4);
      alphas[0] = std::complex<double>(re1, im1);
      alphas[1] = std::complex<double>(re2, im2);
      alphas[2] = std::complex<double>(re3, im3);
      alphas[3] = std::complex<double>(re4, im4);

      double re5, re6, im5, im6, re7, re8, im7, im8;

      // read in the amplitude A's
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
         &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4, &re5, &im5, &re6, &im6, &re7, &im7, &re8, &im8);
      if((valuesWritten != 8 && valuesWritten != 16)
         || (pieceConstFlag == (valuesWritten == 8)) )
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      else if(valuesWritten == 16)
      {
        A.resize(8);
        A[4] = std::complex<double>(re5, im5);
        A[5] = std::complex<double>(re6, im6);
        A[6] = std::complex<double>(re7, im7);
        A[7] = std::complex<double>(re8, im8);
      }
      else
        A.resize(4);

      A[0] = std::complex<double>(re1, im1);
      A[1] = std::complex<double>(re2, im2);
      A[2] = std::complex<double>(re3, im3);
      A[3] = std::complex<double>(re4, im4);

      // read in the B's
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
         &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      B.resize(4);
      B[0] = std::complex<double>(re1, im1);
      B[1] = std::complex<double>(re2, im2);
      B[2] = std::complex<double>(re3, im3);
      B[3] = std::complex<double>(re4, im4);

      fileClose:
      {
       fclose(fid);
      }
    }

    if (fileReadErrorFlag)
    {
     // default parameter values
     std::cout << "Error reading input file, Exiting.\n" << std::endl;
     exit(1);
    }
    else
     std::cout << "Input file successfully read" << std::endl;
  }

  void ElasticProblem::getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                          int const maxSize, int* const endOfFileFlag)
  {
    *endOfFileFlag = 0;
    do
    {
      if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
      {
        *endOfFileFlag = 1;
        break;
      }
      while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
             (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
      {
        nextLinePtr = (nextLinePtr + 1);
      }
    }
    while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
  }

  void ElasticProblem::save_current_state(unsigned int indx)
  {
    // create the output directory if it doesnt exist

    char saved_state_dir[MAXLINE];
    strcpy(saved_state_dir, output_directory);
    strcat(saved_state_dir, "/saved_state");

    // see if the directory exists
    struct stat st;
    if (stat(saved_state_dir, &st) == -1)
         mkdir(saved_state_dir, 0700);


    char index_char[32];
    sprintf(index_char, "%u", indx);

    // Triangulation
    char triag_file[MAXLINE];
    strcpy(triag_file, saved_state_dir);
    strcat(triag_file, "/triag_");
    strcat(triag_file, index_char);
    strcat(triag_file, ".dat");
    std::ofstream triag_out (triag_file);
    boost::archive::text_oarchive triag_ar(triag_out);
    triangulation.save(triag_ar, 1);

    // dof handler
    char dof_file[MAXLINE];
    strcpy(dof_file, saved_state_dir);
    strcat(dof_file, "/dof_");
    strcat(dof_file, index_char);
    strcat(dof_file, ".dat");
    std::ofstream dof_out (dof_file);
    boost::archive::text_oarchive dof_ar(dof_out);
    dof_handler.save(dof_ar, 1);

    // Present Solution
    char present_solution_file[MAXLINE];
    strcpy(present_solution_file, saved_state_dir);
    strcat(present_solution_file, "/present_solution_");
    strcat(present_solution_file, index_char);
    strcat(present_solution_file, ".dat");
    std::ofstream solution_out (present_solution_file);
    boost::archive::text_oarchive solution_ar(solution_out);
    present_solution.save(solution_ar, 1);

    char lambda_file[MAXLINE];
    strcpy(lambda_file, saved_state_dir);
    strcat(lambda_file, "/present_lambda.dat");
    std::ofstream lambda_out(lambda_file);

    lambda_out << present_lambda;
    lambda_out.close();

  }

  void ElasticProblem::load_state(unsigned int indx)
  {
    // create the output directory

    char input_dir_path[MAXLINE];
    strcpy(input_dir_path, output_directory);
    strcat(input_dir_path, "/saved_state");
    struct stat st;
    if (stat(input_dir_path, &st) == -1)
    {
      std::cout << "Could not find the directory : " << input_dir_path << "\nExiting." <<std::endl;
      exit(-1);
    }

    char index_char[32];
    sprintf(index_char, "%u", indx);

    // Triangulation
    char triag_file[MAXLINE];
    strcpy(triag_file, input_dir_path);
    strcat(triag_file, "/triag_");
    strcat(triag_file, index_char);
    strcat(triag_file, ".dat");
    std::ifstream triag_in (triag_file);
    boost::archive::text_iarchive triag_ar(triag_in);
    triangulation.load(triag_ar, 1);

    // df_handler
    dof_handler.distribute_dofs(fe_collection);
    char dof_file[MAXLINE];
    strcpy(dof_file, input_dir_path);
    strcat(dof_file, "/dof_");
    strcat(dof_file, index_char);
    strcat(dof_file, ".dat");
    std::ifstream dof_in (dof_file);
    boost::archive::text_iarchive dof_ar(dof_in);
    dof_handler.load(dof_ar, 1);

    // Present Solution
    char present_solution_file[MAXLINE];
    strcpy(present_solution_file, input_dir_path);
    strcat(present_solution_file, "/present_solution_");
    strcat(present_solution_file, index_char);
    strcat(present_solution_file, ".dat");
    std::ifstream solution_in (present_solution_file);
    boost::archive::text_iarchive solution_ar(solution_in);
    present_solution.load(solution_ar, 1);
    evaluation_point = present_solution;

    // lambda value
    char lambda_file[MAXLINE];
    strcpy(lambda_file, input_dir_path);
    strcat(lambda_file, "/present_lambda.dat");
    std::ifstream lambda_in(lambda_file);
    if (!lambda_in)
    {
      std::cout << "Error reading file: " << lambda_file << " . Exiting." << std::endl;
      exit(-1);
    }

    lambda_in >> std::setprecision(15) >> present_lambda;
    lambda_in.close();

    update_F0(present_lambda);

    fileLoadFlag = true;
  }

  void ElasticProblem::load_intial_tangent()
  {
    char filename[MAXLINE];
    std::cout << "Please enter file that tangent vector is loaded from : " << std::endl;
    std::cin >> filename;
    std::ifstream tangent_in (filename);
    boost::archive::text_iarchive tangent_ar(tangent_in);
    initial_solution_tangent.load(tangent_ar, 1);

    std::cout << "Please Enter the lambda tangent : " << std::endl;
    std::cin >> initial_lambda_tangent;

    // make sure it has length 1
    double scalingVal = initial_solution_tangent.norm_sqr() + initial_lambda_tangent*initial_lambda_tangent;
    scalingVal = 1.0/sqrt(scalingVal);

    initial_solution_tangent *= scalingVal;
    initial_lambda_tangent *= scalingVal;


  }

  void ElasticProblem::renumber_boundary_ids()
  {

    // renumber boundary ids because they have problems being saved for nonuniform mesh.
    typename Triangulation<DIM>::active_cell_iterator cell =
     triangulation.begin_active(), endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int f = 0; f < GeometryInfo<DIM>::faces_per_cell; ++f)
      {

        const Point<DIM> face_center = cell->face(f)->center();
        if (fabs(face_center[0] + domain_dimensions[0]/2.0) < 1e-4)  // left
        {
          cell->face(f)->set_boundary_id (1);
        }
        else if (fabs(face_center[0] - domain_dimensions[0]/2.0) < 1e-4)
        {// right
          cell->face(f)->set_boundary_id (2);
        }
        else if (fabs(face_center[1]) < 1e-6) //bottom
        {
          cell->face(f)->set_boundary_id (3);
        }
        else if (fabs(face_center[1] - 1.0) < 1e-6) //top
        {
         cell->face(f)->set_boundary_id (4);

        }

      }
  }

  void ElasticProblem::print_dof_coords_and_vals(unsigned int indx)
  {
    unsigned int number_dofs = dof_handler.n_dofs();

    std::string fileName = "dofsAndVals_";
    fileName += std::to_string(indx);
    fileName += ".dat";
    std::ofstream dof_out;
    dof_out.open(fileName);

    // get the coords of the dofs
    std::vector<Point<DIM>> support_points(number_dofs);

    map_dofs_to_support_points(dof_handler, support_points);

    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    std::vector<bool> is_x2_comp(number_dofs);

    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);

    for(unsigned int i = 0 ; i < number_dofs; i ++)
    {
      dof_out << "(" << support_points[i](0) << ", " << support_points[i](1) << ")   "  << is_x2_comp[i] << "     " << evaluation_point[i] << std::endl;
    }

    dof_out.close();

  }


  void
  ElasticProblem::map_dofs_to_support_points (const hp::DoFHandler<DIM, DIM>  &dof_handler, std::vector<Point<DIM> > &support_points)
  {
    const unsigned int dim = DIM;
    const unsigned int spacedim = DIM;

    hp::FECollection<dim, spacedim> fe_collection(dof_handler.get_fe());
    hp::QCollection<dim> q_coll_dummy;

    for (unsigned int fe_index = 0; fe_index < fe_collection.size(); ++fe_index)
      {
        // check whether every fe in the collection has support points
        Assert(fe_collection[fe_index].has_support_points(),
               typename FiniteElement<dim>::ExcFEHasNoSupportPoints());
        q_coll_dummy.push_back(
          Quadrature<dim> (
            fe_collection[fe_index].get_unit_support_points()));
      }

    // Now loop over all cells and enquire the support points on each
    // of these. we use dummy quadrature formulas where the quadrature
    // points are located at the unit support points to enquire the
    // location of the support points in real space.
    //
    // The weights of the quadrature rule have been set to invalid
    // values by the used constructor.
    MappingQ1<DIM> mapping;
    hp::MappingCollection<DIM, DIM> mapping1(mapping);
    hp::FEValues<dim, spacedim> hp_fe_values(mapping1, fe_collection,
                                             q_coll_dummy, update_quadrature_points);
    typename hp::DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof_handler.begin_active(), endc = dof_handler.end();

    std::vector<types::global_dof_index> local_dof_indices;
    for (; cell != endc; ++cell)
      // only work on locally relevant cells
      if (cell->is_artificial() == false)
        {
          hp_fe_values.reinit(cell);
          const FEValues<dim, spacedim> &fe_values = hp_fe_values.get_present_fe_values();

          local_dof_indices.resize(cell->get_fe().dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);

          const std::vector<Point<spacedim> > &points =
            fe_values.get_quadrature_points();
          for (unsigned int i = 0; i < cell->get_fe().dofs_per_cell; ++i)
            // insert the values into the map
            support_points[local_dof_indices[i]] = points[i];
        }
  }

}





#endif // COMPRESSEDSTRIPPACABLOCH_CC_
