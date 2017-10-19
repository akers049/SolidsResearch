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



#define MU_VALUE 1.0
#define NU_VALUE 0.33



namespace NeoHookean_Newton
{
  using namespace dealii;


  /****************************************************************
                       Function Definitions
  ****************************************************************/

  template <int dim>
  double NuFunction<dim>::value (const Point<dim>  &p,
                                 const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0, ExcNotImplemented());

    // Put your function for nu. p(0) is x1 value, p(1) is the x2 value

    double nuValue = NU_VALUE + 0.0*(p(1) - 0.5);

    return nuValue;
  }

  template <int dim>
  void NuFunction<dim>::value_list(const std::vector< Point< dim > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = NuFunction<dim>::value(points[i], component);

  }


  template <int dim>
  double MuFunction<dim>::value (const Point<dim>  &p,
                                 const unsigned int  component) const
  {
    // this is a scalar function, so make sure component is zero...
    Assert (component == 0, ExcNotImplemented());

    // Put your function for mu. p(0) is x1 value, p(1) is the x2 value

    double muValue = mu0*exp(kappa*p(1));

    return muValue;
  }

  template <int dim>
  void MuFunction<dim>::value_list(const std::vector< Point< dim > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = MuFunction<dim>::value(points[i], component);
  }



  // computes right hand side values if we were to have body forces. But it just
  // always returns zeros because we don't.
  template <int dim>
  void right_hand_side (const std::vector<Point<dim> > &points,
                        std::vector<Tensor<1, dim> >   &values)
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));
    Assert (dim >= 2, ExcNotImplemented());


    // not imposing body forces or tractions
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      values[point_n][0] = 0.0;
      values[point_n][1] = 0.0;
    }
  }

  // Functions that get the tensors used in the stiffness and rhs calculations.

  template <int dim>
  inline
  Tensor<2,dim>
  ElasticProblem<dim>::get_deformation_gradient(std::vector<Tensor<1,dim> > old_solution_gradient)
  {
    Tensor<2,dim> tmp = F0;

    for (unsigned int i = 0; i < dim; i ++)
      for(unsigned int j = 0; j < dim; j++)
      {
        tmp[i][j] += old_solution_gradient[i][j];
      }

    return tmp;
  }


  template <int dim>
  inline
  Tensor<4,dim>
  ElasticProblem<dim>::get_incremental_moduli_tensor(const double nu,
                            const double mu,
                            Tensor<2,dim> F_inv, double II_F)
  {
    Tensor<4,dim> tmp;

    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        for (unsigned int k=0; k<dim; ++k)
          for (unsigned int l=0; l<dim; ++l)
          {
            tmp[i][j][k][l] = ((i==k) && (j==l) ? mu : 0.0) +
                (mu - (2.0*mu*nu/(1.0-nu))*(II_F*II_F - II_F))*F_inv[j][k]*F_inv[l][i] +
                (4.0*nu*mu/(1.0-nu))*(II_F*II_F - 0.5*II_F)*F_inv[l][k]*F_inv[j][i];
          }

    return tmp;
  }

  template <int dim>
  inline
  Tensor<2,dim>
  ElasticProblem<dim>::get_piola_kirchoff_tensor(const double nu, const double mu,
                                                 Tensor<2, dim> F, Tensor<2,dim> F_inv,
                                                 double II_F)
  {
    Tensor<2, dim> tmp;

    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
      {
        tmp[i][j] = mu*F[i][j] - mu*F_inv[j][i] +
                    (2.0*mu*nu/(1.0- nu))*(II_F*II_F - II_F)*F_inv[j][i];
      }

    return tmp;

  }

  template <int dim>
  inline
  double ElasticProblem<dim>::get_energy(const double nu, const double mu,
      Tensor<2, dim> F, double II_F)
  {
    double I_C = F[1][0]*F[1][0] + F[0][0]*F[0][0] + F[0][1]*F[0][1] + F[1][1]*F[1][1];
    double W = mu*(0.5*(I_C - 2 - log(II_F*II_F)) + (nu/(1.0- nu))*(II_F - 1.0)*(II_F - 1.0));

    return W;
  }


  template <int dim>
  ElasticProblem<dim>::ElasticProblem ()
    :
    dof_handler (triangulation),
    fe (FE_Q<dim>(1), dim)
  {}




  template <int dim>
  ElasticProblem<dim>::~ElasticProblem ()
  {
    dof_handler.clear ();
    delete mu;
  }


  template <int dim>
  void ElasticProblem<dim>::create_mesh()
  {

    // set domain dimensions
    domain_dimensions[0] = number_unit_cells*2.0*(4.0*atan(1.0))/critical_frequency;
    domain_dimensions[1] = 1.0;

    // creates our strip.
    Point<dim> corner1, corner2;
    corner1(0) = -domain_dimensions[0]/2.0;
    corner1(1) =  0.0;
    corner2(0) =  domain_dimensions[0]/2.0;
    corner2(1) =  domain_dimensions[1];
    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);


    /*
    // Now we will refine this mesh
    const int numSections = 2;
    for (int i = 0 ; i < numSections; i ++)
    {
      double section_x2 = (0.5/numSections)*i + 0.5;
      Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
      Triangulation<2>::active_cell_iterator endc = triangulation.end();
      for (; cell!=endc; ++cell)
      {
        for (unsigned int v=0;
             v < GeometryInfo<2>::vertices_per_cell;
             ++v)
          {
            const double x2_pos = (cell->vertex(v))(1);
            if ( x2_pos > section_x2)
              {
                cell->set_refine_flag ();
                break;
              }
          }
      }
      triangulation.execute_coarsening_and_refinement ();
    }*/
  }

  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.

    dof_handler.distribute_dofs (fe);

    // only reinit the present solution if it wasn't already loaded in
    if (present_solution.size() == 0)
      present_solution.reinit (dof_handler.n_dofs());

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    setup_system_constraints();

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    drhs_dlambda.reinit (dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

    update_bloch_wave_constraints(0.1);

    setup_bloch();

    evaluation_point = present_solution;

  }

  template <int dim>
  void ElasticProblem<dim>::setup_system_constraints ()
  {
    constraints.clear ();

    // periodic in x1 on 0 and 1 face (x1 faces) will do the x2 direction with x2 symmetry...
    DoFTools::make_periodicity_constraints<DoFHandler<dim>>(dof_handler, 0, 1, 0, constraints);

    // now do constraints that the average x1 displacement is zero

    const unsigned int   number_dofs = dof_handler.n_dofs();

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> boundary_dof_x1 (number_dofs, false);
    std::vector<bool> boundary_01_dof_x1 (number_dofs, false);

    std::set< types::boundary_id > boundary_id_01;
    boundary_id_01.insert(0);
    boundary_id_01.insert(1);

    DoFTools::extract_boundary_dofs(dof_handler,
                                    x1_mask,
                                    boundary_dof_x1);

    DoFTools::extract_boundary_dofs (dof_handler,
                                     x1_mask,
                                     boundary_01_dof_x1,
                                     boundary_id_01);



    unsigned int first_boundary_dof = 0;
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if ((boundary_dof_x1[i] == true) && (boundary_01_dof_x1[i] == false))
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

    // now do the constraint that the x2 displacements are symmetric about the line x1 = 0
    // and the x1 displacements are antisymmetric

    // get the coords of the dofs
    std::vector<Point<dim>> support_points(number_dofs);
    MappingQ1<dim> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

    // get the vector that tells us it the dof is an x2 component
    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    std::vector<bool> is_x2_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);

    // get the vector that tells us if it is a boundary 2 x2 dof
    // (these are already constrained to zero)

    std::vector<bool> boundary_2_dof_x2 (number_dofs, false);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);

    DoFTools::extract_boundary_dofs(dof_handler,
                                       x2_mask,
                                       boundary_2_dof_x2,
                                       boundary_id_2);



    matched_dofs.resize(number_dofs, -1);
    for(unsigned int i = 0; i < number_dofs; i++)
    {
      if (is_x2_comp[i])
      {
        if ( boundary_2_dof_x2[i] || matched_dofs[i] != -1 || (support_points[i](0) == 0.0))
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
        if (support_points[i](0) == 0.0)
        {
          constraints.add_line(i);
        }
        else if ( boundary_01_dof_x1[i] || matched_dofs[i] != -1 || (support_points[i](0) == 0.0))
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

  template<int dim>
  void ElasticProblem<dim>::setup_bloch()
  {
    // manually do the sparsity stuff for the double sized matrix
    const unsigned int number_dofs = dof_handler.n_dofs();
    DynamicSparsityPattern dsp_bloch(2*number_dofs, 2*number_dofs);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell->get_dof_indices (local_dof_indices);
      constraints_bloch.add_entries_local_to_global (local_dof_indices,
                                                            dsp_bloch,
                                                           /*keep_constrained_dofs*/ true);
      for(unsigned int i = 0; i < dofs_per_cell; i++)
        local_dof_indices[i] += number_dofs;


      constraints_bloch.add_entries_local_to_global (local_dof_indices,
                                                            dsp_bloch,
                                                           /*keep_constrained_dofs*/ true);
    }

    sparsity_pattern_bloch.copy_from(dsp_bloch);

    bloch_matrix.reinit(sparsity_pattern_bloch);
  }

  template<int dim>
  void ElasticProblem<dim>::update_bloch_wave_constraints(double wave_ratio)
  {
   constraints_bloch.clear();

   double k_x1 = 2.0*(4.0*atan(1.0))*wave_ratio;

   // get the boundary 0 and boundary 1 dofs
   const unsigned int   number_dofs = dof_handler.n_dofs();

   std::vector<bool> x1_components = {true, false};
   ComponentMask x1_mask(x1_components);

   std::vector<bool> boundary_1_dof (number_dofs, false);
   std::vector<bool> boundary_0_dof (number_dofs, false);

   std::set< types::boundary_id > boundary_id_0;
   boundary_id_0.insert(0);

   std::set< types::boundary_id > boundary_id_1;
   boundary_id_1.insert(1);

   std::vector<bool> x1_x2_components = {true, true};
   ComponentMask x1_x2_mask(x1_x2_components);

   DoFTools::extract_boundary_dofs(dof_handler,
                                   x1_x2_mask,
                                   boundary_0_dof,
                                   boundary_id_0);

   DoFTools::extract_boundary_dofs (dof_handler,
                                    x1_x2_mask,
                                    boundary_1_dof,
                                    boundary_id_1);

   std::vector<bool> is_x1_comp(number_dofs, false);

   DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

   // get the coords of the dofs
   std::vector<Point<dim>> support_points(number_dofs);
   MappingQ1<dim> mapping;
   DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

   for(unsigned int i = 0; i < number_dofs; i ++)
   {
     if (boundary_0_dof[i])
     {
       double x2_coord = support_points[i](1);
       for (unsigned int j = 0; j < number_dofs; j++)
       {
         if (boundary_1_dof[j] &&
             (fabs(x2_coord - support_points[j](1)) < 1e-12 ) &&
             (is_x1_comp[i] == is_x1_comp[j]) )
         {

           constraints_bloch.add_line (i);
           constraints_bloch.add_entry (i, j, cos(k_x1)); //-cos(k_x1));
           constraints_bloch.add_entry (i, (j + number_dofs) , sin(k_x1)); //-cos(k_x1));

           constraints_bloch.add_line (i + number_dofs);
           constraints_bloch.add_entry ((i + number_dofs), j, -sin(k_x1));
           constraints_bloch.add_entry ((i + number_dofs), (j + number_dofs), cos(k_x1));

           break;
         }
       }
     }
   }

   constraints_bloch.close();

  //  constraints.merge(hanging_node_constraints, ConstraintMatrix::MergeConflictBehavior::right_object_wins);
  }

  template<int dim>
  void ElasticProblem<dim>::set_boundary_values()
  {
    // this sets the boundary values of the solution vector so that the Newton step
    // can use homogeneous direchlet conditions on the set boundaries. It also makes sure
    // that the periodic faces' DoFs start with the same values (zero).

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::map<types::global_dof_index,double> boundary_values;

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 0,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 1,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 2,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 3,
                                                 ZeroFunction<dim>(dim),
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

  template<int dim>
  void ElasticProblem<dim>::update_F0(const double lambda)
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


  template <int dim>
  void ElasticProblem<dim>::add_small_pertubations(double amplitude, bool firstTime)
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

  template <int dim>
  void ElasticProblem<dim>::assemble_system_matrix()
  {
    // Assembling the system matrix. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_matrix = 0.0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_matrix = 0.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        Tensor<2,dim> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, dim> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<4,dim> d2W_dFdF = get_incremental_moduli_tensor(nu_values[q_point],
                                           mu_values[q_point], F_inv, II_F);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {
            const unsigned int component_m = fe.system_to_component_index(m).first;

            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int l=0; l<dim; ++l)
              {
                cell_matrix(n,m) += d2W_dFdF[component_n][j][component_m][l]*
                        fe_values.shape_grad(n, q_point)[j]*fe_values.shape_grad(m, q_point)[l]*fe_values.JxW(q_point);
              }
          }
        }
      }

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

  template <int dim>
  void ElasticProblem<dim>::assemble_system_rhs()
  {
    // Assembling the system rhs. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_rhs = 0.0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    std::vector<Tensor<1, dim> > rhs_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);
      right_hand_side (fe_values.get_quadrature_points(), rhs_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        Tensor<2,dim> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, dim> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<2,dim> dW_dF = get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point],
                                                        F, F_inv, II_F);
        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          for(unsigned int j = 0; j<dim; ++j)
          {
            cell_rhs(n) -= dW_dF[component_n][j]*fe_values.shape_grad(n, q_point)[j]*fe_values.JxW(q_point);
          }

          cell_rhs(n) += fe_values.shape_value(n, q_point)*rhs_values[q_point][component_n]*fe_values.JxW(q_point);
        }
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        system_rhs(local_dof_indices[n]) += cell_rhs(n);

    }

    constraints.condense (system_rhs);

  }

  template <int dim>
  void ElasticProblem<dim>::assemble_system_energy_and_congugate_lambda()
  {
    double lambda_eval = present_lambda;

    system_energy = 0.0;
    congugate_lambda = 0.0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();


    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points,
                                                std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    std::vector<Tensor<1, dim> > rhs_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,dim> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, dim> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<2,dim> dW_dF = get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point],
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

        double W = get_energy(nu_values[q_point], mu_values[q_point], F, II_F);
        system_energy += W*fe_values.JxW(q_point);
      }
    }

  }

  template <int dim>
  void ElasticProblem<dim>::assemble_drhs_dlambda()
  {
    double lambda_eval = present_lambda;
    update_F0(lambda_eval);

    drhs_dlambda = 0.0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    Vector<double>       cell_drhs_dlambda (dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points,
                                                std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();


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

    for (; cell!=endc; ++cell)
    {
      cell_drhs_dlambda = 0.0;
      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,dim> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, dim> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<4,dim> d2W_dFdF = get_incremental_moduli_tensor(nu_values[q_point],
                                           mu_values[q_point], F_inv, II_F);

        Tensor<2,dim> dW_dF_dlambda;
        for (unsigned int i = 0 ; i < dim; i ++)
          for (unsigned int j = 0 ; j < dim; j ++)
              dW_dF_dlambda[i][j] = d2W_dFdF[i][j][0][0]*dlambda1_dlambda  +  d2W_dFdF[i][j][1][1]*dlambda2_dlambda;



        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          for(unsigned int j = 0; j<dim; ++j)
          {
            cell_drhs_dlambda(n) += dW_dF_dlambda[component_n][j]*fe_values.shape_grad(n, q_point)[j]*fe_values.JxW(q_point);
          }
        }
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        drhs_dlambda(local_dof_indices[n]) += cell_drhs_dlambda(n);
    }

   constraints.condense (drhs_dlambda);
  }

  template<int dim>
  void ElasticProblem<dim>::assemble_bloch_matrix()
  {
    bloch_matrix = 0.0;

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

  template <int dim>
  void ElasticProblem<dim>::apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs)
  {
    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      if ((*homogenous_dirichlet_dofs)[i] == true)
        (*rhs)[i] = 0.0;
    }
  }

  template <int dim>
  void ElasticProblem<dim>::apply_boundaries_and_constraints_system_matrix()
  {
    constraints.condense (system_matrix);

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side2_components = {false, true};
    ComponentMask side2_mask(side2_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
                                              ZeroFunction<dim, double>(dim),
                                              boundary_values,
                                              side2_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);

  }

  template <int dim>
  void ElasticProblem<dim>::apply_boundaries_and_constraints_bloch_matrix()
  {
    constraints_bloch.condense(bloch_matrix);

    const unsigned int number_dofs = dof_handler.n_dofs();

    std::vector<bool> side2_components = {false, true};
    ComponentMask x2_mask(side2_components);
    std::vector<bool> boundary_2_dof_x2 (number_dofs, false);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);

    DoFTools::extract_boundary_dofs(dof_handler,
                                      x2_mask,
                                      boundary_2_dof_x2,
                                      boundary_id_2);

    unsigned int m = bloch_matrix.m();
    // set values on the diagonal to the first diagonal element,
    // or 1 if it is nonexistent
    double first_nonzero_diagonal_entry = 1.0;
    for (unsigned int i=0; i<m; ++i)
    {
      if (bloch_matrix.diag_element(i) != 0.0)
        {
          first_nonzero_diagonal_entry = bloch_matrix.diag_element(i);
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
        if((boundary_2_dof_x2[row%number_dofs] || boundary_2_dof_x2[entry->column()%number_dofs])
            && (row != entry->column()))
        {
          entry->value() = 0.0;
        }
        else if(boundary_2_dof_x2[row%number_dofs] && (row == entry->column()) &&
                (bloch_matrix.diag_element(row%number_dofs) == 0.0))
        {
          entry->value() = first_nonzero_diagonal_entry;
        }

      }

    }

   }

  template <int dim>
  void ElasticProblem<dim>::newton_iterate()
  {
    /* This function executes the newton iterations until it converges to a solution
     * or it exceeds the maximum number of iterations.
     */

    double current_residual;
    unsigned int iteration = 0;

    // get the dofs that we will apply dirichlet condition to
    std::vector<bool> boundary_2_dof_x2 (dof_handler.n_dofs(), false);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);


    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    DoFTools::extract_boundary_dofs(dof_handler,
                                       x2_mask,
                                       boundary_2_dof_x2,
                                       boundary_id_2);

    // Starts by getting the residual in the current configuration.
    evaluation_point = present_solution;

    assemble_system_rhs();
    apply_boundaries_to_rhs(&system_rhs, &boundary_2_dof_x2);

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
      line_search_and_add_step_length(current_residual, &boundary_2_dof_x2);

      evaluation_point = present_solution;
      current_residual = system_rhs.l2_norm();
      iteration ++;
    }

    // output iterations for convergance.
    std::cout << "    Converging Iterations : " << iteration << "\n";

  }

  template<int dim>
  void ElasticProblem<dim>::path_follow_PACA_iterate(Vector<double> *solVectorDir,
                                                          double lambdaDir, double ds)
  {
    // The setup stuff, may need to take  in more inputs. But for now we will do with this

    Vector<double> previousSolution = present_solution;
    double previousLambda = present_lambda;
    double current_residual = 0.0;
    unsigned int iteration = 0;

    // Scale the input tangent to unity just incase it wasnt already
//    double scalingVal = solVectorDir->norm_sqr() + lambdaDir*lambdaDir;
//    scalingVal = 1.0/sqrt(scalingVal);
//    lambdaDir *= scalingVal;
//    (*solVectorDir) *= scalingVal;

    // get the dofs that we will apply dirichlet condition to
    std::vector<bool> boundary_2_dof_x2 (dof_handler.n_dofs(), false);

    std::set< types::boundary_id > boundary_id_2;
    boundary_id_2.insert(2);


    std::vector<bool> x2_components = {false, true};
    ComponentMask x2_mask(x2_components);

    DoFTools::extract_boundary_dofs(dof_handler,
                                       x2_mask,
                                       boundary_2_dof_x2,
                                       boundary_id_2);

    // Starts by getting the residual in the initial guess.
    present_lambda += lambdaDir*ds;
    present_solution.add(ds, *solVectorDir);
    evaluation_point = present_solution;

    update_F0(present_lambda);

    assemble_system_rhs();
    apply_boundaries_to_rhs(&system_rhs, &boundary_2_dof_x2);

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
      apply_boundaries_to_rhs(&system_rhs, &boundary_2_dof_x2);

      // Assemble the stiffness matrix
      assemble_system_matrix();
      apply_boundaries_and_constraints_system_matrix();

      // Assemble the drhs_dlambda
      assemble_drhs_dlambda();
      apply_boundaries_to_rhs(&drhs_dlambda, &boundary_2_dof_x2);

      // Get the solution diff (bottom boarder of system mat)
      solution_diff = present_solution;
      solution_diff -= previousSolution;
      apply_boundaries_to_rhs(&solution_diff, &boundary_2_dof_x2);

      // get the bottom corner of system mat
      lambda_diff = present_lambda - previousLambda;

      // get bottom element of the rhs
      rhs_bottom = -0.5*( solution_diff.norm_sqr() + lambda_diff*lambda_diff - ds*ds);

      // Solve it with the boarder matrix solver
      solve_boarder_matrix_system();

      // add the update using the line search
      line_search_and_add_step_length_PACA(current_residual, &boundary_2_dof_x2, &previousSolution, previousLambda, ds);

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

  template<int dim>
  void ElasticProblem<dim>::line_search_and_add_step_length_PACA(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs,
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

  template<int dim>
  void ElasticProblem<dim>::line_search_and_add_step_length(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs)
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

  template <int dim>
  void ElasticProblem<dim>::solve ()
  {


    // direct solver for the system

    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);
    constraints.distribute (newton_update);

  }

  template <int dim>
  void ElasticProblem<dim>::solve_boarder_matrix_system()
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

  template<int dim>
  unsigned int ElasticProblem<dim>::get_system_eigenvalues(double lambda_eval, const int cycle)
  {

    update_F0(lambda_eval);
    evaluation_point = present_solution;
    assemble_system_matrix();
    apply_boundaries_and_constraints_system_matrix();

    // copy current sparse system matrix to a full matrix.
    LAPACKFullMatrix<double> system_matrix_full;
    system_matrix_full.copy_from(system_matrix);

    Vector<double> eigenvalues;
    FullMatrix<double> eigenvectors;

    system_matrix_full.compute_eigenvalues_symmetric(-10, 0.3, 1e-6, eigenvalues, eigenvectors);


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

  template<int dim>
  void ElasticProblem<dim>::get_bloch_eigenvalues(const int cycle, const int step, double wave_ratio)
  {

    evaluation_point = present_solution;
    update_F0(present_lambda);
    assemble_system_matrix();
    assemble_bloch_matrix();

    update_bloch_wave_constraints(wave_ratio);

    apply_boundaries_and_constraints_bloch_matrix();

    // copy current sparse system matrix to a full matrix.

    LAPACKFullMatrix<double> bloch_matrix_full;
    bloch_matrix_full.copy_from(bloch_matrix);


    Vector<double> eigenvalues;
    FullMatrix<double> eigenvectors;

    bloch_matrix_full.compute_eigenvalues_symmetric(-10, 0.3, 1e-6, eigenvalues, eigenvectors);

    std::string filename(output_directory);
    filename += "/bloch_eigenvalues";

    // see if the directory exists
    struct stat st;
    if (stat(filename.c_str(), &st) == -1)
      mkdir(filename.c_str(), 0700);

    // open the file
    filename += "/bloch_eigenvalues-";
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
     outputFile << std::setprecision(15) <<  eigenvalues[i] << std::endl;
    }

    outputFile.close();

  }

  template <int dim>
  void ElasticProblem<dim>::set_unstable_eigenvector_as_initial_tangent(unsigned int index)
  {

    evaluation_point = present_solution;
    assemble_system_matrix();
    apply_boundaries_and_constraints_system_matrix();

    // copy current sparse system matrix to a full matrix.
    LAPACKFullMatrix<double> system_matrix_full;
    system_matrix_full.copy_from(system_matrix);

    Vector<double> eigenvalues;
    FullMatrix<double> eigenvectors;

    // get the eigenvalues
    system_matrix_full.compute_eigenvalues_symmetric(-0.3, 0.3, 1e-6, eigenvalues, eigenvectors);

    bool eigenvector_flag = false;
    unsigned int numNegVals = 0;
    for (unsigned int i = 0 ; i < eigenvalues.size(); i ++)
    {
      double nextEigenVal = eigenvalues[i];
      if (nextEigenVal < 0.0)
      {
        numNegVals ++;

        // set the eigenvector corresponnding the negative eigenvalue specified by index
        if (numNegVals == index)
        {
          eigenvector_flag = true;

          initial_solution_tangent.reinit(dof_handler.n_dofs());
          for (unsigned int j = 0; j < dof_handler.n_dofs(); j++)
            initial_solution_tangent[j] = eigenvectors[j][i];

          break;
        }
      }
    }

    if (eigenvector_flag == false)
    {
      std::cout << "Unable to set unstable eigenvetor because there was not " << index << " negative eigenvalues. Exiting" << std::endl;
      exit(-1);
    }
    // So the eigenvector had bad entries for the constrained dofs, so need to distribute the constraints
    constraints.distribute(initial_solution_tangent);

    // also, normalize the tangent vector
    initial_solution_tangent *= (1.0/initial_solution_tangent.l2_norm());

    // and finally, just set the lambda tangent to 0.0.
    initial_lambda_tangent = 0.0;
  }

  template <int dim>
  double ElasticProblem<dim>::bisect_find_lambda_critical(double lowerBound, double upperBound,
                                                          double tol, unsigned int maxIter)
  {
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

  template <int dim>
  void ElasticProblem<dim>::output_results (const unsigned int cycle) const
  {

    std::vector<std::string> solution_names;
    switch (dim)
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

    DataOut<dim> data_out_totalDisp;

    data_out_totalDisp.attach_dof_handler (dof_handler);


    // Get the total displacement of each of the points.
    // Get the points of the dofs so we can do some shifting...
    std::vector<Point<dim>> support_points(dof_handler.n_dofs());
    MappingQ1<dim> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

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
    data_out_totalDisp.build_patches ();
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

    DataOut<dim> data_out_disp_from_uniform;

    data_out_disp_from_uniform.attach_dof_handler (dof_handler);

    data_out_disp_from_uniform.add_data_vector (present_solution, solution_names);
    data_out_disp_from_uniform.build_patches ();
    data_out_disp_from_uniform.write_vtk (output_disp_from_uniform);

    // Now output the deformed mesh

    // just need to shift the corrdinates of the verticies by the shifted solution vector

    DataOut<dim> deformed_data_out;

    deformed_data_out.attach_dof_handler(dof_handler);
    deformed_data_out.add_data_vector(shifted_solution, solution_names);

    MappingQEulerian<dim> q_mapping(1,  dof_handler, shifted_solution);
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

  template<int dim>
  void ElasticProblem<dim>::output_load_info(std::vector<double> lambda_values,
                                             std::vector<double> energy_values,
                                             std::vector<double> congugate_lambda_values,
                                             std::vector<double> displacement_magnitude,
                                             const unsigned int cycle) const
  {

    // output the lambda value, system energy, and the congugate lambda vales for each step

    std::string filename(output_directory);
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
      load_data_output << std::setprecision(15) << std::setw(25) << displacement_magnitude[i] << std::endl;
    }

    load_data_output.close();
  }

  template<int dim>
  void ElasticProblem<dim>::read_input_file(char* filename)
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

      // read in exponential growth parameter
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg", &kappa);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      mu = new MuFunction<dim>(kappa);

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

      // read in the critical frequency
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %u", &ds, &load_steps);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the critical frequency
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

    // make the output directory
    struct stat st;
    if (stat("./output", &st) == -1)
       mkdir("./output", 0700);

    if (stat(output_directory, &st) == -1)
      mkdir(output_directory, 0700);

  }

  template <int dim>
  void ElasticProblem<dim>::getNextDataLine( FILE* const filePtr, char* nextLinePtr,
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

  template <int dim>
  void ElasticProblem<dim>::save_current_state()
  {
    // create the output directory if it doesnt exist


    char saved_state_dir[MAXLINE];
    strcpy(saved_state_dir, output_directory);
    strcat(saved_state_dir, "/saved_state");

    // see if the directory exists
    struct stat st;
    if (stat(saved_state_dir, &st) == -1)
         mkdir(saved_state_dir, 0700);


    // write the mesh
    char mesh_file[MAXLINE];
    strcpy(mesh_file, saved_state_dir);
    strcat(mesh_file, "/mesh.msh");
    std::ofstream mesh_out (mesh_file);
    GridOut grid_out;
    grid_out.set_flags(GridOutFlags::Msh (true, true));
    grid_out.write_msh (triangulation, mesh_out);

    // Present Solution
    char present_solution_file[MAXLINE];
    strcpy(present_solution_file, saved_state_dir);
    strcat(present_solution_file, "/present_solution.dat");
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

  template <int dim>
  void ElasticProblem<dim>::load_state()
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

    // load the mesh
    char mesh_file[MAXLINE];
    strcpy(mesh_file, input_dir_path);
    strcat(mesh_file, "/mesh.msh");
    std::ifstream mesh_in (mesh_file);
    GridIn<2> grid_in;
    grid_in.attach_triangulation(triangulation);
    grid_in.read_msh(mesh_in);

    // Present Solution
    char present_solution_file[MAXLINE];
    strcpy(present_solution_file, input_dir_path);
    strcat(present_solution_file, "/present_solution.dat");
    std::ifstream solution_in (present_solution_file);
    boost::archive::text_iarchive solution_ar(solution_in);
    present_solution.load(solution_ar, 1);

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
    lambda_in >> present_lambda;
    lambda_in.close();




  }
}

#endif // COMPRESSEDSTRIPPACABLOCH_CC_
