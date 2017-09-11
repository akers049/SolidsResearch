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

#include "NeoHookean_Newton_CompressedStrip.h"

#include <fstream>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>


#define MAXLINE 1024
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
    domain_dimensions[0] = (2.0*(4.0*atan(1.0))/critical_frequency)/10.0;
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
    }
    */
  }

  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.

    dof_handler.distribute_dofs (fe);
    present_solution.reinit (dof_handler.n_dofs());
    evaluation_point = present_solution;

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    setup_constraints();

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

  template <int dim>
  void ElasticProblem<dim>::setup_constraints ()
  {
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             hanging_node_constraints);

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

    std::vector<bool> x1_x2_components = {true, true};
    ComponentMask x1_x2_mask(x1_x2_components);


/*
    unsigned int first_boundary_dof = 0;
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if ((boundary_dof_x1[i] == true) && (boundary_01_dof_x1[i] == false))
      {
        first_boundary_dof = i;
        break;
      }
    }
    hanging_node_constraints.add_line (first_boundary_dof);
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if (i == first_boundary_dof)
        continue;

      if(boundary_dof_x1[i] == true)
        hanging_node_constraints.add_entry (first_boundary_dof, i, -1);

    }
*/
    /*
   std::vector<bool> x1_x2_components = {true, true};
    ComponentMask x1_x2_mask(x1_x2_components);
    DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0,  constraints, x1_x2_mask);*/
    //DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 1,  constraints, x1_x2_mask);

    /*
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
                                  (fabs(x2_coord - support_points[j](1)) < 1e-12))
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
                                  (fabs(x2_coord - support_points[j](1)) < 1e-12))
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
    }*/
    hanging_node_constraints.close();

DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0,  constraints, x1_x2_mask);
constraints.close();

    update_bloch_wave_constraints(1); // just use 1 for now to have correct sparsity stuff.

  }

  template<int dim>
  void ElasticProblem<dim>::update_bloch_wave_constraints(unsigned int periodicity_number)
  {
    constraints_bloch.clear();

    double k_x1 = 2.0*(4.0*atan(1.0))/periodicity_number;

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
  void ElasticProblem<dim>::assemble_system_matrix()
  {
    // Assembling the system matrix. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_matrix = 0.0;
    bloch_matrix = 0.0;

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
      {
        for (unsigned int m=0; m<dofs_per_cell; ++m)
        {
          system_matrix.add (local_dof_indices[n],
                             local_dof_indices[m],
                             cell_matrix(n,m));
        }
      }
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
  void ElasticProblem<dim>::apply_boundaries_and_constraints_bloch_matrix(SparseMatrix<double> *blochMat)
  {
    constraints_bloch.condense(*blochMat);
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

    unsigned int m = blochMat->m();
    // set values on the diagonal to the first diagonal element,
    // or 1 if it is nonexistent
    double first_nonzero_diagonal_entry = 1.0;
    for (unsigned int i=0; i<m; ++i)
    {
      if (blochMat->diag_element(i) != 0.0)
        {
          first_nonzero_diagonal_entry = blochMat->diag_element(i);
          break;
        }
    }
    // now march through matrix, zeroing out rows and columns.
    // If there is a current value on the diagonal of the constrained
    // boundary dof, don't touch it. If there is not one, then we can
    // just set it equal to the first nonzero entry we just found
    for (unsigned int row = 0; row < m; ++row)
    {

      const typename SparseMatrix<double>::iterator end_row = blochMat->end(row);
      for (typename SparseMatrix<double>::iterator entry = blochMat->begin(row);
                    entry != end_row; ++entry)
      {
        if((boundary_2_dof_x2[row%number_dofs] || boundary_2_dof_x2[entry->column()%number_dofs])
            && (row != entry->column()))
        {
          entry->value() = 0.0;
        }
        else if(boundary_2_dof_x2[row%number_dofs] && (row == entry->column()) &&
                (blochMat->diag_element(row%number_dofs) == 0.0))
        {
          entry->value() = first_nonzero_diagonal_entry;
        }

      }

    }

   }

  template<int dim>
  void ElasticProblem<dim>::assemble_bloch_matrix()
  {
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


  template<int dim>
  unsigned int ElasticProblem<dim>::get_system_eigenvalues(double lambda_eval, const int cycle,
                                                     const unsigned int periodicity_number)
  {

    update_F0(lambda_eval);
    assemble_system_matrix();
    assemble_bloch_matrix();

    apply_boundaries_and_constraints_bloch_matrix(&bloch_matrix);

    // copy current sparse system matrix to a full matrix.

    LAPACKFullMatrix<double> bloch_matrix_full;
    bloch_matrix_full.copy_from(bloch_matrix);


    Vector<double> eigenvalues;
    FullMatrix<double> eigenvectors;

    bloch_matrix_full.compute_eigenvalues_symmetric(-1, 0.3, 1e-6, eigenvalues, eigenvectors);
    //bloch_matrix_full.compute_eigenvalues();

    unsigned int number_negative_eigs = 0;
    if (cycle != -1)
    {
      struct stat st;
      if (stat("./output", &st) == -1)
          mkdir("./output", 0700);

      std::ostringstream cycle_str;
      cycle_str << cycle;

      std::string filename = "output/eigenvalues-";
      filename += cycle_str.str();

      std::ofstream outputFile;
      outputFile.open(filename.c_str());

      outputFile << "# eigenvalues of the system matrix" << std::endl;
      outputFile << periodicity_number << std::endl;
      outputFile << std::setprecision(10) << lambda_eval << std::endl;
/*
      for (unsigned int i = 0; i < 2*dof_handler.n_dofs(); i ++)
            {
              for(unsigned int j = 0 ; j < 2*dof_handler.n_dofs(); j ++)
              {

                if (bloch_matrix.el(i,j) > 1e-13)
                  outputFile << std::setprecision(3) << std::setw(5) << (bloch_matrix.el(i,j)) << " ";
                else
                  outputFile<< std::setprecision(3) << std::setw(5) << 0 << " ";
              }
            outputFile << "\n";
            }

      outputFile << "\n\n\n\n";
      for (unsigned int i = 0; i < dof_handler.n_dofs(); i ++)
            {
              for(unsigned int j = 0 ; j < dof_handler.n_dofs(); j ++)
              {

                if (system_matrix.el(i,j) > 1e-13)
                  outputFile << std::setprecision(3) << std::setw(5) << (system_matrix.el(i,j)) << " ";
                else
                  outputFile<< std::setprecision(3) << std::setw(5) << 0 << " ";
              }
            outputFile << "\n";
            }
      outputFile.close();
      exit(-1); */

      outputFile << "\n";
      for (unsigned int i = 0 ; i < eigenvalues.size(); i ++) //eigenvalues.size(); i ++)
      {

        outputFile << std::setprecision(15) << eigenvalues[i] << std::endl;
//          outputFile << std::setprecision(15) << bloch_matrix_full.eigenvalue(i) << std::endl;

      }



     /* for (unsigned int i = 0; i < 2*dof_handler.n_dofs(); i ++)
      {
        for(unsigned int j = 0 ; j < 2*dof_handler.n_dofs(); j ++)
          outputFile << std::setprecision(3) << std::setw(3) << (bloch_matrix.el(i,j)) << " ";
      outputFile << "\n";
      }*/

      outputFile.close();
    }
    else
    {
      for (unsigned int i = 0 ; i < eigenvalues.size(); i ++) //eigenvalues.size(); i ++)
      {
        if (eigenvalues[i] < 0.0)
          number_negative_eigs ++;
      }
    }

    return number_negative_eigs;
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
      // Read in the grid dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u", &grid_dimensions[0], &grid_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

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



    FILE* fid_1;
    fid_1 = std::fopen("lambda_guesses.dat", "r");
    double next_lambda;
    unsigned int next_periodicity;
    getNextDataLine(fid_1, nextLine, MAXLINE, &endOfFileFlag);
    while (1 != endOfFileFlag)
    {
      valuesWritten = sscanf(nextLine, "%u %lg", &next_periodicity, &next_lambda);
      if(valuesWritten != 2)
      {
        std::cout << "Error reading the lambda guess stuff. Exiting \n";
        fclose(fid_1);
        exit(-1);
      }
      lambda_guesses.push_back(next_lambda);
      periodicity_vals.push_back(next_periodicity);

      getNextDataLine(fid_1, nextLine, MAXLINE, &endOfFileFlag);

    }
    fclose(fid_1);
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
  void ElasticProblem<dim>::run ()
  {

    char fileName[MAXLINE];
    std::cout << "Please enter an input file: " << std::endl;
    std::cin >> fileName;
    read_input_file(fileName);

    create_mesh();

    setup_system();

    std::cout << "\n   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;


    std::cout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl << std::endl;


    std::vector<double> lambda_crits;
    for (unsigned int i = 0; i < periodicity_vals.size(); i++)
    {
      update_bloch_wave_constraints(periodicity_vals[i]);
      double lambda_crit = bisect_find_lambda_critical(lambda_guesses[i] - 0.05,
                                                       lambda_guesses[i] + 0.05, 1e-5, 50);

      std::cout << periodicity_vals[i] << ":  " << lambda_crit << "\n";
      lambda_crits.push_back(lambda_crit);

    }
    std::string filename = "output/lambda_crits.dat";

    std::ofstream outputFile;
    outputFile.open(filename.c_str());

    outputFile << "# peridoocity number    lambda_crit" << std::endl;
    for (unsigned int i = 0 ; i < periodicity_vals.size(); i ++)
      outputFile << periodicity_vals[i] << "     " << lambda_crits[i] << std::endl;

    outputFile.close();

  }
}





int main ()
{
  try
    {
      NeoHookean_Newton::ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}


