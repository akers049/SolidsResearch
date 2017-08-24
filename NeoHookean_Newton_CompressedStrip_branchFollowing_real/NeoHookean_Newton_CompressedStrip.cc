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
  }


  template <int dim>
  void ElasticProblem<dim>::create_mesh()
  {

    // set domain dimensions
    domain_dimensions[0] = 2.0*(4.0*atan(1.0))/critical_frequency;
    domain_dimensions[1] = 1.0;

    // creates our strip.
    Point<dim> corner1, corner2;
    corner1(0) = -domain_dimensions[0]/2.0;
    corner1(1) =  0.0;
    corner2(0) =  domain_dimensions[0]/2.0;
    corner2(1) =  domain_dimensions[1];
    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);
  }

  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.

    dof_handler.distribute_dofs (fe);
    present_solution.reinit (dof_handler.n_dofs());

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    setup_constraints();

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

  }

  template <int dim>
  void ElasticProblem<dim>::setup_constraints ()
  {
    constraints.clear ();

    DoFTools::make_hanging_node_constraints (dof_handler, constraints);


    // periodic in x1 and x2 for the 0 and 1 face (x1 faces)
    DoFTools::make_periodicity_constraints<DoFHandler<dim>>(dof_handler, 0, 1, 0, constraints);

    // now do constraints that the average x1 displacement on boundary 2 is zero

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

    // and the vector if it is a boundary 0 or 1 x2 dof. These are already constrained to each other.
    std::vector<bool> boundary_01_dof_x2 (number_dofs, false);
    DoFTools::extract_boundary_dofs (dof_handler,
                                     x2_mask,
                                     boundary_01_dof_x2,
                                     boundary_id_01);


    matched_dofs.resize(number_dofs, -1);
    for(unsigned int i = 0; i < number_dofs; i++)
    {
      if (is_x2_comp[i])
      {
        if (boundary_2_dof_x2[i] || boundary_01_dof_x2[i] || (matched_dofs[i] != -1))
        {
          // these dofs have already been constrained!
          continue;
        }
        else if (support_points[i](0) == 0.0)
        {
          // these are in the "miidle" so we don't constrained them
          continue;
        }
        else
        {
          double x1_coord = support_points[i](0);
          double x2_coord = support_points[i](1);
          for (unsigned int j = 0; j < number_dofs; j++)
          {

            if (i == j)
              continue;

            if (is_x2_comp[j])
            {
              if ((fabs(x1_coord + support_points[j](0)) < 1e-12 ) && (fabs(x2_coord - support_points[j](1)) < 1e-12))
              {
                constraints.add_line (i);
                constraints.add_entry (i, j, 1);
                matched_dofs[i] = j;
                matched_dofs[j] = i;
                continue;
              }

            }
            else
              continue;
            }
          }
        }
      }


    constraints.close ();

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
  void ElasticProblem<dim>::add_first_bif_displacements(double epsilon)
  {
   /* add the bifurcated solution to the current (zero) solution vector.
    * epsilon scales the magnitude of the added bifurcated solution.
    */

    std::vector<Point<dim>> support_points(dof_handler.n_dofs());
    MappingQ1<dim> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

    const unsigned int   number_dofs = dof_handler.n_dofs();

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> is_x1_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    for (unsigned int i = 0; i < number_dofs; i++)
    {
      if(is_x1_comp[i])
      {
        // it is an x1 component

        double v1 = 0.0;
        for (int j = 0; j < 4; j++)
          v1 += amplitudes_v1[j]*exp(charateristic_roots[j]*support_points[i](1));

        present_solution[i] += -epsilon*sin(critical_frequency*support_points[i](0))*v1;
      }
      else
      {
        // it is an x2 component

        double v2 = 0.0;
        for (int j = 0; j < 4; j++)
          v2 += amplitudes_v2[j]*exp(charateristic_roots[j]*support_points[i](1));

        present_solution[i] += epsilon*cos(critical_frequency*support_points[i](0))*v2;

      }
    }

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

    MuFunction<dim>  mu(kappa);
    NuFunction<dim>  nu;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_matrix = 0.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);


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

    constraints.condense (system_matrix);

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side2_components = {false, true};
    ComponentMask side2_mask(side2_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
                                              ZeroFunction<dim>(dim),
                                              boundary_values,
                                              side2_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);

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

    MuFunction<dim> mu(kappa);
    NuFunction<dim> nu;

    std::vector<Tensor<1, dim> > rhs_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);
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

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side2_components = {false, true};
    ComponentMask side2_mask(side2_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
                                              ZeroFunction<dim>(dim),
                                              boundary_values,
                                              side2_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);
  }

  template <int dim>
  void ElasticProblem<dim>::apply_boundaries_to_rhs(Vector<double> rhs, std::vector<bool> homogenous_dirichlet_dofs)
  {
    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      if (homogenous_dirichlet_dofs[i] == true)
        rhs[i] = 0.0;
    }
  }

  template <int dim>
  void ElasticProblem<dim>::assemble_system_energy_and_congugate_lambda(double lambda)
  {

    system_energy = 0.0;
    congugate_lambda = 0.0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points,
                                                std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    NuFunction<dim> nu;
    MuFunction<dim> mu(kappa);

    std::vector<Tensor<1, dim> > rhs_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);


      double dlambda1_dlambda, dlambda2_dlambda, a_nu;
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,dim> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, dim> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<2,dim> dW_dF = get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point],
                                                        F, F_inv, II_F);

        dlambda1_dlambda = -1.0;

        a_nu = (2.0*nu_values[q_point]/(1.0 - nu_values[q_point]));
        dlambda2_dlambda = ( a_nu + sqrt(4 + a_nu + a_nu*a_nu))
                            / ((2 + 2*a_nu)*(1.0 - lambda)*(1.0 - lambda));

        congugate_lambda += (dlambda1_dlambda*dW_dF[0][0] + dlambda2_dlambda*dW_dF[1][1])
                             *fe_values.JxW(q_point);

        double W = get_energy(nu_values[q_point], mu_values[q_point], F, II_F);
        system_energy += W*fe_values.JxW(q_point);
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
    apply_boundaries_to_rhs(system_rhs, boundary_2_dof_x2);

    current_residual = system_rhs.l2_norm();

    // Loops until coverge or go over max iterations
    while((current_residual > tol) &&
             (iteration < maxIter))
    {
      // Assemble the stiffness matrix
      assemble_system_matrix();

      // solve for the newton step
      solve();

      // Find the step length and add it to the current solution.
      // This function also calls assemble_system_rhs() so we don't need to
      // do another rhs call.
      line_search_and_add_step_length(current_residual, boundary_2_dof_x2);

      evaluation_point = present_solution;
      current_residual = system_rhs.l2_norm();
      iteration ++;
    }

    // output iterations for convergance.
    std::cout << "    Converging Iterations : " << iteration << "\n";

  }

  template<int dim>
  void ElasticProblem<dim>::line_search_and_add_step_length(double last_residual, std::vector<bool> homogenous_dirichlet_dofs)
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
      apply_boundaries_to_rhs(system_rhs, homogenous_dirichlet_dofs);

      current_residual = system_rhs.l2_norm();

      if(current_residual < last_residual)
        break;

    }
    present_solution = evaluation_point;
  }

  template <int dim>
  void ElasticProblem<dim>::solve ()
  {

    /*
    SolverControl           solver_control (10000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, newton_update, system_rhs,
              preconditioner);
*/
    SparseDirectUMFPACK  A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult (newton_update, system_rhs);
    constraints.distribute (newton_update);

  }

  template<int dim>
  bool ElasticProblem<dim>::get_system_eigenvalues(double lambda_eval, const int cycle)
  {


    update_F0(lambda_eval);
    assemble_system_matrix();

    // copy current sparse system matrix to a full matrix.
    LAPACKFullMatrix<double> system_matrix_full;
    system_matrix_full.copy_from(system_matrix);

    Vector<double> eigenvalues;
    FullMatrix<double> eigenvectors;

    system_matrix_full.compute_eigenvalues_symmetric(-1.0, 1.0, 1e-6, eigenvalues, eigenvectors);


    bool positive_definite = true;
    if (cycle != -1)
    {
      std::ostringstream cycle_str;
      cycle_str << cycle;

      std::string filename = "output/eigenvalues-";
      filename += cycle_str.str();

      std::ofstream outputFile;
      outputFile.open(filename.c_str());

      outputFile << "# eigenvalues of the system matrix" << std::endl;

      for (unsigned int i = 0 ; i < eigenvalues.size(); i ++)
      {
        double nextEigenVal = eigenvalues[i];

        outputFile << std::setprecision(15) << nextEigenVal << std::endl;

        if (nextEigenVal < 0.0)
          positive_definite = false;

      }

      outputFile << "\nIs positive definite : " << positive_definite << std::endl;
      outputFile.close();
    }
    else
    {
      for (unsigned int i = 0; i < eigenvalues.size(); i++)
      {
        if (eigenvalues[i] < 0.0)
        {
          positive_definite = false;
          break;
        }
      }
    }

    return positive_definite;
  }

  template <int dim>
  double ElasticProblem<dim>::bisect_find_lambda_critical(double lowerBound, double upperBound, double tol, unsigned int maxIter)
  {
    unsigned int N = 1;
    double middleVal = 0.0;
    while ( N < maxIter)
    {
      middleVal = (upperBound + lowerBound)/2.0;

      if ((upperBound - lowerBound)/2 < tol)
        return middleVal;

      N += 1;

      if (get_system_eigenvalues(middleVal, -1) && get_system_eigenvalues(lowerBound, -1))
        lowerBound = middleVal;
      else
        upperBound = middleVal;
    }

    return middleVal;
  }

  template <int dim>
  void ElasticProblem<dim>::output_results (const unsigned int cycle) const
  {

    std::ostringstream cycle_str;
    cycle_str << cycle;

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

    // See if the output file exists. If not, make the directory.
    struct stat st;
    if (stat("./output", &st) == -1)
        mkdir("./output", 0700);

    // output the total displacements. this requires adding in the uniform solution on top of the displacements
    std::string filename0 = "output/total_displacement-";
    filename0 += cycle_str.str();

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

    std::string filename1 = "output/displacement_from_uniform-";
    filename1 += cycle_str.str();

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

    std::string filename2 = "output/deformed_mesh-";
    filename2 += cycle_str.str();
    filename2 += ".vtk";
    std::ofstream output_deformed_mesh(filename2.c_str());
    deformed_data_out.write_vtk(output_deformed_mesh);


  }

  template<int dim>
  void ElasticProblem<dim>::output_load_info(std::vector<double> lambda_values,
                                             std::vector<double> energy_values,
                                             std::vector<double> congugate_lambda_values) const
  {

    // see if output directory exists, if not create it
    struct stat st;
    if (stat("./output", &st) == -1)
        mkdir("./output", 0700);

    // output the lambda value, system energy, and the congugate lambda vales for each step
    std::string filename = "output/load_info.txt";
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
      load_data_output << std::setprecision(15) << std::setw(25) << congugate_lambda_values[i] << std::endl;
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
    charateristic_roots.resize(4);
    amplitudes_v1.resize(4);
    amplitudes_v2.resize(4);

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

      // read in the roots to the characteristic
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg",
          &charateristic_roots[0], &charateristic_roots[1], &charateristic_roots[2], &charateristic_roots[3]);
      if(valuesWritten != 4)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the v1 amplitudes
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg",
          &amplitudes_v1[0], &amplitudes_v1[1], &amplitudes_v1[2], &amplitudes_v1[3]);
      if(valuesWritten != 4)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the v2 amplitudes
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg",
          &amplitudes_v2[0], &amplitudes_v2[1], &amplitudes_v2[2], &amplitudes_v2[3]);
      if(valuesWritten != 4)
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


    // get the critical lambda value
    evaluation_point = present_solution;
    double lambda_c = bisect_find_lambda_critical(critical_lambda_analytical - 0.01, critical_lambda_analytical + 0.01, 1e-6, 50);

    std::cout << "The lambda_c is: " << lambda_c << std::endl;

    update_F0(0);
    newton_iterate();
    output_results (0);

    double lambda_start = lambda_c - 1e-5;
    update_F0(lambda_start);
    newton_iterate();
    output_results (1);

    add_first_bif_displacements(-0.01);
    output_results(999);
    newton_iterate();
    output_results (2);


    unsigned int load_steps = 3000;
    std::vector<double> lambda_values(load_steps);
    std::vector<double> congugate_lambda_values(load_steps);
    std::vector<double> energy_values(load_steps);
   for(unsigned int i = 1; i < load_steps; i ++)
    {

     double lambda = lambda_start - 5e-7*i;
      update_F0(lambda);
      newton_iterate();


    /*  // get energy and congugate lambda value and save them.
            assemble_system_energy_and_congugate_lambda(lambda);
            lambda_values[i] = lambda;
            congugate_lambda_values[i] = congugate_lambda;
            energy_values[i] = system_energy;*/

      if ((i)%50 == 0)
      {
       output_results((i)/50 + 2);
       get_system_eigenvalues(lambda, (i)/50);
      }
    }
   // output_load_info(lambda_values, energy_values, congugate_lambda_values);

    /*
    update_F0(critical_lambda + 0.000002);
        newton_iterate(tol, 50);
        output_results (3);
*/


    // output initial mesh
/*
    std::vector<double> lambda_values(load_steps);
    std::vector<double> congugate_lambda_values(load_steps);
    std::vector<double> energy_values(load_steps);

    double lambda_step = final_lambda/load_steps;
    double lambda = 0.0;
    for(unsigned int i = 0; i < load_steps; i++)
    {

      // update lambda
      lambda +=lambda_step;
      std::cout << "Load step "<< i + 1 << " With loading parameter lambda = " << lambda << std::endl;

      update_F0(lambda);

      newton_iterate(tol, 50);

      // get energy and congugate lambda value and save them.
      assemble_system_energy_and_congugate_lambda(lambda);
      lambda_values[i] = lambda;
      congugate_lambda_values[i] = congugate_lambda;
      energy_values[i] = system_energy;

      std::cout << "    System Energy: " << system_energy << "\n\n";

      // output data if we're on the right step.
      if((i+1)%output_every_n == 0)
        output_results ((i+1)/output_every_n);
    }
*/
    // output_load_info(lambda_values, energy_values, congugate_lambda_values);
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


