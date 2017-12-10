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
#define NU_VALUE 0.3

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

    double nuValue = NU_VALUE;

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
    Assert (p.dimension == 2, ExcNotImplemented())

    return MU_VALUE;
  }

  void MuFunction::value_list(const std::vector< Point< DIM > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = MuFunction::value(points[i], component);
  }



  // computes right hand side values if we were to have body forces. But it just
  // always returns zeros because we don't.
  void ElasticProblem::right_hand_side (const std::vector<Point<DIM> > &points,
                        std::vector<Tensor<1, DIM> >   &values)
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));
    Assert (DIM >= 2, ExcNotImplemented());

    if (problem_2_flag == false)
    {
      for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
      {
//        std::cout << points[point_n][0]  << "  " << points[point_n][1] <<std::endl;
//        if (fabs(points[point_n][0] - domain_dimensions[0]) < 1e-6)
//        {
//          values[point_n][0] = -load_val;
//        }
//        else
          values[point_n][0] = 0.0;

        values[point_n][1] = 0.0;
      }
    }
    else
    {
      for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
      {
        values[point_n][0] = 0.0;
        values[point_n][1] = 0.0;
      }
    }
  }


  // Functions that get the tensors used in the stiffness and rhs calculations.

  inline
  Tensor<4, DIM>
  ElasticProblem::get_the_D(double mu, double nu)
  {
    Tensor<4, DIM> tmp;

    double scalingFactor = mu/(1.0 - nu*nu);
    double val = (1 - nu)/4.0;

    tmp[0][0][0][0] = 1.0;
    tmp[0][0][1][1] = nu;
    tmp[1][1][0][0] = nu;
    tmp[1][1][1][1] = 1.0;
    tmp[0][1][0][1] = val;
    tmp[0][1][1][0] = val;
    tmp[1][0][0][1] = val;
    tmp[1][0][1][0] = val;

    tmp *= scalingFactor;

    return tmp;
  }

  inline
  Tensor<2,DIM>
  ElasticProblem::get_deformation_gradient(std::vector<Tensor<1,DIM> > old_solution_gradient)
  {
    Tensor<2,DIM> tmp;

    for (unsigned int i = 0; i < DIM; i ++)
    {
      for(unsigned int j = 0; j < DIM; j++)
      {
        tmp[i][j] += old_solution_gradient[i][j];
      }
      tmp[i][i] += 1.0;
    }

    return tmp;
  }

  inline
  Tensor<2, DIM> ElasticProblem::get_lagrangian_strain(Tensor<2, DIM> F)
  {
    Tensor<2, DIM> tmp;

    for (unsigned int i = 0; i < DIM; i ++)
    {
      for(unsigned int j = 0; j < DIM; j++)
      {
        for(unsigned int k = 0; k <DIM; k++)
        {
          tmp[i][j] += F[k][i]*F[k][j];
        }
      }
      tmp[i][i] -= 1.0;
    }
    tmp *= 0.5;

    return tmp;
  }


  ElasticProblem::ElasticProblem ()
    :
    dof_handler (triangulation),
    fe (FE_Q<DIM>(1), DIM)
  {}




  ElasticProblem::~ElasticProblem ()
  {
    dof_handler.clear ();
  }


  void ElasticProblem::create_mesh()
  {

    // creates our strip.
    Point<DIM> corner1, corner2;
    corner1(0) =  0;
    corner1(1) =  -domain_dimensions[1]/2.0;
    corner2(0) =  domain_dimensions[0];
    corner2(1) =  domain_dimensions[1]/2.0;
    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);

    // Make sure to renumber the boundaries
    renumber_boundary_ids();

  }

  void ElasticProblem::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.


    dof_handler.distribute_dofs (fe);

    present_solution.reinit (dof_handler.n_dofs());

    setup_system_constraints();

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    drhs_dlambda.reinit (dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);

    if(problem_2_flag)
    {
      corner_dofs.resize(4, 0);
      inside_dofs.resize(4, 0);

      const unsigned int  number_dofs = dof_handler.n_dofs();

      std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
      MappingQ1<DIM> mapping;
      DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

      std::vector<bool> x2_components = {false, true};
      ComponentMask x2_mask(x2_components);

      std::vector<bool> is_x2_comp(number_dofs);

      DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);



      for(unsigned int i = 0; i < number_dofs; i++)
      {
        if( (fabs(support_points[i](0) - domain_dimensions[0]) < 1e-6) &&
            (fabs(fabs(support_points[i](1)) - domain_dimensions[1]/2.0) < 1e-6 ))
          {
            if(!is_x2_comp[i])
            {
              if(fabs(support_points[i](1) + domain_dimensions[1]/2.0) < 1e-6)
                corner_dofs[0] = i;
              else if((fabs(support_points[i](1) - domain_dimensions[1]/2.0) < 1e-6))
                corner_dofs[2] = i;
            }
            else
            {
              if(fabs(support_points[i](1) + domain_dimensions[1]/2.0) < 1e-6)
                corner_dofs[1] = i;
              else if((fabs(support_points[i](1) - domain_dimensions[1]/2.0) < 1e-6))
                corner_dofs[3] = i;
            }

          }

        if( (fabs(support_points[i](0) + domain_dimensions[0]/grid_dimensions[0] - domain_dimensions[0]) < 1e-6) &&
            (fabs(fabs(support_points[i](1)) - domain_dimensions[1]/2.0) < 1e-6 ))
          {
            if(!is_x2_comp[i])
            {
              if(fabs(support_points[i](1) + domain_dimensions[1]/2.0) < 1e-6)
                inside_dofs[0] = i;
              else if((fabs(support_points[i](1) - domain_dimensions[1]/2.0) < 1e-6))
                inside_dofs[2] = i;
            }
            else
            {
              if(fabs(support_points[i](1) + domain_dimensions[1]/2.0) < 1e-6)
                inside_dofs[1] = i;
              else if((fabs(support_points[i](1) - domain_dimensions[1]/2.0) < 1e-6))
                inside_dofs[3] = i;
            }

          }
      }

      for (unsigned int i = 0; i < corner_dofs.size(); i++)
      {
        for (unsigned int j = 0; j < corner_dofs.size(); j++)
          dsp.add(corner_dofs[i], corner_dofs[j]);
      }

    }

    sparsity_pattern.copy_from (dsp);


//    GridTools::distort_random(0.4, triangulation, true);

    system_matrix.reinit (sparsity_pattern);

    evaluation_point = present_solution;

  }

  void ElasticProblem::setup_system_constraints ()
  {

    constraints.clear ();

    // if we do problem 1, need to constrain the average x2 displacement to 0 on boundry
    if (problem_2_flag == false)
    {
      unsigned int number_dofs = dof_handler.n_dofs();

      // now do constraints that the average x1 displacement is zero
      std::vector<bool> x2_components = {false, true};
      ComponentMask x2_mask(x2_components);

      std::vector<bool> boundary_dof_x2 (number_dofs, false);
      std::vector<bool> boundary_12_dof_x2 (number_dofs, false);

      std::set< types::boundary_id > boundary_id_12;
      boundary_id_12.insert(1);
      boundary_id_12.insert(2);

      DoFTools::extract_boundary_dofs(dof_handler,
                                      x2_mask,
                                      boundary_dof_x2);

      DoFTools::extract_boundary_dofs (dof_handler,
                                       x2_mask,
                                       boundary_12_dof_x2,
                                       boundary_id_12);



      unsigned int first_boundary_dof = 0;
      for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
      {
        if ((boundary_dof_x2[i] == true) && (boundary_12_dof_x2[i] == false))
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

        if(boundary_dof_x2[i] == true)
          constraints.add_entry (first_boundary_dof, i, -1);

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

  void ElasticProblem::assemble_system_matrix()
  {
    // Assembling the system matrix. I chose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_matrix = 0.0;

    QGauss<1> quad_2(2);
    QGauss<1> quad_1(1);

    QAnisotropic<DIM> quadrature_problem_2(quad_1, quad_2);


    FEValues<DIM> fe_values (fe, quadrature_problem_2,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
 //   unsigned int   n_q_points    = quadrature_formula.size();
 //   if(problem_2_flag)
    unsigned int n_q_points = quadrature_problem_2.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,DIM>>(DIM));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    FullMatrix<double> E_dot_hat_n(DIM, DIM);
    FullMatrix<double> dE_dU_m(DIM, DIM);
    FullMatrix<double> dE_dot_hat_n_dU_m(DIM, DIM);

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_matrix = 0.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        D = get_the_D(mu_values[q_point], nu_values[q_point]);
        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> E = get_lagrangian_strain(F);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          E_dot_hat_n = 0.0;
          for(unsigned int i = 0; i<DIM; ++i)
            for(unsigned int j = 0; j<DIM; ++j)
            {
              E_dot_hat_n[i][j] = 0.5*(
                  fe_values.shape_grad(n, q_point)[i]*old_solution_gradients[q_point][component_n][j] +
                  fe_values.shape_grad(n, q_point)[j]*old_solution_gradients[q_point][component_n][i] +
                  (i == component_n ? 1.0 : 0.0)*fe_values.shape_grad(n, q_point)[j] +
                  (j == component_n ? 1.0: 0.0)*fe_values.shape_grad(n,q_point)[i]);
            }


          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {
            const unsigned int component_m = fe.system_to_component_index(m).first;

            dE_dU_m = 0.0;
            dE_dot_hat_n_dU_m = 0.0;
            for(unsigned int i = 0; i<DIM; ++i)
              for(unsigned int j = 0; j<DIM; ++j)
              {
                dE_dU_m[i][j] = 0.5*(
                    fe_values.shape_grad(m, q_point)[i]*old_solution_gradients[q_point][component_m][j] +
                    fe_values.shape_grad(m, q_point)[j]*old_solution_gradients[q_point][component_m][i] +
                    (i == component_m ? 1.0 : 0.0)*fe_values.shape_grad(m, q_point)[j] +
                    (j == component_m ? 1.0: 0.0)*fe_values.shape_grad(m,q_point)[i]);

                dE_dot_hat_n_dU_m[i][j] = 0.5*(
                    (component_m == component_n ? 1.0 : 0.0)*
                            fe_values.shape_grad(n, q_point)[i]*fe_values.shape_grad(m, q_point)[j] +
                    (component_n == component_m ? 1.0 : 0.0)*
                            fe_values.shape_grad(m, q_point)[j]*fe_values.shape_grad(n, q_point)[i]
                    );
              }
            for(unsigned int i = 0; i<DIM; ++i)
              for(unsigned int j = 0; j<DIM; ++j)
                for(unsigned int l = 0; l<DIM; ++l)
                  for(unsigned int r = 0; r<DIM; ++r)
                  {

                      cell_matrix[n][m] += D[i][j][l][r]*(
                          dE_dU_m[l][r]*E_dot_hat_n[i][j] +
                          E[l][r]*dE_dot_hat_n_dU_m[i][j] )*fe_values.JxW(q_point);
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

    apply_loads_to_system_matrix();
  }

  void ElasticProblem::assemble_system_rhs()
  {
    // Assembling the system rhs. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

    system_rhs = 0.0;

    QGauss<1> quad_2(2);
    QGauss<1> quad_1(1);

    QAnisotropic<DIM> quadrature_problem_2(quad_1, quad_2);


    FEValues<DIM> fe_values (fe, quadrature_problem_2,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    // unsigned int   n_q_points    = quadrature_formula.size();

      unsigned int n_q_points = quadrature_problem_2.size();

    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,DIM>>(DIM));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    FullMatrix<double> E_dot_hat_n(DIM, DIM);

    std::vector<Tensor<1, DIM> > rhs_values (n_q_points);

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list  (fe_values.get_quadrature_points(), mu_values);
      right_hand_side (fe_values.get_quadrature_points(), rhs_values);



      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        D = get_the_D(mu_values[q_point], nu_values[q_point]);
        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> E = get_lagrangian_strain(F);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          E_dot_hat_n = 0.0;
          for(unsigned int i = 0; i<DIM; ++i)
            for(unsigned int j = 0; j<DIM; ++j)
            {
              E_dot_hat_n[i][j] = 0.5*(
                  fe_values.shape_grad(n, q_point)[i]*old_solution_gradients[q_point][component_n][j] +
                  fe_values.shape_grad(n, q_point)[j]*old_solution_gradients[q_point][component_n][i] +
                  (i == component_n ? 1.0 : 0.0)*fe_values.shape_grad(n, q_point)[j] +
                  (j == component_n ? 1.0: 0.0)*fe_values.shape_grad(n,q_point)[i]);
            }

          for(unsigned int i = 0; i<DIM; ++i)
            for(unsigned int j = 0; j<DIM; ++j)
              for(unsigned int l = 0; l<DIM; ++l)
                for(unsigned int m = 0; m<DIM; ++m)
                {
                  cell_rhs[n] += D[i][j][l][m]*E[l][m]*E_dot_hat_n[i][j]*fe_values.JxW(q_point);
                }

       //   cell_rhs[n] += fe_values.shape_value(n, q_point)*rhs_values[q_point][component_n]*fe_values.JxW(q_point);
        }
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        system_rhs(local_dof_indices[n]) += cell_rhs(n);

    }

    apply_loads_to_rhs();

    constraints.condense (system_rhs);

    system_rhs *= -1.0;

  }

  void ElasticProblem::apply_loads_to_rhs()
  {
    unsigned int number_dofs = dof_handler.n_dofs();

    if(problem_2_flag)
    {
      double u11, u12, u21, u22;
      u11 = evaluation_point[corner_dofs[0]];
      u12 = evaluation_point[corner_dofs[1]];
      u21 = evaluation_point[corner_dofs[2]];
      u22 = evaluation_point[corner_dofs[3]];

      double h = domain_dimensions[1];

      double a0 = load_val/((u21 - u11)*(u21 - u11) + (u22 + h - u12)*(u22 + h - u12));

      system_rhs[corner_dofs[0]] -= a0*(h + u22 - u12);
      system_rhs[corner_dofs[1]] -= a0*(u11 - u21);

      system_rhs[corner_dofs[2]] -= -a0*(h + u22 - u12);
      system_rhs[corner_dofs[3]] -= -a0*(u11 - u21);

    }
    else
    {
      std::vector<bool> boundary_2_dofs_x1 (dof_handler.n_dofs(), false);

      std::vector<bool> boundary_34_dofs_x1(dof_handler.n_dofs(), false);

      std::set< types::boundary_id > boundary_id_2;
      boundary_id_2.insert(2);

      std::set< types::boundary_id > boundary_id_34;
      boundary_id_34.insert(3);
      boundary_id_34.insert(4);

      std::vector<bool> x1_components = {true, false};

      ComponentMask comp_mask(x1_components);

      DoFTools::extract_boundary_dofs(dof_handler,
                                       comp_mask,
                                       boundary_2_dofs_x1,
                                       boundary_id_2);

      DoFTools::extract_boundary_dofs(dof_handler,
                                      comp_mask,
                                      boundary_34_dofs_x1,
                                      boundary_id_34);

      for (unsigned int i = 0; i<number_dofs; i ++)
      {
        if(boundary_2_dofs_x1[i])
        {
          system_rhs[i] -= load_val/(grid_dimensions[1]);

          if(boundary_34_dofs_x1[i])
          {
            system_rhs[i] += load_val/(2*grid_dimensions[1]);
          }
        }
      }

    }
  }

  void ElasticProblem::apply_loads_to_system_matrix()
  {
    if(problem_2_flag)
    {
      double u11, u12, u21, u22;
      unsigned int x11, x12, x21, x22;
      x11 = corner_dofs[0];
      x12 = corner_dofs[1];
      x21 = corner_dofs[2];
      x22 = corner_dofs[3];

      u11 = evaluation_point[x11];
      u12 = evaluation_point[x12];
      u21 = evaluation_point[x21];
      u22 = evaluation_point[x22];
      double h = domain_dimensions[1];

//      std::cout << u11 << " " << u12 << " " << u21 << " " << u22 << std::endl;
//      std::cout << ((u21 - u11)*(u21 - u11) + (u22 + h - u12)*(u22 + h - u12)) << std::endl;



      double a0 = 1.0/((u21 - u11)*(u21 - u11) + (u22 + h - u12)*(u22 + h - u12));
      double a1 = a0*a0;
      a0 *= load_val;
      a1 *= -load_val;

      double dF11_du11 = -a1*2*(u21 - u11)*(u22 + h - u12);
      double dF11_du12 = -a1*2*(u22 + h - u12)*(u22 + h - u12) - a0;
      double dF11_du21 = a1*2*(u21 - u11)*(u22 + h - u12);
      double dF11_du22 =  a1*2*(u22 + h - u12)*(u22 + h - u12) + a0;


      double dF12_du11 = -a1*2*(u21 - u11)*(u11 - u21) + a0;
      double dF12_du12 = -a1*2*(u22 + h - u12)*(u11 - u21);
      double dF12_du21 = a1*2*(u21 - u11)*(u11 - u21) - a0;
      double dF12_du22 = a1*2*(u22 + h - u12)*(u11 - u21);

//      system_matrix.add(x11, x11, -dF11_du11);
//      system_matrix.add(x12, x11, -dF11_du12);
//      system_matrix.add(x21, x11, -dF11_du21);
//      system_matrix.add(x22, x11, -dF11_du22);
//
//      system_matrix.add(x11, x12, -dF12_du11);
//      system_matrix.add(x12, x12, -dF12_du12);
//      system_matrix.add(x21, x12, -dF12_du21);
//      system_matrix.add(x22, x12, -dF12_du22);
//
//      system_matrix.add(x11, x21, dF11_du11);
//      system_matrix.add(x12, x21, dF11_du12);
//      system_matrix.add(x21, x21, dF11_du21);
//      system_matrix.add(x22, x21, dF11_du22);
//
//      system_matrix.add(x11, x22, dF12_du11);
//      system_matrix.add(x12, x22, dF12_du12);
//      system_matrix.add(x21, x22, dF12_du21);
//      system_matrix.add(x22, x22, dF12_du22);

      system_matrix.add(x11, x11, -dF11_du11);
      system_matrix.add(x11, x12, -dF11_du12);
      system_matrix.add(x11, x21, -dF11_du21);
      system_matrix.add(x11, x22, -dF11_du22);

      system_matrix.add(x12, x11, -dF12_du11);
      system_matrix.add(x12, x12, -dF12_du12);
      system_matrix.add(x12, x21, -dF12_du21);
      system_matrix.add(x12, x22, -dF12_du22);

      system_matrix.add(x21, x11, dF11_du11);
      system_matrix.add(x21, x12, dF11_du12);
      system_matrix.add(x21, x21, dF11_du21);
      system_matrix.add(x21, x22, dF11_du22);

      system_matrix.add(x22, x11, dF12_du11);
      system_matrix.add(x22, x12, dF12_du12);
      system_matrix.add(x22, x21, dF12_du21);
      system_matrix.add(x22, x22, dF12_du22);


    }
  }

  void ElasticProblem::apply_boundaries_and_constraints_system_matrix()
  {
    constraints.condense (system_matrix);

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side1_components = {true, false};

    if(problem_2_flag)
      side1_components = {true, true};


    ComponentMask side1_mask(side1_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              ZeroFunction<DIM, double>(DIM),
                                              boundary_values,
                                              side1_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);

  }

  void ElasticProblem::apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs)
  {
    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      if ((*homogenous_dirichlet_dofs)[i] == true)
        (*rhs)[i] = 0.0;
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
    std::vector<bool> homogenous_boundary_dofs (dof_handler.n_dofs(), false);

    std::set< types::boundary_id > boundary_id_1;
    boundary_id_1.insert(1);


    std::vector<bool> homo_components = {true, false};
    if (problem_2_flag)
      homo_components = {true, true};

    ComponentMask comp_mask(homo_components);

    DoFTools::extract_boundary_dofs(dof_handler,
                                       comp_mask,
                                       homogenous_boundary_dofs,
                                       boundary_id_1);

    // Starts by getting the residual in the current configuration.
    evaluation_point = present_solution;

    assemble_system_rhs();
    apply_boundaries_to_rhs(&system_rhs, &homogenous_boundary_dofs);

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
      line_search_and_add_step_length(current_residual, &homogenous_boundary_dofs);

      evaluation_point = present_solution;
      current_residual = system_rhs.l2_norm();
      iteration ++;
    }

    // output iterations for convergance.
    std::cout << "    Converging Iterations : " << iteration << "\n";
    std::cout << "    Residual              : " << current_residual << "\n";

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



  unsigned int ElasticProblem::get_system_eigenvalues(const int cycle)
  {
    // get the system's eigenvalues. I really don't need to have it take in lambda_eval
    // and could just have it use the present_lambda, but its fine. the cycle is the
    // the output file will have appended on its name. -1 for no output. Outputs the number
    // of negative eigenvalues in the range -10 < lambda < 0.3

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


  void ElasticProblem::output_results (const unsigned int cycle) const
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

    DataOut<DIM> data_out_totalDisp;

    data_out_totalDisp.attach_dof_handler (dof_handler);


    // Get the total displacement of each of the points.
    // Get the points of the dofs so we can do some shifting...
    std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
    MappingQ1<DIM> mapping;
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
        shifted_solution[i] = present_solution[i];
      }
      else
      {
        // it is an x2 component
        shifted_solution[i] = present_solution[i];
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

    DataOut<DIM> data_out_disp_from_uniform;

    data_out_disp_from_uniform.attach_dof_handler (dof_handler);

    data_out_disp_from_uniform.add_data_vector (present_solution, solution_names);
    data_out_disp_from_uniform.build_patches ();
    data_out_disp_from_uniform.write_vtk (output_disp_from_uniform);

    // Now output the deformed mesh

    // just need to shift the corrdinates of the verticies by the shifted solution vector

    DataOut<DIM> deformed_data_out;

    deformed_data_out.attach_dof_handler(dof_handler);
    deformed_data_out.add_data_vector(shifted_solution, solution_names);

    MappingQEulerian<DIM> q_mapping(1,  dof_handler, shifted_solution);
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
                                             const unsigned int cycle) const
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
      load_data_output << std::setprecision(15) << std::setw(25) << displacement_magnitude[i] << std::endl;
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
      // Read in which problem we are doing
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      unsigned int problem_number = 0;
      valuesWritten = sscanf(nextLine, "%u", &problem_number);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      if(problem_number == 2)
        problem_2_flag = true;
      else
        problem_2_flag = false;

      // Read in the output name
      char directory_name[MAXLINE];
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%s", directory_name);
      if (valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

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

      // Read in the domain dimensions
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg", &domain_dimensions[0], &domain_dimensions[1]);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the absolute tolerance of newton iteration and max iterations
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg  %u", &tol, &maxIter);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
      }



      // read in the final load steps
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %u", &final_load, &load_steps);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the  output every
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
    dof_handler.distribute_dofs(fe);
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
        if (fabs(face_center[0]) < 1e-4)  // left
        {
          cell->face(f)->set_boundary_id (1);
        }
        else if (fabs(face_center[0] - domain_dimensions[0]) < 1e-4)
        {// right
          cell->face(f)->set_boundary_id (2);
        }
        else if (fabs(face_center[1] + domain_dimensions[1]/2.0) < 1e-6) //bottom
        {
          cell->face(f)->set_boundary_id (3);
        }
        else if (fabs(face_center[1] - domain_dimensions[1]/2.0) < 1e-6) //top
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
    MappingQ1<DIM> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

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

}

#endif // COMPRESSEDSTRIPPACABLOCH_CC_
