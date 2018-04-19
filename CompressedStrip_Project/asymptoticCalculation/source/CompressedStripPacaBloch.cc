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

#define DIM 2



namespace compressed_strip
{
  using namespace dealii;


  template< typename F >
    class gsl_function_pp : public gsl_function {
    public:
    gsl_function_pp(const F& func) : _func(func) {
      function = &gsl_function_pp::invoke;
      params=this;
    }
    private:
    const F& _func;
    static double invoke(double x, void *params) {
      return static_cast<gsl_function_pp*>(params)->_func(x, params);
    }
  };


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
    dof_handler (triangulation),
    fe (FE_Q<DIM>(1), DIM)
  {}




  ElasticProblem::~ElasticProblem ()
  {
    dof_handler.clear ();
    delete mu;
    gsl_integration_workspace_free (w);
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



    if (nonUniform_mesh_flag)
    {
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

      dof_handler.distribute_dofs (fe);
      present_solution.reinit (dof_handler.n_dofs());
    }

    evaluation_point = present_solution;

  }


  void ElasticProblem::get_grad_u1_value_list(const std::vector< Point< DIM > > &  points, std::vector<Tensor<2, DIM>> value_list)
  {

    Tensor<2, DIM> tmp;
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      double x1 = points[point_n](0);
      double x2 = points[point_n](1);

      tmp =  0;

      unsigned int t = 0;
      if(x2 > L1)
        t = 4;

      for(unsigned int i = 0; i < 4; i++)
      {
        double Ae = A[i + t]*exp(alphas[i]*x2);
        double AeB = Ae*B[i];

        tmp[0][0] += Ae;
        tmp[0][1] += Ae*alphas[i];
        tmp[1][0] += AeB;
        tmp[1][1] += AeB*alphas[i];

      }

      tmp[0][0] *= -critical_frequency*cos(critical_frequency*x1);
      tmp[0][1] *= -sin(critical_frequency*x1);
      tmp[0][0] *= -critical_frequency*sin(critical_frequency*x1);
      tmp[0][1] *= cos(critical_frequency*x1);

      value_list[point_n] = tmp;
    }
  }

  void ElasticProblem::get_grad_u2_value_list(const std::vector< Point< DIM > > &  points, std::vector<Tensor<2, DIM>> value_list)
  {
    Tensor<2, DIM> tmp;
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      double x1 = points[point_n](0);
      double x2 = points[point_n](1);

      tmp =  0;

      unsigned int t = 0;
      if(x2 > L1)
        t = 4;

      // homogenous part first
      for(unsigned int i = 0; i < 4; i++)
      {
        double Ce = C[i + t]*exp(r[i]*x2);

        tmp[0][0] += phi1[i]*Ce;
        tmp[0][1] += r[i]*phi1[i]*Ce;
        tmp[1][0] += phi3[i]*Ce;
        tmp[1][1] += r[i]*phi3[i]*Ce;

      }

      // Now do that particular part
      std::vector<double> E1_exp_integral(4);
      std::vector<double> E2_exp_integral(4);

      double a = 0;
      if(x2 > L1)
        a = L1;

      double abserr = 0., relerr = 1.e-7; // requested errors
      double error; // the error estimate

      ElasticProblem* ptr_1 = this;
      auto ptr1 = [=](double x, void *params)->double{return ptr_1->E1exp(x, params);};
      auto ptr2 = [=](double x, void *params)->double{return ptr_1->E2exp(x, params);};

      gsl_function_pp<decltype(ptr1)> Fp1(ptr1);
      gsl_function_pp<decltype(ptr2)> Fp2(ptr2);

      gsl_function *F1 = static_cast<gsl_function*>(&Fp1);
      gsl_function *F2 = static_cast<gsl_function*>(&Fp2);


//      gsl_function F1 = &this->E1exp;
//      gsl_function F2 = &this->E2exp;

      for(unsigned int i = 0; i < 4; i++)
      {
        F1->params = &i;
        F2->params = &i;

        gsl_integration_qag (F1, a, x2, abserr, relerr, np,
            GSL_INTEG_GAUSS15, w, &(E1_exp_integral[i]), &error);
        gsl_integration_qag (F2, a, x2, abserr, relerr, np,
            GSL_INTEG_GAUSS15, w, &(E2_exp_integral[i]), &error);

        E1_exp_integral[i] *= (-1.0/L(1,2,1,2));
        E2_exp_integral[i] *= (-1.0/L(2,2,2,2));
      }

      std::vector<double> exp_integral_eh(4);
      std::vector<double> h(4);
      for(unsigned int i = 0; i < 4; i++)
      {
        exp_integral_eh[i] = exp(r[i]*x2)*
                            (phi_inv_T_2[i]*E1_exp_integral[i] + phi_inv_T_4[i]*E2_exp_integral[i]);

        h[i] = -(1.0/L(1,2,1,2))*phi_inv_T_2[i]*E1(x1) - (1.0/L(2,2,2,2))*phi_inv_T_4[i]*E2(x2);

        tmp[0][0] += phi1[i]*exp_integral_eh[i];
        tmp[0][1] += phi1[i]*(r[i]*exp_integral_eh[i] + h[i]);
        tmp[1][0] += phi3[i]*exp_integral_eh[i];
        tmp[1][1] += phi3[i]*(r[i]*exp_integral_eh[i] + h[i]);
      }


      tmp[1][1] += E2_tilde(x2);

      tmp[0][0] *= 2*w_c*cos(2*w_c*x1);
      tmp[0][1] *= sin(2*w_c*x1);
      tmp[0][0] *= -2*w_c*sin(2*w_c*x1);
      tmp[0][1] *= cos(2*w_c*x1);

      value_list[point_n] = tmp;
    }
  }

  double ElasticProblem::E1(double x2)
  {
    unsigned int t = 0;
    if(x2 > L1)
      t = 4;


    double v1 = 0.0;
    double v1_1 = 0.0;
    double v1_11 = 0.0;
    double v2 = 0.0;
    double v2_1 = 0.0;
    double v2_11 = 0.0;
    for(unsigned int i = 0; i < 4; i++)
    {
      double Ae = A[i + t]*exp(alphas[i]*x2);
      double AeB = Ae*B[i];

      v1 += Ae;
      v1_1 += alphas[i]*Ae;
      v1_11 += alphas[i]*alphas[i]*Ae;

      v2 += AeB;
      v2_1 += alphas[i]*AeB;
      v2_11 += alphas[i]*alphas[i]*AeB;


    }

    double w_c_2 = w_c*w_c;
    double w_c_3 = w_c*w_c_2;

    double out =  -M(1,1,1,1,1,1)*(w_c_3)*v1*v1 + 2*M(1,1,1,1,2,2)*(w_c_2)*v1*v2_1 - M(1,1,2,2,2,2)*(w_c)*v2_1*v2_1 +
            2*M(1,1,1,2,2,1)*(w_c_2)*v1_1*v2 + M(1,1,2,1,2,1)*(w_c_3)*v2*v2 +
            M(1,2,1,1,1,2)*w_c*(v1_1*v1_1 + v1*v1_11) + M(1,2,1,1,2,1)*(w_c_2)*(v1_1*v2 + v1*v2_1) +
            -M(1,2,2,2,1,2)*(v1_11*v2_1 + v1_1*v2_11) - M(1,2,2,2,2,1)*w_c*(v2_1*v2_1 + v2*v2_11);

    if(pieceConstFlag == false)
      out += kappa*(M(1,2,1,1,1,2)*w_c*v1*v1_1 + M(1,2,1,1,2,1)*(w_c_2)*v1*v2 -
                    M(1,2,2,2,1,2)*v1_1*v2_1 - M(1,2,2,2,2,1)*w_c*v2*v2_1);


    return out;
  }

  double ElasticProblem::E2(double x2)
  {

    unsigned int t = 0;
    if(x2 > L1)
      t = 4;


    double v1 = 0.0;
    double v1_1 = 0.0;
    double v1_11 = 0.0;
    double v2 = 0.0;
    double v2_1 = 0.0;
    double v2_11 = 0.0;
    for(unsigned int i = 0; i < 4; i++)
    {
      double Ae = A[i + t]*exp(alphas[i]*x2);
      double AeB = Ae*B[i];

      v1 += Ae;
      v1_1 += alphas[i]*Ae;
      v1_11 += alphas[i]*alphas[i]*Ae;

      v2 += AeB;
      v2_1 += alphas[i]*AeB;
      v2_11 += alphas[i]*alphas[i]*AeB;
    }

    double w_c_2 = w_c*w_c;
    double w_c_3 = w_c*w_c_2;

    double out = 2*(M(2,1,1,1,1,2)*(w_c_2)*v1*v1_1 + M(2,1,1,1,2,1)*(w_c_3)*v1*v2 - M(2,1,2,2,2,1)*(w_c_2)*v2*v2_1 -
                    M(2,1,2,2,1,2)*w_c*v1_1*v2_1) +
                   M(2,2,2,2,2,2)*v2_1*v2_11 - M(2,2,1,1,2,2)*w_c*(v1_1*v2_1 + v1*v2_11) +
                   M(2,2,1,1,1,1)*(w_c_2)*v1*v1_1 - M(2,2,1,2,1,2)*v1_1*v1_11 - M(2,2,2,1,2,1)*(w_c_2)*v2*v2_1 -
                   M(2,2,1,2,2,1)*w_c*(v1_11*v2 + v1_1*v2_1);

    if(pieceConstFlag == false)
      out += kappa*(0.5*M(2,2,2,2,2,2)*v2_1*v2_1 - M(2,2,1,1,2,2)*w_c*v1*v2_1 + 0.5*M(2,2,1,1,1,1)*(w_c_2)*v1*v1 -
                    0.5*M(2,2,1,2,1,2)*v1_1*v1_1 + -0.5*M(2,2,2,1,2,1)*(w_c_2)*v2*v2 - M(2,2,1,2,2,1)*(w_c)*v1_1*v2);


    return out;
  }

  double ElasticProblem::E2_tilde(double x2)
  {

    unsigned int t = 0;
    if(x2 > L1)
      t = 4;


    double v1 = 0.0;
    double v1_1 = 0.0;
    double v2 = 0.0;
    double v2_1 = 0.0;
    for(unsigned int i = 0; i < 4; i++)
    {
      double Ae = A[i + t]*exp(alphas[i]*x2);
      double AeB = Ae*B[i];

      v1 += Ae;
      v1_1 += alphas[i]*Ae;

      v2 += AeB;
      v2_1 += alphas[i]*AeB;
    }

    double w_c_2 = w_c*w_c;

    double out = 0.5*(M(2,2,2,2,2,2)*v2_1*v2_1 - 2*M(2,2,1,1,2,2)*w_c*v1*v2_1 + M(2,2,1,1,1,1)*(w_c_2)*v1*v1 +
               M(2,2,1,2,1,2)*v1_1*v1_1 + M(2,2,2,1,2,1)*(w_c_2)*v2*v2 + 2*M(2,2,1,2,2,1)*w_c*v1_1*v2);

    return out;
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




  void ElasticProblem::assemble_system_energy_and_congugate_lambda()
  {
    double lambda_eval = present_lambda;

    system_energy = 0.0;
    congugate_lambda = 0.0;

    QGauss<DIM>  quadrature_formula(2);

    FEValues<DIM> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();


    std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points,
                                                std::vector<Tensor<1,DIM>>(DIM));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    std::vector<Tensor<1, DIM> > rhs_values (n_q_points);

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

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
      unsigned int nonUniformFlag;
      nonUniformFlag = 0;
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u %u", &grid_dimensions[0], &grid_dimensions[1], &nonUniformFlag);
      if(valuesWritten < 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      else if(valuesWritten == 3)
      {
        if(nonUniformFlag == 1)
          nonUniform_mesh_flag = true;
        else
          nonUniform_mesh_flag = false;
      }
      else
        nonUniform_mesh_flag = false;

      grid_dimensions[0] *= number_unit_cells;

      // read in the absolute tolerance of newton iteration
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg  %u", &tol, &maxIter);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
      }

      // read in exponential growth parameter and possibly the l1 value
      double l1;
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg", &kappa, &l1);
      if(valuesWritten < 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      else if (valuesWritten == 2)
      {
        mu = new MuFunction(kappa, l1);
      }
      else
        mu = new MuFunction(kappa);



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

    // set the domain dimensions
    domain_dimensions[0] = number_unit_cells*2.0*(4.0*atan(1.0))/critical_frequency;
    domain_dimensions[1] = 1.0;

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

}

#endif // COMPRESSEDSTRIPPACABLOCH_CC_
