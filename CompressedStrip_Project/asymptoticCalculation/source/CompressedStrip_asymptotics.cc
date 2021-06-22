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
#include <CompressedStrip_asymptotics.h>
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
      return static_cast<gsl_function_pp*>(params)->_func(x);
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
      present_solution = 0.0;

      vexex_eq.reinit (dof_handler.n_dofs());
      vexex_eq = 0.0;

      u1.reinit (dof_handler.n_dofs());
      u1 = 0.0;
    }

    evaluation_point = present_solution;

  }

  void ElasticProblem::check_W_derivs()
  {

    double eps = 1e-7;
    Tensor<2,DIM> F;
    F[0][0] = 0.2345;
    F[0][1] = 0.1723;
    F[1][0] = 0.0543;
    F[1][1] = 0.12223;

    Tensor<2, DIM> F_inv = invert(F);
    double II_F = determinant(F);

    Tensor<2, DIM> dW_dF = nh.get_piola_kirchoff_tensor(0.33, 1.0, F, F_inv, II_F);
    Tensor<4,DIM> d2W_dF = nh.get_incremental_moduli_tensor(0.33, 1.0, F_inv, II_F);
    Tensor<6,DIM> d3W_dF = nh.get_d3W_dFdFdF(0.33, 1.0, F_inv, II_F);
    Tensor<8,DIM> d4W_dF = nh.get_d4W_dFdFdFdF(0.33, 1.0, F_inv, II_F);


    Tensor<4, DIM> d2_numer;
    // check d3w
    for(unsigned int i = 0; i < DIM ;i ++)
      for(unsigned int j = 0; j < DIM; j++)
      {
        F[i][j] += eps;
        Tensor<2, DIM> F_inv_ = invert(F);
        double II_F_ = determinant(F);
        Tensor<2,DIM> dW_dF_new = nh.get_piola_kirchoff_tensor(0.33, 1.0, F, F_inv_, II_F_);
        F[i][j] -= eps;


        for(unsigned int k = 0; k < DIM ;k ++)
          for(unsigned int l = 0; l < DIM; l++)
              {
                d2_numer[k][l][i][j] = ((dW_dF_new[k][l] - dW_dF[k][l])/eps - d2W_dF[k][l][i][j])/d2W_dF[k][l][i][j];
              }
      }


    Tensor<6, DIM> d3_numer;
    // check d3w
    for(unsigned int i = 0; i < DIM ;i ++)
      for(unsigned int j = 0; j < DIM; j++)
      {
        F[i][j] += eps;
        Tensor<2, DIM> F_inv_ = invert(F);
        double II_F_ = determinant(F);
        F[i][j] -= eps;

        Tensor<4,DIM> d2W_dF_new = nh.get_incremental_moduli_tensor(0.33, 1.0, F_inv_, II_F_);

        for(unsigned int k = 0; k < DIM ;k ++)
          for(unsigned int l = 0; l < DIM; l++)
            for(unsigned int m = 0; m < DIM ;m ++)
              for(unsigned int n = 0; n < DIM; n++)
              {
                d3_numer[k][l][m][n][i][j] = ((d2W_dF_new[k][l][m][n] - d2W_dF[k][l][m][n])/eps - d3W_dF[k][l][m][n][i][j])/d3W_dF[k][l][m][n][i][j];
              }
      }



    // check d4w
    Tensor<8, DIM> d4_numer;
    for(unsigned int i = 0; i < DIM ;i ++)
      for(unsigned int j = 0; j < DIM; j++)
      {
        F[i][j] += eps;
        Tensor<2, DIM> F_inv_ = invert(F);
        double II_F_ = determinant(F);
        F[i][j] -= eps;

        Tensor<6,DIM> d3W_dF_new = nh.get_d3W_dFdFdF(0.33, 1.0, F_inv_, II_F_);

        for(unsigned int k = 0; k < DIM ;k ++)
          for(unsigned int l = 0; l < DIM; l++)
            for(unsigned int m = 0; m < DIM ;m ++)
              for(unsigned int n = 0; n < DIM; n++)
                for(unsigned int p = 0; p < DIM ; p ++)
                  for(unsigned int q = 0; q < DIM; q++)
                  {
                    d4_numer[k][l][m][n][p][q][i][j] = ((d3W_dF_new[k][l][m][n][p][q] - d3W_dF[k][l][m][n][p][q])/eps - d4W_dF[k][l][m][n][p][q][i][j])/(d4W_dF[k][l][m][n][p][q][i][j]);
                  }
      }


    std::cout << "Norm of the difference for d2W : " << d2_numer.norm() << std::endl;
    std::cout << "Norm of the difference for d3W : " << d3_numer.norm() << std::endl;
    std::cout << "Norm of the difference for d4W : " << d4_numer.norm() << std::endl;

  }

  void ElasticProblem::get_grad_u1_value_list(const std::vector< Point< DIM > > &  points, std::vector<Tensor<2, DIM>> & value_list)
  {

    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      Tensor<2, DIM> tmp;
      Tensor<2, DIM, std::complex<double>> tmp_complex;

      double x1 = points[point_n](0);
      double x2 = points[point_n](1);

      unsigned int t = 0;
      if(x2 > L1)
        t = 4;

      for(unsigned int i = 0; i < 4; i++)
      {
        std::complex<double> Ae = A[i + t]*std::exp(alphas[i]*x2);
        std::complex<double> AeB = Ae*B[i];

        tmp_complex[0][0] += Ae;
        tmp_complex[0][1] += Ae*alphas[i];
        tmp_complex[1][0] += AeB;
        tmp_complex[1][1] += AeB*alphas[i];
      }

      tmp[0][0] = -critical_frequency*cos(critical_frequency*x1)*(tmp_complex[0][0]).real();
      tmp[0][1] = -sin(critical_frequency*x1)*(tmp_complex[0][1]).real();
      tmp[1][0] = -critical_frequency*sin(critical_frequency*x1)*(tmp_complex[1][0]).real();
      tmp[1][1] = cos(critical_frequency*x1)*(tmp_complex[1][1]).real();

      value_list[point_n] = tmp;
    }
  }

  void ElasticProblem::get_grad_u2_value_list(const std::vector< Point< DIM > > &  points, std::vector<Tensor<2, DIM>> & value_list)
  {
    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
    {
      Tensor<2, DIM> tmp;
      Tensor<2, DIM, std::complex<double>> tmp_complex;
      double x1 = points[point_n](0);
      double x2 = points[point_n](1);

      unsigned int t = 0;
      if(x2 > L1)
        t = 4;

      // homogenous part first
      for(unsigned int i = 0; i < 4; i++)
      {
        std::complex<double> Ce = C[i + t]*std::exp(r[i]*x2);

        tmp_complex[0][0] += phi1[i]*Ce;
        tmp_complex[0][1] += r[i]*phi1[i]*Ce;
        tmp_complex[1][0] += phi3[i]*Ce;
        tmp_complex[1][1] += r[i]*phi3[i]*Ce;
      }

      // Now do that particular part
      std::vector<double> E1_exp_integral_real(4);
      std::vector<double> E2_exp_integral_real(4);

      std::vector<double> E1_exp_integral_imag(4);
      std::vector<double> E2_exp_integral_imag(4);

      double a = 0;
      if(x2 > L1 && pieceConstFlag == true)
        a = L1;

      double abserr = 1.e-7, relerr = 1.e-7; // requested errors
      double error; // the error estimate

      ElasticProblem* ptr_1 = this;
      auto ptr1 = [=](double x)->double{return ptr_1->E1exp_real(x);};
      auto ptr2 = [=](double x)->double{return ptr_1->E2exp_real(x);};
      auto ptr3 = [=](double x)->double{return ptr_1->E1exp_imag(x);};
      auto ptr4 = [=](double x)->double{return ptr_1->E2exp_imag(x);};

      gsl_function_pp<decltype(ptr1)> Fp1(ptr1);
      gsl_function_pp<decltype(ptr2)> Fp2(ptr2);
      gsl_function_pp<decltype(ptr3)> Fp3(ptr3);
      gsl_function_pp<decltype(ptr4)> Fp4(ptr4);

      gsl_function *F1 = static_cast<gsl_function*>(&Fp1);
      gsl_function *F2 = static_cast<gsl_function*>(&Fp2);
      gsl_function *F3 = static_cast<gsl_function*>(&Fp3);
      gsl_function *F4 = static_cast<gsl_function*>(&Fp4);

      for(unsigned int i = 0; i < 4; i++)
      {

        root_index = i;

        gsl_integration_qag (F1, a, x2, abserr, relerr, np,
            GSL_INTEG_GAUSS31, w, &(E1_exp_integral_real[i]), &error);
        gsl_integration_qag (F2, a, x2, abserr, relerr, np,
            GSL_INTEG_GAUSS31, w, &(E2_exp_integral_real[i]), &error);
        gsl_integration_qag (F3, a, x2, abserr, relerr, np,
            GSL_INTEG_GAUSS31, w, &(E1_exp_integral_imag[i]), &error);
        gsl_integration_qag (F4, a, x2, abserr, relerr, np,
            GSL_INTEG_GAUSS31, w, &(E2_exp_integral_imag[i]), &error);


        E1_exp_integral_real[i] *= (-1.0/L(1,2,1,2));
        E2_exp_integral_real[i] *= (-1.0/L(2,2,2,2));
        E1_exp_integral_imag[i] *= (-1.0/L(1,2,1,2));
        E2_exp_integral_imag[i] *= (-1.0/L(2,2,2,2));
      }

      std::vector<std::complex<double>> exp_integral_eh(4);
      std::vector<std::complex<double>> h(4);

      std::vector<std::complex<double>> E1_exp_integral(4);
      std::vector<std::complex<double>> E2_exp_integral(4);
      for(unsigned int i = 0; i < 4; i++)
      {
        E1_exp_integral[i] = std::complex<double>(E1_exp_integral_real[i], E1_exp_integral_imag[i]);
        E2_exp_integral[i] = std::complex<double>(E2_exp_integral_real[i], E2_exp_integral_imag[i]);
        exp_integral_eh[i] = std::exp(r[i]*x2)*
                             (phi_inv_T_2[i]*E1_exp_integral[i] + phi_inv_T_4[i]*E2_exp_integral[i]);

        h[i] = -(1.0/L(1,2,1,2))*phi_inv_T_2[i]*E1(x2) - (1.0/L(2,2,2,2))*phi_inv_T_4[i]*E2(x2);

        tmp_complex[0][0] += phi1[i]*exp_integral_eh[i];
        tmp_complex[0][1] += phi1[i]*(r[i]*exp_integral_eh[i] + h[i]);
        tmp_complex[1][0] += phi3[i]*exp_integral_eh[i];
        tmp_complex[1][1] += phi3[i]*(r[i]*exp_integral_eh[i] + h[i]);

      }

      tmp[0][0] = 2*w_c*cos(2*w_c*x1)*(tmp_complex[0][0]).real();
      tmp[0][1] = sin(2*w_c*x1)*(tmp_complex[0][1]).real();
      tmp[1][0] = -2*w_c*sin(2*w_c*x1)*(tmp_complex[1][0]).real();
      tmp[1][1] = cos(2*w_c*x1)*(tmp_complex[1][1]).real() - (1.0/L(2,2,2,2))*E2_tilde(x2);

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
      std::complex<double> Ae = A[i + t]*std::exp(alphas[i]*x2);
      std::complex<double> AeB = Ae*B[i];

      v1 += (Ae).real();
      v1_1 += (alphas[i]*Ae).real();
      v1_11 += (alphas[i]*alphas[i]*Ae).real();

      v2 += (AeB).real();
      v2_1 += (alphas[i]*AeB).real();
      v2_11 += (alphas[i]*alphas[i]*AeB).real();

    }

    double w_c_2 = w_c*w_c;
    double w_c_3 = w_c*w_c_2;

    double out =  -M(1,1,1,1,1,1)*(w_c_3)*v1*v1 + 2.0*M(1,1,1,1,2,2)*(w_c_2)*v1*v2_1 - M(1,1,2,2,2,2)*(w_c)*v2_1*v2_1 +
            2.0*M(1,1,1,2,2,1)*(w_c_2)*v1_1*v2 + M(1,1,2,1,2,1)*(w_c_3)*v2*v2 +
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
      std::complex<double> Ae = A[i + t]*std::exp(alphas[i]*x2);
      std::complex<double> AeB = Ae*B[i];

      v1 += (Ae).real();
      v1_1 += (alphas[i]*Ae).real();
      v1_11 += (alphas[i]*alphas[i]*Ae).real();

      v2 += (AeB).real();
      v2_1 += (alphas[i]*AeB).real();
      v2_11 += (alphas[i]*alphas[i]*AeB).real();
    }

    double w_c_2 = w_c*w_c;
    double w_c_3 = w_c*w_c_2;

    double out = 2.0*(M(2,1,1,1,1,2)*(w_c_2)*v1*v1_1 + M(2,1,1,1,2,1)*(w_c_3)*v1*v2 - M(2,1,2,2,2,1)*(w_c_2)*v2*v2_1 -
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
      std::complex<double> Ae = A[i + t]*std::exp(alphas[i]*x2);
      std::complex<double> AeB = Ae*B[i];

      v1 += (Ae).real();
      v1_1 += (alphas[i]*Ae).real();

      v2 += (AeB).real();
      v2_1 += (alphas[i]*AeB).real();
    }

    double w_c_2 = w_c*w_c;

    double out = 0.5*(M(2,2,2,2,2,2)*v2_1*v2_1 - 2.0*M(2,2,1,1,2,2)*w_c*v1*v2_1 + M(2,2,1,1,1,1)*(w_c_2)*v1*v1 +
               M(2,2,1,2,1,2)*v1_1*v1_1 + M(2,2,2,1,2,1)*(w_c_2)*v2*v2 + 2.0*M(2,2,1,2,2,1)*w_c*v1_1*v2);

    return out;
  }
  void ElasticProblem::update_F0(const double lambda)
  {
    // Update the uniform diagonal F0 tensor.

    Assert ((lambda >= 0 && lambda <= 1.0), ExcInternalError());

    double lambda1 = 1.0 - lambda;
    double lambda2;

    double nuVal = nu.get_nu_value();
    // solve quadratic for lambda2
    double a = (1.0 + (2.0*nuVal/(1.0 - nuVal))*lambda1*lambda1);
    double b = (-2.0*nuVal/(1.0 - nuVal))*lambda1;
    double c = -1.0;


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

        double nuVal = nu.get_nu_value();

        a_nu = (2.0*nuVal/(1.0 - nuVal));
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

  void ElasticProblem::assemble_asymptotic_integrals()
  {
    double lambda_eval = present_lambda;

    E_u1u1u1u1 = 0.0;
    E_u2u1u1 = 0.0;
    dEdlambda_u1u1 = 0.0;
    E_u1u1u1 = 0.0;
    E_u2u2 = 0.0;

    QGauss<DIM>  quadrature_formula(5);

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

    std::vector<Tensor<2, DIM>> grad_u1(n_q_points);
    std::vector<Tensor<2, DIM>> grad_u2(n_q_points);

    // get the dF_dlambda
    double dlambda1_dlambda, dlambda2_dlambda, a_nu;

    dlambda1_dlambda = -1.0;

    double nuVal = nu.get_nu_value();
    a_nu = (2.0*nuVal/(1.0 - nuVal));
    double lambda1 = 1.0 - lambda_eval;
    double term1 = 2.0*a_nu*(1.0 + a_nu*lambda1*lambda1);
    double term2 = -4.0*a_nu*a_nu*lambda1*lambda1;
    double term3 = (2.0*a_nu*a_nu*lambda1 + 8.0*a_nu*lambda1)*(1.0 + a_nu * lambda1*lambda1)/
                          sqrt(a_nu*a_nu*lambda1*lambda1 + 4.0 *a_nu*lambda1*lambda1 + 4.0);
    double term4 = -sqrt(a_nu*a_nu*lambda1*lambda1 + 4.0 *a_nu*lambda1*lambda1 + 4.0)*4.0*a_nu*lambda1;

    dlambda2_dlambda = -(term1 + term2 + term3 + term4)/(4.0*(1.0 + a_nu*lambda1*lambda1)*(1.0 + a_nu*lambda1*lambda1));

    Tensor<2, DIM> dF_dlambda;
    dF_dlambda[0][0] = dlambda1_dlambda;
    dF_dlambda[0][1] = 0.0;
    dF_dlambda[1][0] = 0.0;
    dF_dlambda[1][1] = dlambda2_dlambda;

//    std::cout << "\n\n" << lambda_eval << "   " << dlambda2_dlambda << std::endl;
//
//
//    std::vector<Point<DIM> > pp(1);
//    Point<DIM> ppp(0.1248, 0.432);
//    pp[0] = ppp;
//    get_grad_u1_value_list(pp, grad_u2);
//    std::cout << (grad_u2[0])[0][0] << " " <<  (grad_u2[0])[0][1] << "  " << (grad_u2[0])[1][0] << " " <<  (grad_u2[0])[1][1] << std::endl;
//
//    Point<DIM> pp2(0.1248, 0.932);
//    pp[0] = pp2;
//    get_grad_u1_value_list(pp, grad_u2);
//    std::cout << (grad_u2[0])[0][0] << " " <<  (grad_u2[0])[0][1] << "  " << (grad_u2[0])[1][0] << " " <<  (grad_u2[0])[1][1] << std::endl;
//    exit(-1);
    double u1u1 = 0.0;

    typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);

      get_grad_u1_value_list(fe_values.get_quadrature_points(), grad_u1);
      get_grad_u2_value_list(fe_values.get_quadrature_points(), grad_u2);

      std::vector<Vector<double>> u1(n_q_points, Vector<double>(DIM));

      u1_value_list(fe_values.get_quadrature_points(), u1);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        u1u1 += (u1[q_point])*(u1[q_point])*fe_values.JxW(q_point);


        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<4,DIM> d2W_dF = nh.get_incremental_moduli_tensor(nu_values[q_point], mu_values[q_point], F_inv, II_F);
        Tensor<6,DIM> d3W_dF = nh.get_d3W_dFdFdF(nu_values[q_point], mu_values[q_point], F_inv, II_F);
        Tensor<8,DIM> d4W_dF = nh.get_d4W_dFdFdFdF(nu_values[q_point], mu_values[q_point], F_inv, II_F);

        double contrib_E_u1u1u1 = 0.0;
        double contrib_E_u1u1u1u1 = 0.0;
        double contrib_E_u2u1u1 = 0.0;
        double contrib_dEdlambda_u1u1 = 0.0;
        double contrib_E_u2u2 = 0.0;
        for (unsigned int i=0; i<DIM; ++i)
          for (unsigned int j=0; j<DIM; ++j)
            for (unsigned int k=0; k<DIM; ++k)
              for (unsigned int l=0; l<DIM; ++l)
              {
                contrib_E_u2u2 += d2W_dF[i][j][k][l]*(grad_u2[q_point])[i][j]*(grad_u2[q_point])[k][l];
                for (unsigned int m=0; m<DIM; ++m)
                  for (unsigned int n=0; n<DIM; ++n)
                  {
                    contrib_E_u2u1u1 +=
                        d3W_dF[i][j][k][l][m][n]*(grad_u2[q_point])[i][j]*(grad_u1[q_point])[k][l]*(grad_u1[q_point])[m][n];

                    contrib_E_u1u1u1 +=
                        d3W_dF[i][j][k][l][m][n]*(grad_u1[q_point])[i][j]*(grad_u1[q_point])[k][l]*(grad_u1[q_point])[m][n];

                    contrib_dEdlambda_u1u1 +=
                        d3W_dF[i][j][k][l][m][n]*dF_dlambda[m][n]*(grad_u1[q_point])[i][j]*(grad_u1[q_point])[k][l];
                    for (unsigned int p=0; p<DIM; ++p)
                      for (unsigned int q=0; q<DIM; ++q)
                      {
                        contrib_E_u1u1u1u1 +=
                            d4W_dF[i][j][k][l][m][n][p][q]*
                               (grad_u1[q_point])[i][j]*(grad_u1[q_point])[k][l]*
                               (grad_u1[q_point])[m][n]*(grad_u1[q_point])[p][q];

                      }
                  }
              }

        E_u1u1u1u1 += contrib_E_u1u1u1u1*fe_values.JxW(q_point);
        E_u2u1u1 += contrib_E_u2u1u1*fe_values.JxW(q_point);
        dEdlambda_u1u1 += contrib_dEdlambda_u1u1*fe_values.JxW(q_point);

        E_u2u2 += contrib_E_u2u2*fe_values.JxW(q_point);

        E_u1u1u1 += contrib_E_u1u1u1*fe_values.JxW(q_point);
      }
    }
    u1u1 = u1u1/domain_dimensions[0];

    delta_lambda1 = -0.5*E_u1u1u1/dEdlambda_u1u1;
    delta_lambda2 = -(1.0/3.0)*(E_u1u1u1u1 + 3.0*E_u2u1u1)/(dEdlambda_u1u1);

    std::cout << E_u1u1u1u1 << "     " <<  E_u2u1u1 << "     " << dEdlambda_u1u1  << "   " << u1u1 << std::endl;
    // now get the congateLambda2

    delta_conguateLambda2 = 0.0;


    // need to get d2lambda2_dlambdadlambda

    double d2l2_dldl = 0.0;
    double lambda12 = lambda1*lambda1;
    double eta = a_nu;
    double eta2 = eta*eta;
    double phi = eta2 + 4.0*eta;
    double a = 1.0 + eta*lambda12;
    double a_inv = 1.0/a;
    double b = phi*lambda12 + 4.0;
    double bsqrt = sqrt(b);
    double bsqrt_inv = 1.0/bsqrt;

    d2l2_dldl = -3.0*eta2*lambda1*a_inv*a_inv + 4.0*eta2*eta*lambda12*lambda1*a_inv*a_inv*a_inv +
        0.5*phi*bsqrt_inv*a_inv - 0.5*phi*phi*lambda12*a_inv*bsqrt_inv*bsqrt_inv*bsqrt_inv -
        phi*eta*lambda12*bsqrt_inv*a_inv*a_inv - eta*bsqrt*a_inv*a_inv - eta*phi*lambda12*bsqrt_inv*a_inv*a_inv +
        4.0*eta2*lambda12*bsqrt*a_inv*a_inv*a_inv;

    Tensor<2, DIM> dF2_dlambda;
    dF2_dlambda[0][0] = 0.0;
    dF2_dlambda[0][1] = 0.0;
    dF2_dlambda[1][0] = 0.0;
    dF2_dlambda[1][1] = d2l2_dldl;


    typename DoFHandler<DIM>::active_cell_iterator cell2 = dof_handler.begin_active(),
                                                   endc2 = dof_handler.end();
    for (; cell2!=endc2; ++cell2)
    {

      fe_values.reinit (cell2);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu->value_list     (fe_values.get_quadrature_points(), mu_values);

//      get_grad_u1_value_list(fe_values.get_quadrature_points(), grad_u1);
      get_grad_u2_value_list(fe_values.get_quadrature_points(), grad_u2);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
        Tensor<2, DIM> F_inv = invert(F);
        double II_F = determinant(F);

        Tensor<4,DIM> d2W_dF = nh.get_incremental_moduli_tensor(nu_values[q_point], mu_values[q_point], F_inv, II_F);
//        Tensor<2,DIM> dW_dF = nh.get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point], F, F_inv, II_F);

        double tmp = 0.0;
//        tmp +=  delta_lambda2*dW_dF[1][1]*dF2_dlambda[1][1];
        for (unsigned int i=0; i<DIM; ++i)
          for (unsigned int j=0; j<DIM; ++j)
            for (unsigned int k=0; k<DIM; ++k)
              for (unsigned int l=0; l<DIM; ++l)
              {
                tmp += d2W_dF[i][j][k][l]*dF_dlambda[i][j]*((grad_u2[q_point])[k][l] + 1.0*delta_lambda2*dF_dlambda[k][l]);
              }



        delta_conguateLambda2 += tmp*fe_values.JxW(q_point);

      }

    }

    delta_conguateLambda2 += dEdlambda_u1u1;
  }

  void ElasticProblem::assemble_vexex_eq()
   {
     // Assembling the system rhs. I choose to make the rhs and system matrix assemblies separate,
     // because we only do one at a time anyways in the newton method.


    vexex_eq = 0.0;

     QGauss<DIM>  quadrature_formula(10);

     FEValues<DIM> fe_values (fe, quadrature_formula,
                              update_values   | update_gradients |
                              update_quadrature_points | update_JxW_values);

     const unsigned int   dofs_per_cell = fe.dofs_per_cell;
     const unsigned int   n_q_points    = quadrature_formula.size();

     Vector<double>       cell_rhs (dofs_per_cell);

     std::vector<std::vector<Tensor<1,DIM> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,DIM>>(DIM));

     std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

     std::vector<double>     nu_values (n_q_points);
     std::vector<double>     mu_values (n_q_points);

     std::vector<Tensor<2, DIM>> grad_u1(n_q_points);
     std::vector<Tensor<2, DIM>> grad_u2(n_q_points);


     typename DoFHandler<DIM>::active_cell_iterator cell = dof_handler.begin_active(),
                                                    endc = dof_handler.end();
     for (; cell!=endc; ++cell)
     {
       cell_rhs = 0.0;

       fe_values.reinit (cell);

       fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

       nu.value_list (fe_values.get_quadrature_points(), nu_values);
       mu->value_list     (fe_values.get_quadrature_points(), mu_values);

       get_grad_u1_value_list(fe_values.get_quadrature_points(), grad_u1);
       get_grad_u2_value_list(fe_values.get_quadrature_points(), grad_u2);

       for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
       {

         Tensor<2,DIM> F = get_deformation_gradient(old_solution_gradients[q_point]);
         Tensor<2, DIM> F_inv = invert(F);
         double II_F = determinant(F);

         Tensor<4,DIM> d2W_dF = nh.get_incremental_moduli_tensor(nu_values[q_point], mu_values[q_point], F_inv, II_F);
         Tensor<6,DIM> d3W_dF = nh.get_d3W_dFdFdF(nu_values[q_point], mu_values[q_point], F_inv, II_F);

         Tensor<2, DIM> to_contract_with_delta_grad_v;

         for (unsigned int i=0; i<DIM; ++i)
           for (unsigned int j=0; j<DIM; ++j)
             for (unsigned int k=0; k<DIM; ++k)
               for (unsigned int l=0; l<DIM; ++l)
               {
                 to_contract_with_delta_grad_v[i][j] += d2W_dF[k][l][i][j]*(grad_u2[q_point])[k][l];

                 for (unsigned int m=0; m<DIM; ++m)
                   for (unsigned int n=0; n<DIM; ++n)
                   {
                     to_contract_with_delta_grad_v[i][j] += d3W_dF[m][n][k][l][i][j]*(grad_u1[q_point])[k][l]*(grad_u1[q_point])[m][n];
                   }
               }

         for (unsigned int q = 0; q < dofs_per_cell; ++q)
         {
           const unsigned int component_n = fe.system_to_component_index(q).first;

           for(unsigned int p = 0; p<DIM; ++p)
           {
             cell_rhs(q) += to_contract_with_delta_grad_v[component_n][p]*(fe_values.shape_grad(q, q_point)[p])*fe_values.JxW(q_point);

           }

         }
       }

       cell->get_dof_indices (local_dof_indices);

       for (unsigned int n=0; n<dofs_per_cell; ++n)
         vexex_eq(local_dof_indices[n]) += cell_rhs(n);

     }

     double u1_dot_vexex_eq = 0.0;

     std::vector<Vector<double>> solution_vals(n_q_points,
                                               Vector<double>(DIM));

     std::vector<Vector<double>> u1_1(n_q_points, Vector<double>(DIM));

     typename DoFHandler<DIM>::active_cell_iterator cell2 = dof_handler.begin_active(),
                                                  endc2 = dof_handler.end();
     for (; cell2!=endc2; ++cell2)
     {
       fe_values.reinit (cell2);

       fe_values.get_function_values(vexex_eq, solution_vals);

       u1_value_list(fe_values.get_quadrature_points(), u1_1);

       for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
       {
         u1_dot_vexex_eq += solution_vals[q_point]*(u1_1[q_point])*fe_values.JxW(q_point);
 //        u1u1 += (u1[q_point])*(u1[q_point])*fe_values.JxW(q_point);
       }

     }

//     std::cout << "COWS : " << u1_dot_vexex_eq << "\n\n\n\n\n\n\n";

     get_u1();
     u1 *= u1_dot_vexex_eq;
//     vexex_eq -= u1;

     unsigned int number_dofs = dof_handler.n_dofs();

      // get boundary 3 x2 dofs

      std::vector<bool> x2_components = {false, true};
      ComponentMask x2_mask(x2_components);

      std::vector<bool> boundary_3_dof_x2(number_dofs, false);

      std::set< types::boundary_id > boundary_id_3;
      boundary_id_3.insert(3);

      DoFTools::extract_boundary_dofs(dof_handler, x2_mask, boundary_3_dof_x2, boundary_id_3);

      for(unsigned int i = 0; i < number_dofs; i++)
      {
//        if(boundary_3_dof_x2[i])
//            vexex_eq[i] = 0.0;

      }
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

  void ElasticProblem::get_u1()
  {
    unsigned int number_dofs = dof_handler.n_dofs();

    // get the coords of the dofs
    std::vector<Point<DIM>> support_points(number_dofs);
    MappingQ1<DIM> mapping;
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> is_x1_comp(number_dofs, false);

    DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);

    for(unsigned int i = 0; i < number_dofs; i++)
    {
      double x1 = support_points[i](0);
      double x2 = support_points[i](1);

      unsigned int t = 0;
      if(x2 > L1)
        t = 4;

      double currentVal = 0.0;
      for(unsigned int i = 0; i < 4; i++)
      {
        std::complex<double> Ae = A[i + t]*std::exp(alphas[i]*x2);
        std::complex<double> AeB = Ae*B[i];

        if(is_x1_comp[i])
          currentVal += (Ae).real();
        else
          currentVal += (AeB).real();

      }

      if(is_x1_comp[i])
        u1[i] = -sin(w_c*x1)*currentVal;
      else
        u1[i] = cos(w_c*x1)*currentVal;
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
    // Now output the displacements from uniform solution
     std::string filename1("output");
     filename1 += "/vexex_eq";

     // see if the directory exists...
     struct stat st;
     if (stat(filename1.c_str(), &st) == -1)
       mkdir(filename1.c_str(), 0700);

     filename1 += "/vexex_eq-";
     filename1 += std::to_string(cycle);
     filename1 += ".vtk";
     std::ofstream output_vexex_eq (filename1.c_str());

     DataOut<DIM> data_out_vexex_eq;

     data_out_vexex_eq.attach_dof_handler (dof_handler);

     data_out_vexex_eq.add_data_vector (vexex_eq, solution_names);
     data_out_vexex_eq.build_patches ();
     data_out_vexex_eq.write_vtk (output_vexex_eq);

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
//      char num_cells_str[MAXLINE];
//      sprintf(num_cells_str, "%u", number_unit_cells);
//      strcat(directory_name, "_");
//      strcat(directory_name, num_cells_str);

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
        nu.set_nu_value(nuVal);

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
        nu.set_nu_value(nuVal);

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

      // read in the phi1
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
          &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      phi1.resize(4);
      phi1[0] = std::complex<double>(re1, im1);
      phi1[1] = std::complex<double>(re2, im2);
      phi1[2] = std::complex<double>(re3, im3);
      phi1[3] = std::complex<double>(re4, im4);

      // read in the phi3
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
          &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      phi3.resize(4);
      phi3[0] = std::complex<double>(re1, im1);
      phi3[1] = std::complex<double>(re2, im2);
      phi3[2] = std::complex<double>(re3, im3);
      phi3[3] = std::complex<double>(re4, im4);

      // read in the phi_inv_T_2
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
          &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      phi_inv_T_2.resize(4);
      phi_inv_T_2[0] = std::complex<double>(re1, im1);
      phi_inv_T_2[1] = std::complex<double>(re2, im2);
      phi_inv_T_2[2] = std::complex<double>(re3, im3);
      phi_inv_T_2[3] = std::complex<double>(re4, im4);

      // read in the phi_inv_T_4
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
          &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      phi_inv_T_4.resize(4);
      phi_inv_T_4[0] = std::complex<double>(re1, im1);
      phi_inv_T_4[1] = std::complex<double>(re2, im2);
      phi_inv_T_4[2] = std::complex<double>(re3, im3);
      phi_inv_T_4[3] = std::complex<double>(re4, im4);

      // read in the r's
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg %lg %lg %lg %lg",
          &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4);
      if(valuesWritten != 8)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      r.resize(4);
      r[0] = std::complex<double>(re1, im1);
      r[1] = std::complex<double>(re2, im2);
      r[2] = std::complex<double>(re3, im3);
      r[3] = std::complex<double>(re4, im4);

      r_real.resize(4);
      r_real[0] = re1;
      r_real[1] = re2;
      r_real[2] = re3;
      r_real[3] = re4;

      r_imag.resize(4);
      r_imag[0] = im1;
      r_imag[1] = im2;
      r_imag[2] = im3;
      r_imag[3] = im4;


      // read in the amplitude C's
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
        C.resize(8);
        C[4] = std::complex<double>(re5, im5);
        C[5] = std::complex<double>(re6, im6);
        C[6] = std::complex<double>(re7, im7);
        C[7] = std::complex<double>(re8, im8);
      }
      else
        C.resize(4);

      C[0] = std::complex<double>(re1, im1);
      C[1] = std::complex<double>(re2, im2);
      C[2] = std::complex<double>(re3, im3);
      C[3] = std::complex<double>(re4, im4);

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

    // set our constants and stuff
    set_present_lambda(critical_lambda_analytical);
    Tensor<2, DIM> F_inv = invert(F0);
    double II_F = determinant(F0);
    L_tensor = nh.get_incremental_moduli_tensor(nu.get_nu_value(), 1.0, F_inv, II_F);
    M_tensor = nh.get_d3W_dFdFdF(nu.get_nu_value(), 1.0, F_inv, II_F);


//    grid_dimensions[0] = 20;
//    grid_dimensions[1] = 20;

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
