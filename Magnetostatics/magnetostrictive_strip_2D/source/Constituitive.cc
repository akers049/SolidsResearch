#ifndef CONSTITUITIVE_CC_
#define CONSTITUITIVE_CC_

#include "Constituitive.h"

#define DIM 2
#define MU0 1.0

using namespace dealii;

  // Compressible Neo hookean

  inline
  double Compressible_NeoHookean::get_energy(const double nu, const double mu,
                              Tensor<2, DIM> F, double II_F)
  {
    double I_C = F[1][0]*F[1][0] + F[0][0]*F[0][0] + F[0][1]*F[0][1] + F[1][1]*F[1][1];
    double W = mu*(0.5*(I_C - 2 - log(II_F*II_F)) + (nu/(1.0- nu))*(II_F - 1.0)*(II_F - 1.0));

    return W;
  }


  inline
  Tensor<2,DIM> Compressible_NeoHookean::get_piola_kirchoff_tensor(const double nu, const double mu,
                                            Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                            double II_F)
  {
    Tensor<2, DIM> tmp;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
      {
        tmp[i][j] = mu*F[i][j] - mu*F_inv[j][i] +
                    (2.0*mu*nu/(1.0- nu))*(II_F*II_F - II_F)*F_inv[j][i];
      }

    return tmp;

  }

  inline
  Tensor<4,DIM> Compressible_NeoHookean::get_incremental_moduli_tensor(const double nu,
                                                const double mu, Tensor<2,DIM> F_inv,
                                                double II_F)
  {

    Tensor<4,DIM> tmp;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int k=0; k<DIM; ++k)
          for (unsigned int l=0; l<DIM; ++l)
          {
            tmp[i][j][k][l] = ((i==k) && (j==l) ? mu : 0.0) +
                (mu - (2.0*mu*nu/(1.0-nu))*(II_F*II_F - II_F))*F_inv[j][k]*F_inv[l][i] +
                (4.0*nu*mu/(1.0-nu))*(II_F*II_F - 0.5*II_F)*F_inv[l][k]*F_inv[j][i];
          }

    return tmp;
  }

  // Magnetic one

 double Compressible_NH_Langevin::get_energy(const bool rho, const double nu, const double mu,
     const double suseptibility, const double msf,
     Tensor<2,DIM> F, Tensor<1, DIM> B)
 {
   double II_F = determinant(F);
   double I_C = F[1][0]*F[1][0] + F[0][0]*F[0][0] + F[0][1]*F[0][1] + F[1][1]*F[1][1];

   double W = 0.5*(1.0/(II_F*MU0))*((F*B)*(F*B));

   if(rho == true)
   {
     // elastic part
     W += mu*(0.5*(I_C - 2 - log(II_F*II_F)) + (nu/(1.0- nu))*(II_F - 1.0)*(II_F - 1.0));

     if(fabs(suseptibility) >= 1e-6)
     {
       // magnetic part

       double alpha = MU0*msf/(3.0*suseptibility);
       double over_alpha = 1.0/alpha;
       double norm_sqr_FB = (F*B)*(F*B);
       double mag_b = sqrt(norm_sqr_FB)/II_F;

       W += II_F*alpha*msf*(log(over_alpha*mag_b) - log(sinh(over_alpha*mag_b)));
     }
   }

   return W;
 }

  Tensor<1,DIM> Compressible_NH_Langevin::get_dW_dB(const double nu, const double mu,
                                  const double suseptibility, const double msf,
                                  Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                  double II_F, Tensor<1, DIM> B)
  {
    Tensor<1, DIM> tmp;
    Tensor<1, DIM> FB = F*B;

    tmp = 2.0*msf*suseptibility*(transpose(F)*FB);

    return tmp;

  }

  Tensor<2,DIM> Compressible_NH_Langevin::get_d2W_dBdB(const double nu, const double mu,
                                  const double suseptibility, const double msf,
                                  Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                  double II_F, Tensor<1, DIM> B)
  {
    Tensor<2, DIM> tmp;
    tmp = 2.0*msf*suseptibility*(transpose(F)*F);

    return tmp;

  }

  Tensor<2,DIM> Compressible_NH_Langevin::get_piola_kirchoff_tensor(const double nu, const double mu,
                                          const double suseptibility, const double msf,
                                          Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                          double II_F, Tensor<1, DIM> B)
  {
    Tensor<2, DIM> tmp;
    Tensor<1, DIM> FB = F*B;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
      {
        tmp[i][j] = mu*F[i][j] - mu*F_inv[j][i] +
                    (2.0*mu*nu/(1.0- nu))*(II_F*II_F - II_F)*F_inv[j][i] +
                    2.0*msf*suseptibility*FB[i]*B[j];
      }

    return tmp;

  }

  Tensor<4,DIM> Compressible_NH_Langevin::get_incremental_moduli_tensor(const double nu, const double mu,
                                                      const double suseptibility, const double msf,
                                                      Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                                      double II_F, Tensor<1, DIM> B)
  {
    Tensor<4,DIM> tmp;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int k=0; k<DIM; ++k)
          for (unsigned int l=0; l<DIM; ++l)
          {
            tmp[i][j][k][l] = ((i==k) && (j==l) ? mu : 0.0) +
                (mu - (2.0*mu*nu/(1.0-nu))*(II_F*II_F - II_F))*F_inv[j][k]*F_inv[l][i] +
                (4.0*nu*mu/(1.0-nu))*(II_F*II_F - 0.5*II_F)*F_inv[l][k]*F_inv[j][i] +
                (k == i ? 1.0 : 0.0)*2.0*msf*suseptibility*B[j]*B[l];
          }

    return tmp;
  }

  Tensor<3, DIM> Compressible_NH_Langevin::get_d2W_dFdB(const double nu, const double mu,
                                          const double suseptibility, const double msf,
                                          Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                          double II_F, Tensor<1, DIM> B)
  {
    Tensor<3, DIM> tmp;
    Tensor<1, DIM> FB = F*B;


    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int k=0; k<DIM; ++k)
        {
          tmp[i][j][k] = 2.0*msf*suseptibility*(F[i][k]*B[j] + (j == k ? 1.0 : 0.0)*FB[i]);
        }

    return tmp;
  }

  void Compressible_NH_Langevin::get_first_derivs(const bool rho, const double nu, const double mu,
      const double suseptibility, const double msf,
      Tensor<2,DIM> F, Tensor<1, DIM> B,  Tensor<1, DIM> *dW_dB, Tensor<2, DIM> *dW_dF)
  {

    Tensor<2, DIM> F_inv = invert(F);
    double II_F = determinant(F);
    double over_Jmu0 = 1.0/(II_F*MU0);
    Tensor<1, DIM> FB = F*B;
    double norm_sqr_FB = (FB)*(FB);


    (*dW_dB) = over_Jmu0*transpose(F)*FB;
    for(unsigned int i = 0; i < DIM; i ++)
      for(unsigned int j = 0; j < DIM; j ++)
      {
        (*dW_dF)[i][j] = over_Jmu0*(-0.5*norm_sqr_FB*F_inv[j][i] + FB[i]*B[j]);
      }

    if(rho == true)
    {
      // elastic part
      for (unsigned int i=0; i<DIM; ++i)
        for (unsigned int j=0; j<DIM; ++j)
        {
          (*dW_dF)[i][j] += mu*F[i][j] - mu*F_inv[j][i] +
                      (2.0*mu*nu/(1.0- nu))*(II_F*II_F - II_F)*F_inv[j][i];
        }

      if(fabs(suseptibility) >= 1e-6)
      {
       // magnetic part

        double alpha = MU0*msf/(3.0*suseptibility);
        double over_alpha = 1.0/alpha;
        double mag_FB = sqrt(norm_sqr_FB);
        double mag_b = mag_FB/II_F;

        double over_mag_FB_J = 1.0/(mag_FB*II_F);

        double eta = II_F*alpha*msf*(1.0/mag_b - over_alpha/tanh(over_alpha*mag_b));

        double psi = II_F*alpha*msf*(log(over_alpha*mag_b) - log(sinh(over_alpha*mag_b)));

        *dW_dB += eta*over_mag_FB_J*transpose(F)*FB;


        for (unsigned int i=0; i<DIM; ++i)
          for (unsigned int j=0; j<DIM; ++j)
          {
            (*dW_dF)[i][j] += eta*(-mag_b*F_inv[j][i] + over_mag_FB_J*FB[i]*B[j]) + psi*F_inv[j][i];
          }
      }
    }
  }

  void Compressible_NH_Langevin::get_second_derivs(const bool rho, const double nu, const double mu,
          const double suseptibility, const double msf,
          Tensor<2,DIM> F, Tensor<1, DIM> B, Tensor<2, DIM> *d2W_dBdB, Tensor<3, DIM> *d2W_dFdB, Tensor<4, DIM> *d2W_dFdF)
  {
    Tensor<2, DIM> F_inv = invert(F);
    double II_F = determinant(F);
    double over_Jmu0 = 1.0/(II_F*MU0);
    Tensor<1, DIM> FB = F*B;
    double norm_sqr_FB = (FB)*(FB);
    Tensor<2, DIM> C = transpose(F)*F;
    Tensor<1, DIM> CB = C*B;

    *d2W_dBdB = over_Jmu0*C;

    for(unsigned int i = 0; i < DIM; i ++)
      for(unsigned int j = 0; j < DIM; j ++)
        for(unsigned int k = 0; k < DIM; k ++)
        {
          (*d2W_dFdB)[i][j][k] = over_Jmu0*(F[i][k]*B[j] + (j == k ? 1.0 : 0.0)*FB[i] - CB[k]*F_inv[j][i]);
          for(unsigned int l = 0; l < DIM; l ++)
          {
            (*d2W_dFdF)[i][j][k][l] = over_Jmu0*( 0.5*norm_sqr_FB*(F_inv[l][k]*F_inv[j][i] + F_inv[j][k]*F_inv[l][i]) -
                                  FB[k]*B[l]*F_inv[j][i] - FB[i]*F_inv[l][k]*B[j] + (i ==k ? 1.0 : 0.0)*B[l]*B[j]);
          }
        }

    if(rho == true)
    {
      // elastic part
      for (unsigned int i=0; i<DIM; ++i)
        for (unsigned int j=0; j<DIM; ++j)
          for (unsigned int k=0; k<DIM; ++k)
            for (unsigned int l=0; l<DIM; ++l)
            {
              (*d2W_dFdF)[i][j][k][l]  += ((i==k) && (j==l) ? mu : 0.0) +
                  (mu - (2.0*mu*nu/(1.0-nu))*(II_F*II_F - II_F))*F_inv[j][k]*F_inv[l][i] +
                  (4.0*nu*mu/(1.0-nu))*(II_F*II_F - 0.5*II_F)*F_inv[l][k]*F_inv[j][i];

            }

      if(fabs(suseptibility) >= 1e-6)
      {
        // magnetic part

        double alpha = MU0*msf/(3.0*suseptibility);
        double over_alpha = 1.0/alpha;
        double mag_FB = sqrt(norm_sqr_FB);
        double mag_b = mag_FB/II_F;

        double over_mag_FB_J = 1.0/(mag_FB*II_F);
        double over_norm_sqr_FB = 1.0/norm_sqr_FB;
        double cotanh_alpha_mag_b = 1.0/tanh(over_alpha*mag_b);
        double eta = II_F*alpha*msf*(1.0/mag_b - over_alpha*cotanh_alpha_mag_b);
        double zeta = II_F*alpha*msf*(-1.0/(mag_b*mag_b) + over_alpha*over_alpha*(cotanh_alpha_mag_b*cotanh_alpha_mag_b - 1.0));
        double psi = II_F*alpha*msf*(log(over_alpha*mag_b) - log(sinh(over_alpha*mag_b)));

        Tensor<1, DIM> F_transpose_FB = transpose(F)*FB;

        for (unsigned int i=0; i<DIM; ++i)
          for (unsigned int j=0; j<DIM; ++j)
          {
            double dmag_b_dF_i_j = -mag_b*F_inv[j][i] + over_mag_FB_J*FB[i]*B[j];

            (*d2W_dBdB)[i][j] += over_mag_FB_J*(zeta*over_mag_FB_J*F_transpose_FB[i]*F_transpose_FB[j] +
                                  eta*(-over_norm_sqr_FB*F_transpose_FB[i]*F_transpose_FB[j] + C[i][j]));

            for(unsigned int k = 0; k < DIM; k ++)
            {
              (*d2W_dFdB)[i][j][k] += over_mag_FB_J*(zeta*F_transpose_FB[k]*dmag_b_dF_i_j +
                                       eta*(-over_norm_sqr_FB*F_transpose_FB[k]*FB[i]*B[j] +
                                             (F[i][k]*B[j] + (j == k ? 1.0 : 0.0)*FB[i])));
              for(unsigned int l = 0; l < DIM; l++)
              {
                double dpsi_dF_k_l = eta*(-mag_b*F_inv[l][k] + over_mag_FB_J*FB[k]*B[l]) + psi*F_inv[l][k];
                double dmag_b_dF_k_l = -mag_b*F_inv[l][k] + over_mag_FB_J*FB[k]*B[l];
                (*d2W_dFdF)[i][j][k][l] += zeta*dmag_b_dF_i_j*dmag_b_dF_k_l +
                                           eta*over_mag_FB_J*( -FB[k]*B[l]*F_inv[j][i] + norm_sqr_FB*F_inv[j][k]*F_inv[l][i]
                                                 -over_norm_sqr_FB*FB[k]*B[l]*FB[i]*B[j] + (i == k ? 1.0 : 0.0)*B[l]*B[j]) +
                                           dpsi_dF_k_l*F_inv[j][i] -  psi*F_inv[j][k]*F_inv[l][i];
              }
            }
          }
      }
    }

  }


  // The linear one

  inline
  double LinearLagrangian::get_energy(Tensor<4, DIM> D,
                              Tensor<2, DIM> E)
  {
    double W = 0.0;


    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int k=0; k<DIM; ++k)
          for (unsigned int l=0; l<DIM; ++l)
          {
            W += 0.5*D[i][j][k][l]*E[k][l]*E[i][j];
          }


    return W;
  }

  inline
  Tensor<2,DIM> LinearLagrangian::get_piola_kirchoff_tensor(Tensor<4, DIM> D,
                                            Tensor<2,DIM> F, Tensor<2,DIM> E)
  {

    Tensor<2, DIM> tmp;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int n=0; n<DIM; ++n)
          for (unsigned int l=0; l<DIM; ++l)
            for (unsigned int m=0; m<DIM; ++m)
            {
              tmp[i][j] += F[i][j]*D[j][l][m][n]*E[m][n];
            }



    return tmp;
  }

  inline
  Tensor<4,DIM> LinearLagrangian::get_incremental_moduli_tensor(Tensor<4, DIM> D, Tensor<2,DIM> F, Tensor<2,DIM> E)
  {
    Tensor<4, DIM> tmp;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int p=0; p<DIM; ++p)
          for (unsigned int q=0; q<DIM; ++q)
            for (unsigned int m=0; m<DIM; ++m)
              for (unsigned int n=0; n<DIM; ++n)
              {
                tmp[i][j][p][q] += (i == p ? 1.0:0.0)*D[j][q][m][n]*E[m][n] +
                                     F[i][m]*D[j][m][n][p]*F[n][q] + F[i][m]*D[j][m][p][n]*F[n][q];
              }




    return tmp;
  }

  inline
  Tensor<2, DIM> LinearLagrangian::get_lagrangian_strain(Tensor<2, DIM> F)
  {
    Tensor<2,DIM> E;
    for (unsigned int i=0; i<DIM; ++i)
    {
      for (unsigned int j=0; j<DIM; ++j)
      {
        for (unsigned int m=0; m<DIM; ++m)
        {
          E[i][j] += F[i][m]*F[j][m];
        }
      }
      E[i][i] += -1.0;
    }

    return E;
  }

  inline
  Tensor<4,DIM> LinearLagrangian::get_D(double mu, double nu)
  {
    Tensor<4, DIM> tmp;

    double scalingFactor = mu/(1.0 - nu*nu);
    double val = (1 - nu)/4.0;

    tmp[1][1][1][1] = 1.0;
    tmp[1][1][2][2] = nu;
    tmp[2][2][1][1] = nu;
    tmp[2][2][2][2] = 1.0;
    tmp[1][2][1][2] = val;
    tmp[1][2][2][1] = val;
    tmp[2][1][1][2] = val;
    tmp[2][1][2][1] = val;

    tmp *= scalingFactor;

    return tmp;
  }

#endif
