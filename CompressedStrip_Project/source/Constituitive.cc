#ifndef CONSTITUITIVE_CC_
#define CONSTITUITIVE_CC_

#include "Constituitive.h"

#define DIM 2

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

  inline
  Tensor<6, DIM> Compressible_NeoHookean::get_d3W_dFdFdF(const double nu, const double mu,
                                        Tensor<2,DIM> F_inv, double II_F)
  {
    Tensor<6, DIM> tmp;

    const double nu_frac = nu/(1.0 - nu);
    const double II_C = II_F*II_F;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int k=0; k<DIM; ++k)
          for (unsigned int l=0; l<DIM; ++l)
           for (unsigned int m=0; m<DIM; ++m)
             for (unsigned int n=0; n<DIM; ++n)
             {
               tmp[i][j][k][l][m][n] =
                   (2.0*nu_frac*(II_C - II_F) - 1.0)*(F_inv[j][m]*F_inv[n][k]*F_inv[l][i] +
                                                      F_inv[j][k]*F_inv[l][m]*F_inv[n][i]) +
                   -4.0*nu_frac*(II_C - 0.5*II_F)*(F_inv[n][m]*F_inv[j][k]*F_inv[l][i] +
                                                   F_inv[l][m]*F_inv[n][k]*F_inv[j][i] +
                                                   F_inv[l][k]*F_inv[j][m]*F_inv[n][i]) +
                    8.0*nu_frac*(II_C - 0.25*II_F)*F_inv[n][m]*F_inv[l][k]*F_inv[j][i];
             }

    tmp *= mu;

    return tmp;
  }

  inline
  Tensor<8, DIM> Compressible_NeoHookean::get_d4W_dFdFdFdF(const double nu, const double mu,
                                        Tensor<2,DIM> F_inv, double II_F)
  {
    Tensor<8, DIM> tmp;

    const double nu_frac = mu*nu/(1.0 - nu);
    const double II_C = II_F*II_F;

    for (unsigned int i=0; i<DIM; ++i)
      for (unsigned int j=0; j<DIM; ++j)
        for (unsigned int k=0; k<DIM; ++k)
          for (unsigned int l=0; l<DIM; ++l)
            for (unsigned int m=0; m<DIM; ++m)
              for (unsigned int n=0; n<DIM; ++n)
                for (unsigned int p=0; p<DIM; ++p)
                  for (unsigned int q=0; q<DIM; ++q)
                  {
                    tmp[i][j][k][l][m][n][p][q] =
                       (mu - 2.0*nu_frac*(II_C - II_F))*(F_inv[j][p]*F_inv[q][m]*F_inv[n][k]*F_inv[l][i] +
                                                          F_inv[j][m]*F_inv[n][p]*F_inv[q][k]*F_inv[l][i] +
                                                          F_inv[j][m]*F_inv[n][k]*F_inv[l][p]*F_inv[q][i] +
                                                          F_inv[j][p]*F_inv[q][k]*F_inv[l][m]*F_inv[n][i] +
                                                          F_inv[j][k]*F_inv[l][p]*F_inv[q][m]*F_inv[n][i] +
                                                          F_inv[j][k]*F_inv[l][m]*F_inv[n][p]*F_inv[q][i]) +
                       4.0*nu_frac*(II_C - 0.5*II_F)*(F_inv[j][m]*F_inv[n][k]*F_inv[l][i]*F_inv[q][p] +
                                                      F_inv[j][k]*F_inv[l][m]*F_inv[n][i]*F_inv[q][p] +
                                                      F_inv[n][p]*F_inv[q][m]*F_inv[j][k]*F_inv[l][i] +
                                                      F_inv[n][m]*F_inv[j][p]*F_inv[q][k]*F_inv[l][i] +
                                                      F_inv[n][m]*F_inv[j][k]*F_inv[l][p]*F_inv[q][i] +
                                                      F_inv[l][p]*F_inv[q][m]*F_inv[n][k]*F_inv[j][i] +
                                                      F_inv[l][m]*F_inv[n][p]*F_inv[q][k]*F_inv[j][i] +
                                                      F_inv[l][m]*F_inv[n][k]*F_inv[j][p]*F_inv[q][i] +
                                                      F_inv[l][p]*F_inv[q][k]*F_inv[j][m]*F_inv[n][i] +
                                                      F_inv[l][k]*F_inv[j][p]*F_inv[q][m]*F_inv[n][i] +
                                                      F_inv[l][k]*F_inv[j][m]*F_inv[n][p]*F_inv[q][i]) +
                       -8.0*nu_frac*(II_C - 0.25*II_F)*(F_inv[n][m]*F_inv[j][k]*F_inv[l][i]*F_inv[q][p] +
                                                        F_inv[l][m]*F_inv[n][k]*F_inv[j][i]*F_inv[q][p] +
                                                        F_inv[l][k]*F_inv[j][m]*F_inv[n][i]*F_inv[q][p] +
                                                        F_inv[n][p]*F_inv[q][m]*F_inv[l][k]*F_inv[j][i] +
                                                        F_inv[n][m]*F_inv[l][p]*F_inv[q][k]*F_inv[j][i] +
                                                        F_inv[n][m]*F_inv[l][k]*F_inv[j][p]*F_inv[q][i]) +
                       16.0*nu_frac*(II_C - 0.125*II_F)*F_inv[n][m]*F_inv[l][k]*F_inv[j][i]*F_inv[q][p];
                 }


    return tmp;
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
