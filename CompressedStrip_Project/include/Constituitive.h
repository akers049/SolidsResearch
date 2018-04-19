#ifndef CONSTITUITIVE_H_
#define CONSTITUITIVE_H_

#include <iostream>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>

#define DIM 2

using namespace dealii;

  class Compressible_NeoHookean
  {
    public:
    virtual ~Compressible_NeoHookean(){};

    virtual double get_energy(const double nu, const double mu,
                              Tensor<2, DIM> F, double II_F);
    virtual Tensor<2,DIM> get_piola_kirchoff_tensor(const double nu, const double mu,
                                            Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                            double II_F);
    virtual Tensor<4,DIM> get_incremental_moduli_tensor(const double nu,
                                                const double mu, Tensor<2,DIM> F_inv,
                                                double II_F);
    virtual Tensor<6, DIM> get_d3W_dFdFdF(const double nu, const double mu,
                                          Tensor<2,DIM> F_inv, double II_F);
    virtual Tensor<8, DIM> get_d4W_dFdFdFdF(const double nu, const double mu,
                                            Tensor<2,DIM> F_inv, double II_F);
  };

  class LinearLagrangian
  {
    public:
    virtual ~LinearLagrangian(){};


    double get_energy(Tensor<4, DIM> D, Tensor<2, DIM> E);

    Tensor<2,DIM> get_piola_kirchoff_tensor(Tensor<4, DIM> D,
                                            Tensor<2,DIM> F, Tensor<2,DIM> E);

    Tensor<4,DIM> get_incremental_moduli_tensor(Tensor<4, DIM> D, Tensor<2,DIM> F,
                                              Tensor<2,DIM> E);

    Tensor<2, DIM> get_lagrangian_strain(Tensor<2, DIM> F);

    Tensor<4,DIM> get_D(double mu, double nu);

  };
#endif
