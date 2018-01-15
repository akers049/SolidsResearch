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

  };

  class Compressible_NH_Langevin
  {
    public:
    virtual ~Compressible_NH_Langevin(){};

    virtual double get_energy(const bool rho, const double nu, const double mu,
        const double suseptibility, const double msf,
        Tensor<2,DIM> F, Tensor<1, DIM> B);

    virtual Tensor<1,DIM> get_dW_dB(const double nu, const double mu,
                                    const double suseptibility, const double msf,
                                    Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                    double II_F, Tensor<1, DIM> B);

    virtual Tensor<2,DIM> get_d2W_dBdB(const double nu, const double mu,
                                    const double suseptibility, const double msf,
                                    Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                    double II_F, Tensor<1, DIM> B);

    virtual Tensor<2,DIM> get_piola_kirchoff_tensor(const double nu, const double mu,
                                            const double suseptibility, const double msf,
                                            Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                            double II_F, Tensor<1, DIM> B);

    virtual Tensor<4,DIM> get_incremental_moduli_tensor(const double nu, const double mu,
                                                        const double suseptibility, const double msf,
                                                        Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                                        double II_F, Tensor<1, DIM> B);

    virtual Tensor<3, DIM> get_d2W_dFdB(const double nu, const double mu,
                                            const double suseptibility, const double msf,
                                            Tensor<2,DIM> F, Tensor<2,DIM> F_inv,
                                            double II_F, Tensor<1, DIM> B);

    void get_first_derivs(const bool rho, const double nu, const double mu,
        const double suseptibility, const double msf,
        Tensor<2,DIM> F, Tensor<1, DIM> B,  Tensor<1, DIM> *dW_dB, Tensor<2, DIM> *dW_dF);

    void get_second_derivs(const bool rho, const double nu, const double mu,
        const double suseptibility, const double msf,
        Tensor<2,DIM> F, Tensor<1, DIM> B, Tensor<2, DIM> *d2W_dBdB, Tensor<3, DIM> *d2W_dFdB, Tensor<4, DIM> *d2W_dFdF);

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
