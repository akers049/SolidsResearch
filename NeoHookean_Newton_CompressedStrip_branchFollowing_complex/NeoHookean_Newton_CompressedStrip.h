/*
 * NeoHookean_Newton_CompressedStrip.h
 *
 *  Created on: Aug 10, 2017
 *      Author: andrew
 */

#ifndef NEOHOOKEAN_NEWTON_COMPRESSEDSTRIP_H_
#define NEOHOOKEAN_NEWTON_COMPRESSEDSTRIP_H_

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/lapack_full_matrix.h>

namespace NeoHookean_Newton
{
  using namespace dealii;
  /****************************************************************
                       Class Declarations
  ****************************************************************/

  /****  ElasticProblem  *****
   * This is the primary class used, with all the dealii stuff
   */
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem();
    ~ElasticProblem();
    void run ();

  private:
    Tensor<2,dim> get_deformation_gradient(std::vector<Tensor<1,dim> > old_solution_gradient);
    Tensor<4,dim> get_incremental_moduli_tensor(const double nu,
                                                const double mu, Tensor<2,dim> F_inv,
                                                double II_F);
    Tensor<2,dim> get_piola_kirchoff_tensor(const double nu, const double mu,
                                            Tensor<2,dim> F, Tensor<2,dim> F_inv,
                                            double II_F);

    double get_energy(const double nu, const double mu,
                      Tensor<2, dim> F, double II_F);


    void create_mesh();
    void setup_system ();
    void setup_constraints();
    void set_boundary_values();
    void update_F0(const double lambda);
    void add_small_pertubations(double amplitude, bool firstTime);
    void add_first_bif_displacements(double epsilon);
    void assemble_system_matrix();
    void assemble_system_rhs();
    void apply_boundaries_to_rhs(Vector<double> rhs, std::vector<bool> homogenous_dirichlet_dofs);
    void assemble_system_energy_and_congugate_lambda(double lambda);
    void newton_iterate();
    void line_search_and_add_step_length(double current_residual, std::vector<bool> homogenous_dirichlet_dofs);
    void solve();
    bool get_system_eigenvalues(double lambda_eval, const int cycle);
    double bisect_find_lambda_critical(double lowerBound, double upperBound,
                                       double tol, unsigned int maxIter);
    void output_results(const unsigned int cycle) const;
    void output_load_info(std::vector<double> lambda_values,
                          std::vector<double> energy_values,
                          std::vector<double> congugate_lambda_values) const;
    void read_input_file(char* filename);
    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                            int const maxSize, int* const endOfFileFlag);


    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    ConstraintMatrix     constraints;
    std::vector<IndexSet>    owned_partitioning;
    std::vector<IndexSet>    relevant_partitioning;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       present_solution;
    Vector<double>       evaluation_point;
    Vector<double>       newton_update;
    Vector<double>       system_rhs;

    Vector<double>       system_eigenvalues;

    double               system_energy = 0.0;
    double               congugate_lambda = 0.0;

    Tensor<2,dim>        F0;


    std::vector<unsigned int>  grid_dimensions;
    std::vector<double> domain_dimensions;
    double tol = 0.0;
    unsigned int maxIter = 0;
    double kappa = 0.0;
    double critical_lambda_analytical = 0.0;
    double critical_frequency = 0.0;
    std::vector<std::complex<double>> characteristic_roots;
    std::vector<std::complex<double>> amplitudes_v1;
    std::vector<std::complex<double>> amplitudes_v2;

    std::vector<int> matched_dofs;

  };


  /****  NuFunction *****
   * This is a dealii Function class for the nu function.
   */

  template <int dim>
  class NuFunction : public Function<dim>
  {

  public:
    NuFunction () : Function<dim>() {}
    virtual ~NuFunction (){}

    // const double PI = std::atan(1.0)*4;

    virtual double value (const Point<dim> &p,
                          const unsigned int  component = 0) const;

    virtual void value_list(const std::vector< Point< dim > > &  points,
                             std::vector< double > &   values,
                             const unsigned int  component = 0 )   const;

  };


  /****  NuFunction *****
   * This is a dealii Function class for the mu function.
   */

  template <int dim>
  class MuFunction : public Function<dim>
  {

  public:
    MuFunction (double expontential_growth_param) : Function<dim>()
    {
      kappa = expontential_growth_param;
    }
    virtual ~MuFunction (){}

    // const double PI = std::atan(1.0)*4;

    virtual double value (const Point<dim> &p,
                          const unsigned int  component = 0) const;

    virtual void value_list(const std::vector< Point< dim > > &  points,
                             std::vector< double > &   values,
                             const unsigned int  component = 0 )   const;

    double kappa;
    const double mu0 = 1.0;

  };

}


#endif /* NEOHOOKEAN_NEWTON_COMPRESSEDSTRIP_H_ */
