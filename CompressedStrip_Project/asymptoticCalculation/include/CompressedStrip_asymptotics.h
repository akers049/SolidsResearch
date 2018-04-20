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
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

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


#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <gsl/gsl_integration.h>

#include "Constituitive.h"


#define MAXLINE 1024
#define DIM 2

namespace compressed_strip
{
  using namespace dealii;
  /****************************************************************
                       Class Declarations
  ****************************************************************/

  /****  NuFunction *****
    * This is a dealii Function class for the nu function.
    */

   class NuFunction : public Function<DIM>
   {

   public:
     NuFunction () : Function<DIM>() {}
     virtual ~NuFunction (){}

     // const double PI = std::atan(1.0)*4;

     virtual double value (const Point<DIM> &p,
                           const unsigned int  component = 0) const;

     virtual void value_list(const std::vector< Point< DIM > > &  points,
                              std::vector< double > &   values,
                              const unsigned int  component = 0 )   const;

   };


   /****  NuFunction *****
    * This is a dealii Function class for the mu function.
    */

   class MuFunction : public Function<DIM>
   {

   public:
     MuFunction (double mu_1, double l1 = 0.0) : Function<DIM>()
     {
       kappa = mu_1;

       if ( l1 != 0.0)
       {
         Assert ((l1 > 0 && l1 < 1.0), ExcNotImplemented());
         pieceConstFlag = true;
         L1 = l1;
       }
       else
         pieceConstFlag = false;

     }
     virtual ~MuFunction (){}

     // const double PI = std::atan(1.0)*4;

     virtual double value (const Point<DIM> &p,
                           const unsigned int  component = 0) const;

     virtual void value_list(const std::vector< Point< DIM > > &  points,
                              std::vector< double > &   values,
                              const unsigned int  component = 0 )   const;

     bool pieceConstFlag;
     double L1;
     double kappa;
     const double mu0 = 1.0;

   };


  /****  ElasticProblem  *****
   * This is the primary class used, with all the dealii stuff
   */
  class ElasticProblem
  {
  public:
    ElasticProblem();
    ~ElasticProblem();

    void create_mesh();
    void setup_system ();

    void assemble_system_energy_and_congugate_lambda();
    void assemble_asymptotic_integrals();

    void update_F0(const double lambda);

    void read_input_file(char* filename);

    // get methods for important constants
    double get_present_lambda(){return present_lambda;};
    void set_present_lambda(double lambda_val)
              { present_lambda = lambda_val;
                update_F0(present_lambda); };
    unsigned int get_n_dofs(){return dof_handler.n_dofs();};
    unsigned int get_number_active_cells(){return triangulation.n_active_cells();};
    unsigned int get_number_unit_cells(){return number_unit_cells;};

    Vector<double>       present_solution;
    Vector<double>       evaluation_point;
    double               initial_lambda_tangent = 0.0;


    double               system_energy = 0.0;
    double               E_u1u1u1 = 0.0;
    double               E_u1u1u1u1 = 0.0;
    double               E_u2u1u1 = 0.0;
    double               dEdlambda_u1u1 = 0.0;
    double               congugate_lambda = 0.0;

    double critical_lambda_analytical = 0.0;
    double critical_frequency = 0.0;


  private:

    void get_grad_u1_value_list(const std::vector< Point< DIM > > &  points, std::vector<Tensor<2, DIM>> & value_list);


    void get_grad_u2_value_list(const std::vector< Point< DIM > > &  points, std::vector<Tensor<2, DIM>> & value_list);

    double E1(double x2);
    double E2(double x2);
    double E2_tilde(double x2);

    double E1exp_real(double x2)
    {
      unsigned int i = root_index;
      return( exp(-r_real[i]*x2)*cos(-r_imag[i]*x2)*E1(x2));
    };

    double E1exp_imag(double x2)
    {
      unsigned int i = root_index;
      return( exp(-r_real[i]*x2)*sin(-r_imag[i]*x2)*E1(x2));
    };

    double E2exp_real(double x2)
    {
      unsigned int i = root_index;
      return( exp(-r_real[i]*x2)*cos(-r_imag[i]*x2)*E2(x2));
    };
    double E2exp_imag(double x2)
    {
      unsigned int i = root_index;
      return( exp(-r_real[i]*x2)*sin(-r_imag[i]*x2)*E2(x2));
    };
    double L(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
    {
      return L_tensor[i-1][j-1][k-1][l-1];
    };
    double M(unsigned int i, unsigned int j, unsigned int k, unsigned int l,unsigned int m, unsigned int n)
    {
      return M_tensor[i-1][j-1][k-1][l-1][m-1][n-1];
    };

    Tensor<2,DIM> get_deformation_gradient(std::vector<Tensor<1,DIM> > old_solution_gradient);

    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                            int const maxSize, int* const endOfFileFlag);

    void renumber_boundary_ids();


    Triangulation<DIM,DIM>   triangulation;
    DoFHandler<DIM>      dof_handler;

    FESystem<DIM>        fe;

    std::vector<IndexSet>    owned_partitioning;
    std::vector<IndexSet>    relevant_partitioning;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Compressible_NeoHookean nh;

    double               present_lambda = 0.0;

    Tensor<2,DIM>        F0;

    std::vector<unsigned int>  grid_dimensions;
    bool nonUniform_mesh_flag = false;
    std::vector<double> domain_dimensions;
    double kappa = 0.0;

    double L1 = 1.0;

    unsigned int number_unit_cells = 1;

    std::vector<std::complex<double>> A;
    std::vector<std::complex<double>> B;
    std::vector<std::complex<double>> alphas;
    std::vector<std::complex<double>> r;

    std::vector<double> r_real;
    std::vector<double> r_imag;

    std::vector<std::complex<double>> phi1;
    std::vector<std::complex<double>> phi3;
    std::vector<std::complex<double>> C;
    std::vector<std::complex<double>> phi_inv_T_2;
    std::vector<std::complex<double>> phi_inv_T_4;

    double w_c = 0.0;

    Tensor<4, DIM>  L_tensor;
    Tensor<6, DIM>  M_tensor;

    bool pieceConstFlag = false;
    bool fileLoadFlag = false;

    char output_directory[MAXLINE];

    NuFunction nu;
    MuFunction *mu = NULL;


    unsigned int root_index = 0;

    // numerical integration stuff:
    unsigned long int np = 1000; // work area size
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (np);

  };
}

#endif /* NEOHOOKEAN_NEWTON_COMPRESSEDSTRIP_H_ */
