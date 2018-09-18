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
#include <deal.II/lac/arpack_solver.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/slepc_solver.h>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Constituitive.h"
#include "mapping_q_eulerian_hp.h"


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

     void set_nu_value(double nuVal){nuValue = nuVal;}
     double get_nu_value(){return nuValue;}

     double nuValue = 0.33;

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

  class FE_solution : public Function<DIM>
  {
  public:
    FE_solution(hp::DoFHandler<DIM, DIM>& dof_handler,
                Vector<double>& solution) : Function<DIM>(DIM)
    {
      dof_handler_ptr = &dof_handler;
      present_solution = &solution;
    }

    virtual ~FE_solution(){}

    virtual double value (const Point<DIM> &p,
                          const unsigned int  component) const;

    void reinit_solution(Vector<double>& solution)
    {
      present_solution = &solution;
    }

    hp::DoFHandler<DIM, DIM> *dof_handler_ptr;
    Vector<double> *present_solution;
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

    void update_bloch_wave_constraints(double wave_ratio);

    void update_F0(const double lambda);
    void add_small_pertubations(double amplitude, bool firstTime);

    void newton_iterate();

    void path_follow_PACA_iterate(Vector<double> *solVectorDir,
                                         double lambdaGuess, double ds);

    unsigned int get_system_eigenvalues(double lambda_eval, const int cycle);
    double get_bloch_eigenvalues(const int cycle, const int step, double wave_ratio, unsigned int indx, bool first);
    void set_unstable_eigenvector_as_initial_tangent(unsigned int indx);

    double bisect_find_lambda_critical(double lowerBound, double upperBound,
                                      double tol, unsigned int maxIter);
    void output_results(const unsigned int cycle);
    void output_load_info(std::vector<double> lambda_values,
                         std::vector<double> energy_values,
                         std::vector<double> congugate_lambda_values,
                         std::vector<double> displacement_magnitude,
                         const unsigned int cycle,
                         std::vector<bool> stable = std::vector<bool> (1, true)) const;

    void read_input_file(char* filename);
    void read_asymptotic_input_file(char* filename);

    void save_current_state(unsigned int indx);
    void load_state(unsigned int indx);
    void load_intial_tangent();
    void set_boundary_values();
    void print_dof_coords_and_vals(unsigned int indx);

    void interpolate_solution(FE_solution* solution_func_ptr)
    {
      VectorTools::interpolate(dof_handler, *solution_func_ptr, present_solution);
      constraints.distribute(present_solution);
    }

    // get methods for important constants
    double get_present_lambda(){return present_lambda;};
    void set_present_lambda(double lambda_val)
              { present_lambda = lambda_val;
                update_F0(present_lambda); };
    double get_ds(){return ds;};
    double get_first_bif_amplitude(){
      compute_first_bif_amplitude();
      return excee;};
    unsigned int get_output_every(){return output_every;};
    unsigned int get_load_steps(){return load_steps;};
    unsigned int get_n_dofs(){return dof_handler.n_dofs();};
    unsigned int get_number_active_cells(){return triangulation.n_active_cells();};
    unsigned int get_number_unit_cells(){return number_unit_cells;};
    FE_solution* get_fe_solution_function_ptr()
    {
      FE_solution_function = new FE_solution(dof_handler, present_solution);
      return FE_solution_function;
    }


    dealii::Vector<double>       present_solution;
    Vector<double>       evaluation_point;
    Vector<double>       initial_solution_tangent;
    double               initial_lambda_tangent = 0.0;


    double               system_energy = 0.0;
    double               congugate_lambda = 0.0;

    double u1u1 = 0.0;

    double critical_lambda_analytical = 0.0;
    double critical_frequency = 0.0;
    double w_c = 0.0;     // This is the same as critical_frequency...


  private:
    void init_fe_and_quad_collection();

    Tensor<2,DIM> get_deformation_gradient(std::vector<Tensor<1,DIM> > old_solution_gradient);

    void setup_system_constraints();
    void make_periodicity_constraints();
    void make_ave_x1_constraints();
    void make_symmetry_constraints();
    void setup_bloch_matched_dofs();
    void setup_bloch ();

    void assemble_system_matrix();
    void assemble_system_rhs();
    void assemble_drhs_dlambda();
    void assemble_bloch_matrix();

    void compute_first_bif_amplitude();
    void u1_value_list(const std::vector< Point< DIM > > &  points, std::vector<Vector<double>> & value_list);

    void apply_boundaries_and_constraints_system_matrix();
    void apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs);
    void apply_boundaries_and_constraints_bloch_matrix();


    void line_search_and_add_step_length(double current_residual, std::vector<bool> *homogenous_dirichlet_dofs);

    void line_search_and_add_step_length_PACA(double last_residual, std::vector<bool> *homogenous_dirichlet_dofs,
                                              Vector<double> *previousSolution, double previousLambda, double ds);
    void solve();
    void solve_boarder_matrix_system();

    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                            int const maxSize, int* const endOfFileFlag);

    void   map_dofs_to_support_points (const hp::DoFHandler<DIM, DIM>  &dof_handler, std::vector<Point<DIM> > &support_points);



    void renumber_boundary_ids();


    Triangulation<DIM,DIM>   triangulation;
    hp::DoFHandler<DIM, DIM>      dof_handler;
    hp::FECollection<DIM, DIM>    fe_collection;
    hp::QCollection<DIM>     quadrature_collection;

    unsigned int last_q_index = 0;

    ConstraintMatrix     constraints;
    ConstraintMatrix     constraints_bloch;
    ConstraintMatrix     bloch_hanging_node_constraints;


    std::vector<IndexSet>    owned_partitioning;
    std::vector<IndexSet>    relevant_partitioning;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    SparsityPattern      sparsity_pattern_bloch;
    SparseMatrix<double> bloch_matrix;

    PETScWrappers::SparseMatrix          system_matrix_petsc;
    PETScWrappers::SparseMatrix             stiffness_matrix;
    std::vector<PETScWrappers::MPI::Vector> eigenfunctions_;
    std::vector<double>                     eigenvalues_;

    Compressible_NeoHookean nh;

    double               present_lambda = 0.0;
    double               lambda_update= 0.0;

    Vector<double>       newton_update;
    Vector<double>       system_rhs;
    Vector<double>       drhs_dlambda;
    Vector<double>       solution_diff;
    double               lambda_diff = 0.0;
    double               rhs_bottom = 0.0;

    Tensor<2,DIM>        F0;

    std::vector<unsigned int>  grid_dimensions;
    bool nonUniform_mesh_flag = false;
    std::vector<double> domain_dimensions;
    double tol = 0.0;
    unsigned int maxIter = 0;
    double kappa = 0.0;

    double ds = 0.0;
    unsigned int load_steps = 0;
    unsigned int output_every = 0;
    unsigned int number_unit_cells = 1;
    unsigned int n_qpoints_x = 3;
    unsigned int n_qpoints_y = 3;
    unsigned int number_sections = 0;
    std::vector<unsigned int> section_polynomial_degrees;
    std::vector<unsigned int> section_FE_id;
    std::vector<unsigned int> FE_id_polynomial_degree;

    unsigned int max_degree = 0;

    bool fileLoadFlag = false;

    bool pieceConstFlag = false;
    std::vector<std::complex<double>> A;
    std::vector<std::complex<double>> B;
    std::vector<std::complex<double>> alphas;
    double L1 = 1.0;
    double excee = 0.0;


    char output_directory[MAXLINE];

    NuFunction nu;
    MuFunction *mu = NULL;
    FE_solution *FE_solution_function = NULL;

    bool bloch_matched_flag = false;
    std::vector<int> bloch_boundry_1_dof_index;
    std::vector<int> bloch_boundry_2_dof_index;
    std::vector<int> matched_dofs;

  };
}

#endif /* NEOHOOKEAN_NEWTON_COMPRESSEDSTRIP_H_ */
