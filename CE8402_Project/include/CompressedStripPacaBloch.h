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
     MuFunction () : Function<DIM>() {}

     virtual ~MuFunction (){}

     // const double PI = std::atan(1.0)*4;

     virtual double value (const Point<DIM> &p,
                           const unsigned int  component = 0) const;

     virtual void value_list(const std::vector< Point< DIM > > &  points,
                              std::vector< double > &   values,
                              const unsigned int  component = 0 )   const;

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

    void newton_iterate();

    unsigned int get_system_eigenvalues(const int cycle);

    void output_results(const unsigned int cycle) const;
    void output_load_info(std::vector<double> lambda_values,
                         std::vector<double> energy_values,
                         std::vector<double> congugate_lambda_values,
                         std::vector<double> displacement_magnitude,
                         const unsigned int cycle)  const;

    void read_input_file(char* filename);

    void save_current_state(unsigned int indx);
    void load_state(unsigned int indx);

    void set_boundary_values();
    void print_dof_coords_and_vals(unsigned int indx);

    // get methods for important constants

    double get_ds(){return ds;};
    unsigned int get_output_every(){return output_every;};
    double       get_final_load(){return final_load;};
    unsigned int get_load_steps(){return load_steps;};
    unsigned int get_n_dofs(){return dof_handler.n_dofs();};
    unsigned int get_number_active_cells(){return triangulation.n_active_cells();};

    double get_load_val(){ return load_val;};
    void   set_load_val(double val){ load_val = val;};
    void   increment_load_val(double inc){load_val += inc;};

    Vector<double>       present_solution;
    Vector<double>       evaluation_point;
    Vector<double>       initial_solution_tangent;
    double               initial_lambda_tangent = 0.0;


    double               system_energy = 0.0;
    double               congugate_lambda = 0.0;




  private:
    void right_hand_side (const std::vector<Point<DIM> > &points,
                          std::vector<Tensor<1, DIM> >   &values);

    Tensor<4, DIM> get_the_D(double mu, double nu);
    Tensor<2,DIM> get_deformation_gradient(std::vector<Tensor<1,DIM> > old_solution_gradient);
    Tensor<2, DIM> get_lagrangian_strain(Tensor<2, DIM> F);

    void setup_system_constraints();

    void assemble_system_matrix();
    void assemble_system_rhs();
    void assemble_drhs_dlambda();

    void apply_boundaries_and_constraints_system_matrix();
    void apply_boundaries_to_rhs(Vector<double> *rhs, std::vector<bool> *homogenous_dirichlet_dofs);
    void apply_loads_to_rhs();
    void apply_loads_to_system_matrix();


    void line_search_and_add_step_length(double current_residual, std::vector<bool> *homogenous_dirichlet_dofs);

    void solve();

    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                            int const maxSize, int* const endOfFileFlag);

    void renumber_boundary_ids();

    Triangulation<DIM,DIM>   triangulation;
    DoFHandler<DIM>      dof_handler;

    FESystem<DIM>        fe;

    ConstraintMatrix     constraints;

    std::vector<IndexSet>    owned_partitioning;
    std::vector<IndexSet>    relevant_partitioning;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Compressible_NeoHookean nh;


    Vector<double>       newton_update;
    Vector<double>       system_rhs;
    Vector<double>       drhs_dlambda;
    Vector<double>       solution_diff;

    Tensor<4, DIM> D;

    std::vector<unsigned int>  grid_dimensions;
    std::vector<double> domain_dimensions;
    double tol = 0.0;
    unsigned int maxIter = 0;

    double ds = 0.0;
    double final_load = 0.0;
    unsigned int load_steps = 0;
    unsigned int output_every = 0;

    double load_val = 0.0;

    bool problem_2_flag = false;

    char output_directory[MAXLINE];

    NuFunction nu;
    MuFunction mu;

    std::vector<unsigned int> corner_dofs;
    std::vector<unsigned int> inside_dofs;

  };
}

#endif /* NEOHOOKEAN_NEWTON_COMPRESSEDSTRIP_H_ */
