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



#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/solver_cg.h>
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

#include <deal.II/lac/solver_minres.h>

#include <fstream>
#include <iostream>

#define LENGTH 3.0
#define HEIGHT 1.0
#define LAMBDA 0.0013
#define MU_VALUE 1.0
#define NU_VALUE 0.33

namespace NeoHookean_Newton
{
  using namespace dealii;

    template <int dim>
    inline
    Tensor<2,dim>
    get_deformation_gradient(std::vector<Tensor<1,dim> > old_solution_gradient)
	  {
      Tensor<2,dim> tmp;

      for (unsigned int i = 0; i < dim; i ++)
        for(unsigned int j = 0; j < dim; j++)
        {
        	tmp[i][j] = old_solution_gradient[i][j];
        	if(i == j)
        	  tmp[i][j] += 1.0;
        }

      return tmp;
	  }


    template <int dim>
    inline
    Tensor<4,dim>
    get_incremental_moduli_tensor(const double nu,
    		                      const double mu,
    		                      std::vector<Tensor<1,dim> > old_solution_gradient)
    {
    	Tensor<4,dim> tmp;

    	Tensor<2, dim> F  = get_deformation_gradient(old_solution_gradient);
    	double II_F = determinant(F);
    	Tensor<2, dim> F_inv = invert(F);
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            for (unsigned int l=0; l<dim; ++l)
            {
              tmp[i][j][k][l] = ((i==k) && (j==l) ? mu : 0.0) +
                  (mu - (2.0*mu*nu/(1.0-nu))*(II_F*II_F - II_F))*F_inv[j][k]*F_inv[l][i] +
                  (4.0*nu*mu/(1.0-nu))*(II_F*II_F - 0.5*II_F)*F_inv[l][k]*F_inv[j][i];
            }

      return tmp;
    }

    template <int dim>
    inline
    Tensor<2,dim>
    get_piola_kirchoff_tensor(const double nu, const double mu,
        std::vector<Tensor<1,dim> > old_solution_gradient)
    {
      Tensor<2, dim> tmp;
      Tensor<2, dim> F  = get_deformation_gradient(old_solution_gradient);
      double II_F = determinant(F);
      Tensor<2, dim> F_inv = invert(F);
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
        {
          tmp[i][j] = mu*F[i][j] - mu*F_inv[j][i] +
                      (2.0*mu*nu/(1.0- nu))*(II_F*II_F - II_F)*F_inv[j][i];
        }

      return tmp;

    }

    template <int dim>
    inline
    double get_energy(const double nu, const double mu,
        std::vector<Tensor<1,dim> > old_solution_gradient)
    {
      Tensor<2, dim> F  = get_deformation_gradient(old_solution_gradient);
      double II_F = determinant(F);
      double I_C = F[1][0]*F[1][0] + F[0][0]*F[0][0] + F[0][1]*F[0][1] + F[1][1]*F[1][1];
      double W = mu*(0.5*(I_C - 2 - log(II_F*II_F)) + (nu/(1- nu))*(II_F - 1)*(II_F - 1));

      return W;
    }

    template <int dim>
    void right_hand_side (const std::vector<Point<dim> > &points,
                          std::vector<Tensor<1, dim> >   &values)
    {
      Assert (values.size() == points.size(),
              ExcDimensionMismatch (values.size(), points.size()));
      Assert (dim >= 2, ExcNotImplemented());


      // not imposing body forces or tractions
      for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
      {
        values[point_n][0] = 0.0;
        values[point_n][1] = 0.0;
      }

    }

  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem();
    ~ElasticProblem();
    void run ();

  private:
    void create_mesh();
    void setup_system (const bool initial_step);
    void add_small_pertubations(double amplitude);
    void assemble_system_matrix();
    void assemble_system_rhs();
    double assemble_system_energy();
    void compute_and_compare_first_numer_deriv(double epsilon);
    void compute_and_compare_second_numer_deriv(double epsilon);


    void output_results(const unsigned int cycle) const;


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

  };



  template <int dim>
  ElasticProblem<dim>::ElasticProblem ()
    :
    dof_handler (triangulation),
    fe (FE_Q<dim>(1), dim)
  {}




  template <int dim>
  ElasticProblem<dim>::~ElasticProblem ()
  {
    dof_handler.clear ();
  }


  template <int dim>
  void ElasticProblem<dim>::create_mesh()
  {

	  GridGenerator::hyper_cube (triangulation, 0.0, 1.0, false);

	  triangulation.refine_global(1);
  }

  template <int dim>
  void ElasticProblem<dim>::setup_system (const bool initial_step)
  {

	  if (initial_step)
	  {
      dof_handler.distribute_dofs (fe);
      present_solution.reinit (dof_handler.n_dofs());

      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

      constraints.clear ();

      DoFTools::make_hanging_node_constraints (dof_handler,
                                                   constraints);
      constraints.close ();
	  }

	  evaluation_point = present_solution;
    system_rhs.reinit (dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

  }

  template <int dim>
  void ElasticProblem<dim>::add_small_pertubations(double amplitude)
  {
    std::srand(5); //some seed

    for(unsigned int i = 0; i < dof_handler.n_dofs(); i++)
    {
      present_solution[i] += (2.0*(std::rand()/double(RAND_MAX)) - 1.0)*amplitude;
    }

    evaluation_point = present_solution;

  }

  template <int dim>
  void ElasticProblem<dim>::assemble_system_matrix()
  {
    system_matrix = 0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    ConstantFunction<dim> nu(NU_VALUE), mu(MU_VALUE);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<4,dim> d2W_dFdF = get_incremental_moduli_tensor(nu_values[q_point],
                                           mu_values[q_point], old_solution_gradients[q_point]);

        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          for (unsigned int m = 0; m < dofs_per_cell; ++m)
          {
            const unsigned int component_m = fe.system_to_component_index(m).first;

            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                for (unsigned int k=0; k<dim; ++k)
                  for (unsigned int l=0; l<dim; ++l)
                  {
                    cell_matrix(n,m) += (d2W_dFdF[i][j][k][l]*
                        ((i==component_n) && (k == component_m) ?
                            fe_values.shape_grad(n, q_point)[j]*fe_values.shape_grad(m, q_point)[l] : 0.0))*fe_values.JxW(q_point);
                  }
          }
        }
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
            for (unsigned int m=0; m<dofs_per_cell; ++m)
            {
              system_matrix.add (local_dof_indices[n],
                                 local_dof_indices[m],
                                 cell_matrix(n,m));
            }
    }

  }

  template <int dim>
  void ElasticProblem<dim>::assemble_system_rhs()
  {

    system_rhs = 0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    ConstantFunction<dim> nu(NU_VALUE), mu(MU_VALUE);

    std::vector<Tensor<1, dim> > rhs_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);
      right_hand_side (fe_values.get_quadrature_points(), rhs_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        Tensor<2,dim> dW_dF = get_piola_kirchoff_tensor(nu_values[q_point], mu_values[q_point], old_solution_gradients[q_point]);
        for (unsigned int n = 0; n < dofs_per_cell; ++n)
        {
          const unsigned int component_n = fe.system_to_component_index(n).first;

          for(unsigned int i = 0; i<dim; ++i)
            for(unsigned int j = 0; j<dim; ++j)
            {
              cell_rhs(n) -= dW_dF[i][j]*(i==component_n ? fe_values.shape_grad(n, q_point)[j] : 0.0)*fe_values.JxW(q_point);
            }

          cell_rhs(n) += fe_values.shape_value(n, q_point)*rhs_values[q_point][component_n]*fe_values.JxW(q_point);
        }
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int n=0; n<dofs_per_cell; ++n)
        system_rhs(local_dof_indices[n]) += cell_rhs(n);

      }
  }


  template <int dim>
  double ElasticProblem<dim>::assemble_system_energy()
  {

    double system_energy = 0;

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,dim>>(dim));

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>     nu_values (n_q_points);
    std::vector<double>     mu_values (n_q_points);

    ConstantFunction<dim> nu(NU_VALUE), mu(MU_VALUE);

    std::vector<Tensor<1, dim> > rhs_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      cell_rhs = 0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(evaluation_point, old_solution_gradients);

      nu.value_list (fe_values.get_quadrature_points(), nu_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);


      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        double W = get_energy(nu_values[q_point], mu_values[q_point], old_solution_gradients[q_point]);
        system_energy += W*fe_values.JxW(q_point);
      }
    }

    return system_energy;
  }

  template <int dim>
  void ElasticProblem<dim>::compute_and_compare_first_numer_deriv(double epsilon)
  {

    const int numberDofs = (const int) dof_handler.n_dofs();

    evaluation_point = present_solution;
    assemble_system_rhs();
    system_rhs *= -1.0;
    Vector<double> residual = system_rhs;

    double energy = assemble_system_energy();

    Vector<double> numericalResidual(dof_handler.n_dofs());
    numericalResidual = 0;

    for(int i = 0; i < numberDofs; i++)
    {
      evaluation_point = present_solution;
      evaluation_point[i] += epsilon;
      double perturbedEnergy = assemble_system_energy();
      numericalResidual[i] = (perturbedEnergy - energy)/epsilon;
    }

    std::ofstream out("firstDeriv.dat");
    out << ""  << std::endl;
    out << "Residual Value              Difference from numerical"  << std::endl;
    Vector<double> residualDifference(dof_handler.n_dofs());
    for(unsigned int i = 0;  i < dof_handler.n_dofs(); i ++)
      residualDifference[i] = residual[i] - numericalResidual[i];

    for(unsigned int i = 0; i < dof_handler.n_dofs(); i++)
      out << residual[i] << "     " << residualDifference[i] << std::endl;

    out << "\n\nL2 norm of the difference :" << residualDifference.l2_norm() << std::endl;
    out.close();

  }

  template <int dim>
  void ElasticProblem<dim>::compute_and_compare_second_numer_deriv(double epsilon)
  {
    int numberDofs = (int) dof_handler.n_dofs();

    evaluation_point = present_solution;
    assemble_system_rhs();
    assemble_system_matrix();
    system_rhs *= -1.0;
    Vector<double> residual = system_rhs;
    std::vector<std::vector<double>> numericalStiffness(numberDofs);
    for(int i = 0; i < numberDofs; i++)
      numericalStiffness[i].resize(numberDofs);

    for(int j = 0; j < numberDofs; j++)
    {
      evaluation_point = present_solution;
      evaluation_point[j] += epsilon;
      assemble_system_rhs();
      system_rhs *= -1.0;
      for(int i = 0; i < numberDofs; i++)
      {
        numericalStiffness[i][j] = (system_rhs[i] - residual[i])/epsilon;
      }
    }

    std::ofstream out("secondDeriv.dat");
    out << "" << std::endl;
    out << "Stiffness Matrix" << std::endl;


    for(int i = 0; i < numberDofs; i++)
    {
      for(int j = 0; j < numberDofs; j++)
      {
        out << system_matrix.el(i,j) << " ";
      }
      out << std::endl;
    }

    out << "Numerical" << std::endl;

    double diffNorm = 0;
    for(int i = 0; i < numberDofs; i++)
    {
      for(int j = 0; j < numberDofs; j++)
      {
        double diff;
        diff = (system_matrix.el(i,j) - numericalStiffness[i][j]);
        diffNorm += diff*diff;
        out << diff << " ";
      }
      out << std::endl;
    }
    diffNorm = sqrt(diffNorm);

    out << std::endl;
    out << "Norm of the difference: " << diffNorm << std::endl;

    out.close();

  }


  template <int dim>
  void ElasticProblem<dim>::output_results (const unsigned int cycle) const
  {
    std::string filename = "solution-";
    filename += ('0' + cycle);
    Assert (cycle < 100, ExcInternalError());

    filename += ".vtk";
    std::ofstream output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);



    std::vector<std::string> solution_names;
    switch (dim)
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
      }

    data_out.add_data_vector (present_solution, solution_names);
    data_out.build_patches ();
    data_out.write_vtk (output);

    // Now output the deformed mesh

    // this stuff was written by Krishanu Sen

     DataOut<dim> deformed_data_out;
     std::vector<DataComponentInterpretation::DataComponentInterpretation>
     data_component_interpretation(dim,DataComponentInterpretation::component_is_part_of_vector);

     std::vector<std::string> solution_name(dim, "displacement");

     deformed_data_out.attach_dof_handler(dof_handler);

     deformed_data_out.add_data_vector(present_solution,
            solution_name,
            DataOut<dim>::type_dof_data,
            data_component_interpretation);

     Vector<double> soln(present_solution.size());
     for (unsigned int i = 0; i < soln.size(); ++i)
       soln(i) = present_solution(i);

     MappingQEulerian<dim> q_mapping(1,  dof_handler, soln);  //Notice that I am hardcoding the degree
     deformed_data_out.build_patches(q_mapping, 1);      // Degree helps to create a refined mesh on the output, but with interpolated values, i.e. the result will not be better than what it is computed by deal ii (it is just post-processing)
    ///////////////////////////////////////////////////////////////////

     std::string deformed_filename = "deformed_mesh-";
     deformed_filename += ('0' + cycle);
     deformed_filename += ".vtk";
     std::ofstream output_deformed_mesh(deformed_filename.c_str());
     deformed_data_out.write_vtk(output_deformed_mesh);


  }



  template <int dim>
  void ElasticProblem<dim>::run ()
  {

      create_mesh();

      setup_system(true);


      add_small_pertubations(0.1724);

      std::cout << "   Number of active cells:       "
                << triangulation.n_active_cells()
                << std::endl;


      std::cout << "   Number of degrees of freedom: "
                << dof_handler.n_dofs()
                << std::endl;

      double energy = assemble_system_energy();

      std::cout << "System Energy: " << energy << std::endl;

      double epsilon = 1e-6;
      compute_and_compare_first_numer_deriv(epsilon);
      compute_and_compare_second_numer_deriv(epsilon);



  }
}


int main ()
{
  try
    {
      NeoHookean_Newton::ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
