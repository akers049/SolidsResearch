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
#define MU_VALUE 1.0
#define NU_VALUE 0.33

namespace NeoHookean_Newton
{
  using namespace dealii;

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
                                  const double mu,
                                  std::vector<Tensor<1,dim> > old_solution_gradient);
    Tensor<2,dim> get_piola_kirchoff_tensor(const double nu, const double mu,
            std::vector<Tensor<1,dim> > old_solution_gradient);


    void create_mesh();
    void setup_system (const bool initial_step);
    void set_boundary_values();
    void update_F0(const double lambda);
    void add_small_pertubations(double amplitude);
    void assemble_system_matrix();
    void assemble_system_rhs();
    void newton_iterate(const double tolerance,
                        const unsigned int max_iteration);
    void line_search_and_add_step_length(double current_residual);
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;
    void print_average_deformation_gradient();


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

    Tensor<2,dim>        F0;

  };

  template <int dim>
  class NuFunction : public Function<dim>
  {

  public:
    NuFunction () : Function<dim>() {}
    ~NuFunction (){}

    const double PI = std::atan(1.0)*4;

    virtual double value (const Point<dim> &p,
                          const unsigned int  component = 0) const;

    virtual void value_list(const std::vector< Point< dim > > &  points,
                             std::vector< double > &   values,
                             const unsigned int  component = 0 )   const;

  };

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues (const double beta) : Function<dim>(dim)
    {
      deformationTensor[0][0] = 1 + beta*0.0;
      deformationTensor[0][1] = 0 + beta*1.5;
      deformationTensor[1][0] = 0 + beta*0.0;
      deformationTensor[1][1] = 1 + beta*0.0;

    }

    ~BoundaryValues(){};

    Tensor<2,dim> deformationTensor;

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
  };

  template <int dim>
  double NuFunction<dim>::value (const Point<dim>  &p,
                                 const unsigned int  component) const
  {


    double nuValue = 0;
    if (component == 0)
    {
      //nuValue = NU_VALUE + 0.25*cos(p(0)*2.0*PI/(LENGTH/4.0));
      nuValue = NU_VALUE + 0.2*(p(1) - 0.5);
    }

    return nuValue;
  }

  template <int dim>
  void NuFunction<dim>::value_list(const std::vector< Point< dim > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = NuFunction<dim>::value(points[i], component);

  }

  template <int dim>
  double
  BoundaryValues<dim>::value (const Point<dim>  &p,
                              const unsigned int component) const
  {

    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    double boundaryValue = 0.0;


    if(component == 0)
    {
      boundaryValue = deformationTensor[0][0]*p(0) + deformationTensor[0][1]*p(1) - p(0);
    }
    else if(component == 1)
    {
      boundaryValue = deformationTensor[1][0]*p(0) + deformationTensor[1][1]*p(1) - p(1);
    }

    return boundaryValue;
  }


  template <int dim>
  void
  BoundaryValues<dim>::vector_value (const Point<dim> &p,
                                     Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = BoundaryValues<dim>::value (p, c);
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
  inline
  Tensor<2,dim>
  ElasticProblem<dim>::get_deformation_gradient(std::vector<Tensor<1,dim> > old_solution_gradient)
  {
    Tensor<2,dim> tmp = F0;

    for (unsigned int i = 0; i < dim; i ++)
      for(unsigned int j = 0; j < dim; j++)
      {
        tmp[i][j] += old_solution_gradient[i][j];
      }



    return tmp;
  }


  template <int dim>
  inline
  Tensor<4,dim>
  ElasticProblem<dim>::get_incremental_moduli_tensor(const double nu,
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
  ElasticProblem<dim>::get_piola_kirchoff_tensor(const double nu, const double mu,
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

    Point<dim> corner1, corner2;
    corner1(0) = -LENGTH/2.0;
    corner1(1) =  0.0;
    corner2(0) =  LENGTH/2.0;
    corner2(1) =  HEIGHT;
    GridGenerator::hyper_rectangle (triangulation, corner1, corner2, true);

	  triangulation.refine_global(5);
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

      output_results (0);


      constraints.clear ();

      DoFTools::make_hanging_node_constraints (dof_handler, constraints);

      DoFTools::make_periodicity_constraints<DoFHandler<dim>>(dof_handler, 0, 1, 0, constraints);

      // now do constraints that the average x1 displacement on boundary 2 is zero

      std::vector<bool> x1_components = {true, false};
      ComponentMask x1_mask(x1_components);

      std::vector<bool> boundary_2_dof_0 (dof_handler.n_dofs(), false);
      std::vector<bool> boundary_0_dof_0 (dof_handler.n_dofs(), false);
      std::vector<bool> boundary_1_dof_0 (dof_handler.n_dofs(), false);

      std::set< types::boundary_id > boundary_id_2;
      boundary_id_2.insert(2);

      std::set< types::boundary_id > boundary_id_0;
      boundary_id_0.insert(0);

      std::set< types::boundary_id > boundary_id_1;
      boundary_id_1.insert(1);


      DoFTools::extract_boundary_dofs (dof_handler,
                                       x1_mask,
                                       boundary_2_dof_0,
                                       boundary_id_2);

      DoFTools::extract_boundary_dofs (dof_handler,
                                       x1_mask,
                                       boundary_0_dof_0,
                                       boundary_id_0);

      DoFTools::extract_boundary_dofs (dof_handler,
                                       x1_mask,
                                       boundary_1_dof_0,
                                       boundary_id_1);


      unsigned int first_boundary_dof = 0;
      for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
      {
        if ((boundary_2_dof_0[i] == true) && (boundary_1_dof_0[i] == false) && (boundary_0_dof_0[i] == false))
        {
           first_boundary_dof = i;
          break;
        }
      }
      constraints.add_line (first_boundary_dof);
      for (unsigned int i=(first_boundary_dof + 1); i<dof_handler.n_dofs(); ++i)
      {
        if ((boundary_2_dof_0[i] == true) && (boundary_1_dof_0[i] == false) && (boundary_0_dof_0[i] == false))
          constraints.add_entry (first_boundary_dof, i, -1);

      }

      constraints.close ();

	  }

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);

  }

  template<int dim>
  void ElasticProblem<dim>::set_boundary_values()
  {
    // this sets the boundary values of the solution vector so that the Newton step
    // can use homogeneous direchlet conditions on the set boundaries.

    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 2,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 0,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);
    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 1,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        present_solution,
                                        system_rhs);
    for (std::map<types::global_dof_index, double>::const_iterator
         p = boundary_values.begin();
         p != boundary_values.end(); ++p)
      present_solution(p->first) = p->second;

  }

  template<int dim>
  void ElasticProblem<dim>::update_F0(const double lambda)
  {
    Assert ((lambda >= 0 && lambda <= 1.0), ExcInternalError());

    double lambda1 = 1 - lambda;
    double lambda2;

    // solve quadratic for lambda2
    double a = (1.0 + (2.0*NU_VALUE/(1.0 - NU_VALUE))*lambda1*lambda1);
    double b = (-2.0*NU_VALUE/(1.0 - NU_VALUE))*lambda1;
    double c = -1;

    lambda2 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
    if(lambda2 < 1.0)
      lambda2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);


    std::cout << "Lambda 1 : " << lambda1 << std::endl;
    std::cout << "Lambda 2 : " << lambda2 << std::endl;

    F0[0][0] = lambda1;
    // F0[1][1] = lambda1/2;
    F0[1][1] = lambda2;
    F0[0][1] = 0.0;
    F0[1][0] = 0.0;
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

    ConstantFunction<dim>  mu(MU_VALUE);
    NuFunction<dim> nu;

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

    constraints.condense (system_matrix);

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side2_components = {false, true};
    ComponentMask side2_mask(side2_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
                                              ZeroFunction<dim>(dim),
                                              boundary_values,
                                              side2_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);

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

    ConstantFunction<dim> mu(MU_VALUE);
    NuFunction<dim> nu;

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

    constraints.condense (system_rhs);

    std::map<types::global_dof_index,double> boundary_values;

    std::vector<bool> side2_components = {false, true};
    ComponentMask side2_mask(side2_components);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
                                              ZeroFunction<dim>(dim),
                                              boundary_values,
                                              side2_mask);

    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);
  }


  template <int dim>
  void ElasticProblem<dim>::newton_iterate(const double tolerance,
                                           const unsigned int max_iteration)
  {

    double current_residual;
    unsigned int iteration = 0;


    evaluation_point = present_solution;
    assemble_system_rhs();
    current_residual = system_rhs.l2_norm();

    while((current_residual > tolerance) && (iteration <
        max_iteration))
    {
      assemble_system_matrix();
      solve();
      line_search_and_add_step_length(current_residual);

      evaluation_point = present_solution;
      assemble_system_rhs();
      current_residual = system_rhs.l2_norm();
      iteration ++;
    }

    std::cout << "    Converging Iterations : " << iteration << "\n";

  }

  template<int dim>
  void ElasticProblem<dim>::line_search_and_add_step_length(double last_residual)
  {
    double current_residual;
    for(double alpha = 1.0; alpha > 1e-5; alpha *=0.5)
    {
      evaluation_point = present_solution;
      evaluation_point.add(alpha, newton_update);
      assemble_system_rhs();
      current_residual = system_rhs.l2_norm();

      if(current_residual < last_residual)
        break;

    }
    present_solution = evaluation_point;
  }

  template <int dim>
  void ElasticProblem<dim>::solve ()
  {
    SolverControl           solver_control (10000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, newton_update, system_rhs,
              preconditioner);

    constraints.distribute (newton_update);

  }



  template <int dim>
  void ElasticProblem<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(2),
                                        typename FunctionMap<dim>::type(),
                                        present_solution,
                                        estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);

    triangulation.execute_coarsening_and_refinement ();
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
        break;
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
  void ElasticProblem<dim>::print_average_deformation_gradient()
  {
    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<std::vector<Tensor<1,dim> > > old_solution_gradients(n_q_points, std::vector<Tensor<1,dim>>(dim));
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();

    Tensor<2,dim, double> aveDeformationGradient;
    double count = 0;
    for (; cell!=endc; ++cell)
    {
      count = count + 1.0;;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(present_solution, old_solution_gradients);

      Tensor<2, dim, double> F  = get_deformation_gradient(old_solution_gradients[0]);
      aveDeformationGradient += F;
    }

    aveDeformationGradient = aveDeformationGradient/count;

    count = 0;
    Tensor<2,dim, double> deformationGradientSTD;
    for (; cell!=endc; ++cell)
    {
      count = count + 1.0;

      fe_values.reinit (cell);

      fe_values.get_function_gradients(present_solution, old_solution_gradients);

      Tensor<2, dim , double> F  = get_deformation_gradient(old_solution_gradients[0]);
      for (int i = 0; i <dim; i++)
        for(int j = 0; j < dim; j++)
          deformationGradientSTD[i][j] += (F[i][j] - aveDeformationGradient[i][j])*(F[i][j] - aveDeformationGradient[i][j]);
    }

    for (int i = 0; i <dim; i++)
     for(int j = 0; j < dim; j++)
       deformationGradientSTD[i][j] = sqrt(deformationGradientSTD[i][j]/(count - 1.0));


    std::cout << "Average Deformation Gradient:\n";
    for (int i = 0; i <dim; i++)
    {
      for(int j = 0; j < dim; j++)
      {
        std::cout << std::setprecision(15) << aveDeformationGradient[i][j] << "  ";
      }
      std::cout << "\n";
    }

    std::cout << "\nDeformation Gradient std:\n";
    for (int i = 0; i <dim; i++)
    {
      for(int j = 0; j < dim; j++)
      {
        std::cout << std::setprecision(15) << deformationGradientSTD[i][j] << "  ";
      }
      std::cout << "\n";
    }

  }


  template <int dim>
  void ElasticProblem<dim>::run ()
  {


        	create_mesh();

          setup_system(true);




      setup_system(true);

      std::cout << "   Number of active cells:       "
                << triangulation.n_active_cells()
                << std::endl;


      std::cout << "   Number of degrees of freedom: "
                << dof_handler.n_dofs()
                << std::endl;


      double tol = 1e-10;



        double lambda = 0.01;

        add_small_pertubations(0.01);
        update_F0(lambda);
        set_boundary_values();

        newton_iterate(tol, 50);


        output_results (1);
        print_average_deformation_gradient();
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
