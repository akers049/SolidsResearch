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

#define MU_VALUE 1.0
#define NU_VALUE 0.33






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
                                  const double mu,
                                  std::vector<Tensor<1,dim> > old_solution_gradient);
    Tensor<2,dim> get_piola_kirchoff_tensor(const double nu, const double mu,
            std::vector<Tensor<1,dim> > old_solution_gradient);


    void create_mesh(std::vector<double> domain_dimensions, std::vector<unsigned int> grid_dimensions);
    void setup_system ();
    void set_boundary_values();
    void update_F0(const double lambda);
    void add_small_pertubations(double amplitude);
    void assemble_system_matrix();
    void assemble_system_rhs();
    void newton_iterate(const double tolerance,
                        const unsigned int max_iteration);
    void line_search_and_add_step_length(double current_residual);
    void solve();
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

    Tensor<2,dim>        F0;

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

    const double PI = std::atan(1.0)*4;

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
    MuFunction () : Function<dim>() {}
    virtual ~MuFunction (){}

    const double PI = std::atan(1.0)*4;

    virtual double value (const Point<dim> &p,
                          const unsigned int  component = 0) const;

    virtual void value_list(const std::vector< Point< dim > > &  points,
                             std::vector< double > &   values,
                             const unsigned int  component = 0 )   const;

  };

  /****  UnfiromDeformation *****
   * This is a class used ONLY by the output results to trasform
   * coordinates for outputing the deformed mesh...
   */
  template <int dim>
  class UniformDeformation
  {

  public:
    UniformDeformation(Tensor<dim, dim> F_init) : F0(F_init){}
    Point<dim> operator() (const Point<dim> &p) const { return Point<2>(F0[0][0]*p(0), F0[1][1]*p(1));}

  private:
    Tensor<dim, dim>   F0;

  };


  /****************************************************************
                       Function Definitions
  ****************************************************************/

  template <int dim>
  double NuFunction<dim>::value (const Point<dim>  &p,
                                 const unsigned int  component) const
  {
    double nuValue = NU_VALUE + 0.2*(p(1) - 0.5);

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
  double MuFunction<dim>::value (const Point<dim>  &p,
                                 const unsigned int  component) const
  {
    double muValue = MU_VALUE;

    return muValue;
  }

  template <int dim>
  void MuFunction<dim>::value_list(const std::vector< Point< dim > > &  points,
                           std::vector< double > &   values,
                           const unsigned int  component ) const
  {
    for(unsigned int i = 0; i < points.size(); i++)
      values[i] = MuFunction<dim>::value(points[i], component);
  }



  // computes right hand side values if we were to have body forces. But it just
  // always returns zeros because we don't.
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

  // Functions that get the tensors used in the stiffness and rhs calculations.

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
  void ElasticProblem<dim>::create_mesh(std::vector<double> domain_dimensions, std::vector<unsigned int> grid_dimensions)
  {
    // creates our strip.

    Point<dim> corner1, corner2;
    corner1(0) = -domain_dimensions[0]/2.0;
    corner1(1) =  0.0;
    corner2(0) =  domain_dimensions[0]/2.0;
    corner2(1) =  domain_dimensions[1];
    GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, true);
  }

  template <int dim>
  void ElasticProblem<dim>::setup_system ()
  {
    // Sets up system. Makes the constraint matrix, and reinitializes the
    // vectors and matricies used throughout to be the proper size.

    dof_handler.distribute_dofs (fe);
    present_solution.reinit (dof_handler.n_dofs());

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

    constraints.clear ();

    DoFTools::make_hanging_node_constraints (dof_handler, constraints);

    DoFTools::make_periodicity_constraints<DoFHandler<dim>>(dof_handler, 0, 1, 0, constraints);

    // now do constraints that the average x1 displacement on boundary 2 is zero

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::vector<bool> boundary_dof_x1 (dof_handler.n_dofs(), false);
    std::vector<bool> boundary_0_dof_x1 (dof_handler.n_dofs(), false);
    std::vector<bool> boundary_1_dof_x1 (dof_handler.n_dofs(), false);


    std::set< types::boundary_id > boundary_id_0;
    boundary_id_0.insert(0);

    std::set< types::boundary_id > boundary_id_1;
    boundary_id_1.insert(1);


    DoFTools::extract_boundary_dofs(dof_handler,
                                    x1_mask,
                                    boundary_dof_x1);

    DoFTools::extract_boundary_dofs (dof_handler,
                                     x1_mask,
                                     boundary_0_dof_x1,
                                     boundary_id_0);

    DoFTools::extract_boundary_dofs (dof_handler,
                                     x1_mask,
                                     boundary_1_dof_x1,
                                     boundary_id_1);


    unsigned int first_boundary_dof = 0;
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if ((boundary_dof_x1[i] == true) && (boundary_0_dof_x1[i] == false) && (boundary_1_dof_x1[i] == false))
      {
        first_boundary_dof = i;
        break;
      }
    }
    constraints.add_line (first_boundary_dof);
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      if (i == first_boundary_dof)
        continue;

      if(boundary_dof_x1[i] == true)
        constraints.add_entry (first_boundary_dof, i, -1);

    }

    constraints.close ();


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
    // can use homogeneous direchlet conditions on the set boundaries. It also makes sure
    // that the periodic faces' DoFs start with the same values (zero).

    std::vector<bool> x1_components = {true, false};
    ComponentMask x1_mask(x1_components);

    std::map<types::global_dof_index,double> boundary_values;

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 0,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 1,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 2,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values);

    VectorTools::interpolate_boundary_values (   dof_handler,
                                                 3,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values,
                                                 x1_mask);

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
    // Update the uniform diagonal F0 tensor.

    Assert ((lambda >= 0 && lambda <= 1.0), ExcInternalError());

    double lambda1 = 1 - lambda;
    double lambda2;

    // solve quadratic for lambda2
    double a = (1.0 + (2.0*NU_VALUE/(1.0 - NU_VALUE))*lambda1*lambda1);
    double b = (-2.0*NU_VALUE/(1.0 - NU_VALUE))*lambda1;
    double c = -1;


    // make sure we get the root that is greater than 1.0
    lambda2 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
    if(lambda2 < 1.0)
      lambda2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);


    std::cout << "Lambda 1 : " << lambda1 << std::endl;
    std::cout << "Lambda 2 : " << lambda2 << std::endl;

    F0[0][0] = lambda1;
    F0[1][1] = lambda2;
    F0[0][1] = 0.0;
    F0[1][0] = 0.0;
  }


  template <int dim>
  void ElasticProblem<dim>::add_small_pertubations(double amplitude)
  {
    // This is just to make the solution vector non-zero so it doesn't start at the right answer.

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
    // Assembling the system matrix. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

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
    // Assembling the system rhs. I choose to make the rhs and system matrix assemblies separate,
    // because we only do one at a time anyways in the newton method.

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
    /* This function executes the newton iterations until it converges to a solution
     * or it exceeds the maximum number of iterations.
     */

    double current_residual;
    unsigned int iteration = 0;


    // Starts by getting the residual in the current configuration.
    evaluation_point = present_solution;
    assemble_system_rhs();
    current_residual = system_rhs.l2_norm();

    // Loops until coverge or go over max iterations
    while((current_residual > tolerance) &&
             (iteration < max_iteration))
    {
      // Assemble the stiffness matrix
      assemble_system_matrix();

      // solve for the newton step
      solve();

      // Find the step length and add it to the current solution.
      // This function also calls assemble_system_rhs() so we don't need to
      // do another rhs call.
      line_search_and_add_step_length(current_residual);

      evaluation_point = present_solution;
      current_residual = system_rhs.l2_norm();
      iteration ++;
    }

    // output iterations for convergance.
    std::cout << "    Converging Iterations : " << iteration << "\n";

  }

  template<int dim>
  void ElasticProblem<dim>::line_search_and_add_step_length(double last_residual)
  {
   /* this function makes sure that the step sizes we are taking with
    * the newton iteration are making the residual get smaller.
    * Something very similar is used in the dealii step 57?. for now
    * it doesn't do much but we might need to give it more capabilites in the future.
    */

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
  void ElasticProblem<dim>::output_results (const unsigned int cycle) const
  {

    // So this first part is pretty much taken directly from one of the dealii
    // steps to output the data in vtk format
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

    // Need to make a copy of the current mesh, add the uniform deformation,
    // then add the displacements.

    Triangulation<dim> triangulation_deformed;
    triangulation_deformed.copy_triangulation(triangulation);

    DoFHandler<dim> dof_handler_deformed(triangulation_deformed);

    UniformDeformation<dim> uniform_deform(F0);

    GridTools::transform(uniform_deform, triangulation_deformed);

    dof_handler_deformed.distribute_dofs(fe);

    DataOut<dim> deformed_data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim,DataComponentInterpretation::component_is_part_of_vector);

    std::vector<std::string> solution_name(dim, "displacement");
    deformed_data_out.attach_dof_handler(dof_handler_deformed);

    deformed_data_out.add_data_vector(present_solution,
                                      solution_name,
                                      DataOut<dim>::type_dof_data,
                                      data_component_interpretation);

    MappingQEulerian<dim> q_mapping(1,  dof_handler_deformed, present_solution);
    deformed_data_out.build_patches(q_mapping, 1);

    std::string deformed_filename = "deformed_mesh-";
    deformed_filename += ('0' + cycle);
    deformed_filename += ".vtk";
    std::ofstream output_deformed_mesh(deformed_filename.c_str());
    deformed_data_out.write_vtk(output_deformed_mesh);


  }

  template <int dim>
  void ElasticProblem<dim>::run ()
  {


    // These are the lengths of the domain in the x1 and x2 directions.
    double L = 3.0;
    double H = 1.0;

    // These are the dimensions of the grid in the x1 and x2 directions.
    unsigned int x1_grid_dimensions = 4;
    unsigned int x2_grid_dimensions = 3;

    double final_lambda = 0.3;
    int load_steps = 8;

    double tol = 1e-10;


    std::vector<double> domain_dimensions(dim);
    domain_dimensions[0] = L;
    domain_dimensions[1] = H;

    std::vector<unsigned int>  grid_dimensions(dim);
    grid_dimensions[0] = x1_grid_dimensions;
    grid_dimensions[1] = x2_grid_dimensions;

    create_mesh(domain_dimensions, grid_dimensions);

    setup_system();

    std::cout << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;


    std::cout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;




    // small pertubations are added to the non-constrained DoFs. This is so that in the isotropic case
    // when the expected solution vector is zero, it doesn't just start at the solution and has to iterate...
    add_small_pertubations(0.01);

    set_boundary_values();


    double lambda_step = final_lambda/load_steps;
    double lambda = 0.0;
    for(int i = 0; i < load_steps; i++)
    {
      lambda +=lambda_step;
      update_F0(lambda);
      newton_iterate(tol, 50);
      output_results (i);
    }

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


