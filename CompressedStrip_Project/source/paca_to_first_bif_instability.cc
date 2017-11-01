#include "CompressedStripPacaBloch.h"
#include <fstream>

using namespace dealii;
int main ()
{

  compressed_strip::ElasticProblem ep;

  char fileName[MAXLINE];
  std::cout << "Please enter an input file: " << std::endl;
  std::cin >> fileName;
  ep.read_input_file(fileName);

  ep.create_mesh();

  ep.setup_system();

  std::cout << "\n   Number of active cells:       "
            << ep.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << ep.get_n_dofs()
            << std::endl << std::endl;


  // get the critical lambda value

  double lambda_c = ep.bisect_find_lambda_critical(ep.critical_lambda_analytical - 0.05,
                                               ep.critical_lambda_analytical + 0.05, 1e-6, 50);

  std::cout << "The lambda_c is: " << lambda_c << std::endl;

  ep.set_present_lambda(0.0);
  ep.newton_iterate();
  ep.output_results (0);
  ep.assemble_system_energy_and_congugate_lambda();


  double lambda_start = lambda_c + 1e-6;
  ep.set_present_lambda(lambda_start);
  ep.newton_iterate();
  ep.output_results (1);

  // set the eigenvector for the unstable mode
  ep.set_unstable_eigenvector_as_initial_tangent(1);
 // ep.initial_solution_tangent *= -1.0;

  Vector<double> previous_solution = ep.present_solution;
  double previous_lambda = lambda_start;


  ep.path_follow_PACA_iterate( &(ep.initial_solution_tangent), ep.initial_lambda_tangent, ep.get_ds());
  ep.output_results (2);

  // define variables for the tangent to next start point.
  Vector<double> solution_tangent;
  double lambda_tangent;
  double scalingVal;

  std::vector<double> lambda_values;
  std::vector<double> congugate_lambda_values;
  std::vector<double> energy_values;
  std::vector<double> displacement_magnitude;
  lambda_values.push_back(0.0);
  congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
  energy_values.push_back( ep.system_energy/ep.get_number_unit_cells());
  displacement_magnitude.push_back(ep.present_solution.l2_norm()/ep.get_number_unit_cells());


  unsigned int num_negative_eigs = 0;
  unsigned int prev_num_negative_eigs = 0;
  unsigned int step_number = 0;

  for(unsigned int i = 1; i < ep.get_load_steps(); i ++)
  {

    step_number ++;
   // get the differences between past and this solution
   solution_tangent = ep.present_solution;
   solution_tangent -= previous_solution;
   lambda_tangent = ep.get_present_lambda() - previous_lambda;

   // now scale to length 1.0
   scalingVal = 1.0/ep.get_ds();
   solution_tangent *= scalingVal;
   lambda_tangent *= scalingVal;

   previous_lambda = ep.get_present_lambda();
   previous_solution = ep.present_solution;

   ep.path_follow_PACA_iterate(&(solution_tangent), lambda_tangent, ep.get_ds());
   std::cout << std::setprecision(15) << "    lambda = " << ep.get_present_lambda() << std::endl;

   if ((i % ep.get_output_every()) == 0)
   {
      ep.output_results(i/ep.get_output_every() + 2);
   // ep.get_system_eigenvalues(ep.present_lambda, i/ep.get_output_every());
   }

   // get energy and congugate lambda value and save them.
   ep.assemble_system_energy_and_congugate_lambda();
   lambda_values.push_back(ep.get_present_lambda());
   congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
   energy_values.push_back(ep.system_energy/ep.get_number_unit_cells());
   displacement_magnitude.push_back(ep.present_solution.l2_norm()/(sqrt(1.0*ep.get_number_unit_cells())));

   prev_num_negative_eigs = num_negative_eigs;
   num_negative_eigs = ep.get_system_eigenvalues(ep.get_present_lambda(), i);
   std::cout << "    Number negative Eigenvalues : " << num_negative_eigs << std::endl;
   std::cout << "    Step Number : " << step_number << std::endl;
   if ((i > 40 && num_negative_eigs == (prev_num_negative_eigs + 1)))
   {
     std::cout << "\n Eigenvalue Crossing Found. Outputting current state and stopping" << std::endl;
     break;
   }
  }
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 1);

  ep.save_current_state(0);

  return 0;
}
