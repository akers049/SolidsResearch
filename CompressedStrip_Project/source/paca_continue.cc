#include "CompressedStripPacaBloch.h"
#include "CompressedStripPacaBloch.cc"


using namespace dealii;
int main ()
{

  NeoHookean_Newton::ElasticProblem<2> ep;

  char fileName[MAXLINE];
  std::cout << "Please enter an input file: " << std::endl;
  std::cin >> fileName;
  ep.read_input_file(fileName);

  // read it the current state
  ep.load_state();

  ep.setup_system();

  std::cout << "\n   Number of active cells:       "
            << ep.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << ep.get_n_dofs()
            << std::endl << std::endl;

  Vector<double> previous_solution = ep.present_solution;

  double cows = ep.get_system_eigenvalues(ep.present_lambda, 123456);
  std::cout << "    Number negative Eigenvalues : " << cows << std::endl;

  ep.set_unstable_eigenvector(ep.present_lambda, 1);

  double lambda_tangent = 0.0;
  double scalingVal = sqrt(1 - lambda_tangent*lambda_tangent);
  ep.unstable_eigenvector *= scalingVal;

//  ep.present_solution = ep.unstable_eigenvector;
//  ep.output_results(123457);
//   exit(-1);

  ep.present_lambda -= 5e-4;
  ep.update_F0(ep.present_lambda);
  ep.newton_iterate();


  double previous_lambda = ep.present_lambda;
  previous_solution = ep.present_solution;

  ep.path_follow_PACA_iterate(ep.present_solution, ep.present_lambda,
                          ep.unstable_eigenvector, lambda_tangent, ep.get_ds());

  // trying something
  ep.unstable_eigenvector = ep.present_solution;
  ep.unstable_eigenvector -= previous_solution;
  ep.unstable_eigenvector *= -1.0;

  ep.path_follow_PACA_iterate(ep.present_solution, ep.present_lambda, ep.unstable_eigenvector, 0.0, 2.0*ep.get_ds());

  // define variables for the tangent to next start point.
  Vector<double> solution_tangent;

  std::vector<double> lambda_values;
  std::vector<double> congugate_lambda_values;
  std::vector<double> energy_values;
  std::vector<double> displacement_magnitude;

 // double num_negative_eigs = 0;
  for(unsigned int i = 1; i < ep.get_load_steps(); i ++)
  {

   // get the differences between past and this solution
   solution_tangent = ep.present_solution;
   solution_tangent -= previous_solution;
   lambda_tangent = ep.present_lambda - previous_lambda;

   // now scale to length 1.0
   scalingVal = 1.0/ep.get_ds();
   solution_tangent *= scalingVal;
   lambda_tangent *= scalingVal;

   previous_lambda = ep.present_lambda;
   previous_solution = ep.present_solution;

   ep.path_follow_PACA_iterate(ep.present_solution, ep.present_lambda, solution_tangent, lambda_tangent, ep.get_ds());
   std::cout << std::setprecision(15) << "    lambda = " << ep.present_lambda << std::endl;

   if ((i % ep.get_output_every()) == 0)
   {
     ep.output_results(i/ep.get_output_every() + 4000);
    // get_system_eigenvalues(present_lambda, i/output_every);
   }

//   num_negative_eigs = ep.get_system_eigenvalues(ep.present_lambda, i+1000);
//   std::cout << "    Number negative Eigenvalues : " << num_negative_eigs << std::endl;

   // get energy and congugate lambda value and save them. Make sure to scale by number
   // of unit cells...
   ep.assemble_system_energy_and_congugate_lambda(ep.present_lambda);
   lambda_values.push_back(ep.present_lambda);
   congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
   energy_values.push_back(ep.system_energy/ep.get_number_unit_cells());
   displacement_magnitude.push_back(ep.present_solution.l2_norm()/ep.get_number_unit_cells());


  }
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 3);

  return 0;
}
