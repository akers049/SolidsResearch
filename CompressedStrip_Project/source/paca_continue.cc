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
  double lambda_start = 0.364731946335432 ;
  double previous_lambda = lambda_start;

  double lambda_tangent = 0.0;
  double scalingVal = sqrt(1 - lambda_tangent*lambda_tangent);
  ep.unstable_eigenvector *= scalingVal;

  ep.path_follow_PACA_iterate(ep.present_solution, lambda_start,
                          ep.unstable_eigenvector, 0.0, ep.get_ds());


  // define variables for the tangent to next start point.
  Vector<double> solution_tangent;

  std::vector<double> lambda_values;
  std::vector<double> congugate_lambda_values;
  std::vector<double> energy_values;
  std::vector<double> displacement_magnitude;

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
     ep.output_results(i/ep.get_output_every() + 1000);
    // get_system_eigenvalues(present_lambda, i/output_every);
   }

   // get energy and congugate lambda value and save them.
   ep.assemble_system_energy_and_congugate_lambda(ep.present_lambda);
   lambda_values.push_back(ep.present_lambda);
   congugate_lambda_values.push_back(ep.congugate_lambda);
   energy_values.push_back(ep.system_energy);
   displacement_magnitude.push_back(ep.present_solution.l2_norm());


  }
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 5);

  return 0;
}
