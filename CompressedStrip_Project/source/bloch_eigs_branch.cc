/*
 * bloch_eigs_first_bif.cc
 *
 *  Created on: Nov 1, 2017
 *      Author: andrew
 */

#include "CompressedStripPacaBloch.h"


using namespace dealii;
int main ()
{

  compressed_strip::ElasticProblem ep;

  char fileName[MAXLINE];
  std::cout << "Please enter an input file: " << std::endl;
  std::cin >> fileName;
  ep.read_input_file(fileName);

  unsigned int state_index;
  std::cout << "Please enter an the index of the state: " << std::endl;
  std::cin >> state_index;

  ep.load_state(state_index);

  ep.setup_system();

  std::cout << "\n   Number of active cells:       "
            << ep.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << ep.get_n_dofs()
            << std::endl << std::endl;


  Vector<double> previous_solution = ep.present_solution;

  unsigned int number_negative_eigs = ep.get_system_eigenvalues(ep.get_present_lambda(), -1);

  std::cout << "    Number negative Eigenvalues : " << number_negative_eigs << std::endl;

  ep.set_unstable_eigenvector_as_initial_tangent(number_negative_eigs);

  ep.initial_lambda_tangent = 0.5;
  double scalingVal = sqrt(1 - ep.initial_lambda_tangent*ep.initial_lambda_tangent);
  ep.initial_solution_tangent *= -1.0*scalingVal;

  double previous_lambda = ep.get_present_lambda();
  previous_solution = ep.present_solution;

  ep.path_follow_PACA_iterate(&(ep.initial_solution_tangent), ep.initial_lambda_tangent, ep.get_ds());

  // define variables for the tangent to next start point.
  Vector<double> solution_tangent;
  double lambda_tangent;

  std::vector<double> lambda_values(ep.get_load_steps());
  std::vector<double> congugate_lambda_values(ep.get_load_steps());
  std::vector<double> energy_values(ep.get_load_steps());
  std::vector<double> displacement_magnitude(ep.get_load_steps());
  lambda_values[0] = 0.0;
  congugate_lambda_values[0] = ep.congugate_lambda;
  energy_values[0] = ep.system_energy;


  for(unsigned int i = 1; i < ep.get_load_steps(); i ++)
  {
    std::cout << "    Step Number : " << i << std::endl;
    if(i%30 == 0)
    {
      for(unsigned int j = 0; j < 50; j++)
      {
        double wave_ratio = j*0.01;
        ep.get_bloch_eigenvalues(j, i/30, wave_ratio, 1);
      }
    }



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
    std::cout << std::setprecision(15) << "    Lambda = " << ep.get_present_lambda() << std::endl;


    // get energy and congugate lambda value and save them.
    ep.assemble_system_energy_and_congugate_lambda();
    lambda_values[i] = ep.get_present_lambda();
    congugate_lambda_values[i] = ep.congugate_lambda/ep.get_number_unit_cells();
    energy_values[i] = ep.system_energy/ep.get_number_unit_cells();
    displacement_magnitude[i] = ep.present_solution.l2_norm()/(sqrt(1.0*ep.get_number_unit_cells()));

  }
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 10);
}




