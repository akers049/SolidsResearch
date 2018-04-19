/*
 * bloch_eigs_first_bif.cc
 *
 *  Created on: Sep 27, 2017
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

  ep.create_mesh();

  ep.setup_system();

  std::cout << "\n   Number of active cells:       "
            << ep.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << ep.get_n_dofs()
            << std::endl << std::endl;


  // get the critical lambda value
  ep.evaluation_point = ep.present_solution;

  double lambda_c = ep.bisect_find_lambda_critical(ep.critical_lambda_analytical - 0.05,
                                               ep.critical_lambda_analytical + 0.05, 1e-6, 50);

  std::cout << "The lambda_c is: " << lambda_c << std::endl;

  ep.update_F0(0);
  ep.newton_iterate();
  ep.output_results (0);
  ep.assemble_system_energy_and_congugate_lambda();


  double lambda_start = lambda_c + 1e-6;
  ep.update_F0(lambda_start);
  ep.newton_iterate();
  ep.output_results (1);

  // set the eigenvector for the unstable mode
  ep.set_present_lambda(lambda_start);
  ep.set_unstable_eigenvector_as_initial_tangent(1);

  Vector<double> previous_solution = ep.present_solution;
  double previous_lambda = lambda_start;

  ep.set_present_lambda(lambda_start - 1e-6);

  ep.path_follow_PACA_iterate(&(ep.initial_solution_tangent), ep.initial_lambda_tangent, ep.get_ds());
  ep.output_results (2);

  // define variables for the tangent to next start point.
  Vector<double> solution_tangent;
  double lambda_tangent;
  double scalingVal;

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
        ep.get_bloch_eigenvalues(j, i/30, wave_ratio, 0);
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
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 0);
}




