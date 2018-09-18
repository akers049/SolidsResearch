/*
 * bloch_eigs_first_bif.cc
 *
 *  Created on: Sep 27, 2017
 *      Author: andrew
 */

#include "CompressedStripPacaBloch.h"
#include <time.h>
#include <sys/time.h>

using namespace dealii;

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

int main (int argc, char** argv)
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


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

  std::cout << std::setprecision (16)  << "The lambda_c is: " << lambda_c << std::endl;

  double wall0 = get_wall_time();

  ep.update_F0(0);
  ep.newton_iterate();
  ep.output_results (0);
  ep.assemble_system_energy_and_congugate_lambda();

  std::vector<double> lambda_values;
  std::vector<double> congugate_lambda_values;
  std::vector<double> energy_values;
  std::vector<double> displacement_magnitude;
  std::vector<bool>   stability;

  for(unsigned int i = 0; i < 100; i ++)
  {
    ep.set_present_lambda((i/50.0)*lambda_c);
    lambda_values.push_back(ep.get_present_lambda());
    ep.assemble_system_energy_and_congugate_lambda();
    congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
    energy_values.push_back(ep.system_energy/ep.get_number_unit_cells());
    displacement_magnitude.push_back(0.0);

    if (i > 50)
      stability.push_back(true);
    else
      stability.push_back(false);
  }

  double lambda_start = lambda_c + 1e-5;
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

  unsigned int num_blochs = 40;

  bool current_stability = false;
  for(unsigned int i = 1; i < ep.get_load_steps(); i ++)
  {
    std::cout << "    Step Number : " << i << std::endl;
    if(i%ep.get_output_every() == 0 )
    {
      current_stability = true;
      for(unsigned int j = 0; j < num_blochs; j++)
      {
        double wave_ratio = j*(0.5/num_blochs);
        double lowestEig = ep.get_bloch_eigenvalues(j, i/ep.get_output_every(), wave_ratio, 0, (j == 0 ? true : false));

        if(lowestEig < -1e-6)
          current_stability = false;

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
    lambda_values.push_back(ep.get_present_lambda());
    congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
    energy_values.push_back(ep.system_energy/ep.get_number_unit_cells());
    stability.push_back(current_stability);
    std::cout << "    Stable : " << current_stability << std::endl;
    displacement_magnitude.push_back(ep.present_solution.l2_norm()/(sqrt(1.0*ep.get_number_unit_cells())));

//    displacement_magnitude.push_back(ep.get_first_bif_amplitude());

  }

  double wall1 = get_wall_time();
  std::cout << "\n\n\nWall Time = " << wall1 - wall0 << std::endl;
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 0, stability);
}




