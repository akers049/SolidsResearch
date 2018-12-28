#include "CompressedStripPacaBloch.h"
#include <fstream>

using namespace dealii;
int main (int argc, char** argv)
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  compressed_strip::ElasticProblem ep;

  char fileName[MAXLINE];
  std::cout << "Please enter an input file: " << std::endl;
  std::cin >> fileName;

  unsigned int indx;
  std::cout << "Please enter an the index of the state: " << std::endl;
  std::cin >> indx;

  ep.read_input_file(fileName);

  // read in the current state
  ep.load_state(indx);

  ep.setup_system();

  std::cout << "\n   Number of active cells:       "
            << ep.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << ep.get_n_dofs()
            << std::endl << std::endl;
  std::cout << "   Starting at lambda = " << ep.get_present_lambda() << std::endl;



  Vector<double> previous_solution = ep.present_solution;

  unsigned int number_negative_eigs = ep.get_system_eigenvalues(ep.get_present_lambda(), -1);
  std::cout << "    Number negative Eigenvalues : " << number_negative_eigs << std::endl;

  double scalingVal;
  if(ep.initial_lambda_tangent == 0.0)
  {
    ep.initial_lambda_tangent = 0.3;
    scalingVal = sqrt(1 - ep.initial_lambda_tangent*ep.initial_lambda_tangent);
    ep.initial_solution_tangent *= scalingVal;
  }

  double previous_lambda = ep.get_present_lambda();
  previous_solution = ep.present_solution;
  ep.path_follow_PACA_iterate((ep.initial_solution_tangent), ep.initial_lambda_tangent, ep.get_ds());
  std::cout << std::setprecision(15) << "    lambda = " << ep.get_present_lambda() << std::endl;

//   trying something
//  if (ep.get_number_unit_cells()%2 == 0)
//  {
//    ep.initial_solution_tangent = ep.present_solution;
//    ep.initial_solution_tangent -= previous_solution;
//    ep.initial_solution_tangent *= -1.0;
//
//    ep.path_follow_PACA_iterate(&(ep.initial_solution_tangent), 0.0, 2.0*ep.get_ds());
//  }

  // define variables for the tangent to next start point.
  Vector<double> solution_tangent;

  std::vector<double> lambda_values;
  std::vector<double> congugate_lambda_values;
  std::vector<double> energy_values;
  std::vector<double> displacement_magnitude;

  double lambda_tangent = 0.0;

  unsigned int num_negative_eigs = 0;
  unsigned int prev_num_negative_eigs = 0;
  unsigned int step_number = 0;

  bool eig_crossing = false;
  for(unsigned int i = 1; i < ep.get_load_steps(); i ++)
  {
    step_number++;

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

   ep.path_follow_PACA_iterate((solution_tangent), lambda_tangent, ep.get_ds());
   std::cout << std::setprecision(15) << "    lambda = " << ep.get_present_lambda() << std::endl;
   std::cout << "    Step Number: " << i << std::endl;

   if ((i % ep.get_output_every()) == 0)
   {
     ep.output_results(i/ep.get_output_every() + indx*1000);
   }


   // get energy and congugate lambda value and save them. Make sure to scale by number
   // of unit cells...
   ep.assemble_system_energy_and_congugate_lambda();
   lambda_values.push_back(ep.get_present_lambda());
   congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
   energy_values.push_back(ep.system_energy/ep.get_number_unit_cells());
   displacement_magnitude.push_back(ep.present_solution.l2_norm()/(sqrt(1.0*ep.get_number_unit_cells())));

   prev_num_negative_eigs = num_negative_eigs;

//   num_negative_eigs = ep.get_system_eigenvalues(ep.get_present_lambda(), i);
   std::cout << "    Number negative Eigenvalues : " << num_negative_eigs << std::endl;
   if ((i > 10 && num_negative_eigs == (prev_num_negative_eigs + 1)))
   {
     std::cout << "\n Eigenvalue Crossing Found. Outputting current state and stopping" << std::endl;
     eig_crossing = true;
     break;
   }


  }
  if(eig_crossing)
  {
    double lambda_c = ep.bisect_find_lambda_critical(previous_lambda,
                                              ep.get_present_lambda(), 1e-6, 50, true);
    double lambda_start = lambda_c + 1e-6;
    ep.set_present_lambda(lambda_start);
    ep.newton_iterate();
    std::cout << std::setprecision (16) << "The lambda_c is: " << lambda_c << std::endl;

    unsigned int number_negative_eigs = ep.get_system_eigenvalues(ep.get_present_lambda(), -1);
    ep.set_unstable_eigenvector_as_initial_tangent(number_negative_eigs);
    ep.initial_lambda_tangent = 0.0;
  }
  else
  {
    ep.initial_solution_tangent = ep.present_solution;
    ep.initial_solution_tangent -= previous_solution;
    ep.initial_lambda_tangent = ep.get_present_lambda() - previous_lambda;
  }

  ep.save_current_state(indx+1);
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, indx+1);

  return 0;
}
