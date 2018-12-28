#include "CompressedStripPacaBloch.h"
#include <fstream>

using namespace dealii;
int main (int argc, char** argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  compressed_strip::ElasticProblem ep_loaded;

  char fileName[MAXLINE];
  std::cout << "Please enter an input file for the loaded mesh solution: " << std::endl;
  std::cin >> fileName;
  unsigned int indx;
  std::cout << "Please enter an the index of the state: " << std::endl;
  std::cin >> indx;
  char fileName_new[MAXLINE];
  std::cout << "Please enter an input file for the new mesh: " << std::endl;
  std::cin >> fileName_new;


  ep_loaded.read_input_file(fileName);
  ep_loaded.load_state(indx);
  ep_loaded.setup_system();

  std::cout << "Inputted System : " << std::endl;
  std::cout << "\n   Number of active cells:       "
            << ep_loaded.get_number_active_cells()
            << std::endl;

  std::cout << "   Number of degrees of freedom: "
            << ep_loaded.get_n_dofs()
            << std::endl << std::endl;
  std::cout << "   Starting at lambda = " << ep_loaded.get_present_lambda() << std::endl;

  compressed_strip::ElasticProblem ep_new;
  ep_new.read_input_file(fileName_new);
  ep_new.create_mesh();
  ep_new.setup_system();

  std::cout << "New System : " << std::endl;
  std::cout << "\n   Number of active cells:       "
            << ep_new.get_number_active_cells()
            << std::endl;

  std::cout << "   Number of degrees of freedom: "
            << ep_new.get_n_dofs()
            << std::endl << std::endl;

  ep_new.set_present_lambda(ep_loaded.get_present_lambda());
  std::cout << "   Starting at lambda = " << ep_new.get_present_lambda() << std::endl;

  compressed_strip::FE_solution *fe_solution_ptr = ep_loaded.get_fe_solution_function_ptr();
  ep_new.interpolate_solution(fe_solution_ptr);

  compressed_strip::FE_solution *fe_tangent_ptr = ep_loaded.get_fe_tangent_ptr();
  ep_new.interpolate_tangent(fe_tangent_ptr);
  ep_new.initial_lambda_tangent = ep_loaded.initial_lambda_tangent;


  ep_new.newton_iterate();

  ep_loaded.output_results(100000);
  ep_new.output_results(100000);

  ep_new.save_current_state(indx);

  double error = ep_new.compute_difference(fe_solution_ptr, ep_loaded.get_max_degree());
  double sol1_norm = ep_loaded.compute_solution_l2_norm();
  double sol2_norm = ep_new.compute_solution_l2_norm();



  std::cout << "Loaded Solution Norm  : " << sol1_norm << std::endl;
  std::cout << "New Solution Norm     : " << sol2_norm << std::endl;

  std::cout << "\nAbsolute difference : " << error << std::endl;
  std::cout << "Relative Difference   : " << error/sol1_norm << std::endl;




  return 0;
}
