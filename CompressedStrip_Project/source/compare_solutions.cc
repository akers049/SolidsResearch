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
  ep_loaded.read_input_file(fileName);

  // read it the current state
  ep_loaded.load_state(1);

  ep_loaded.setup_system();

  std::cout << "\n   Number of active cells:       "
            << ep_loaded.get_number_active_cells()
            << std::endl;


  std::cout << "   Number of degrees of freedom: "
            << ep_loaded.get_n_dofs()
            << std::endl << std::endl;
  std::cout << "   Starting at lambda = " << ep_loaded.get_present_lambda() << std::endl;

  compressed_strip::ElasticProblem ep_new;
  std::cout << "Please enter an input file for the new mesh: " << std::endl;
  std::cin >> fileName;
  ep_new.read_input_file(fileName);

  ep_new.create_mesh();

  ep_new.setup_system();


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

  ep_new.newton_iterate();



  return 0;
}
