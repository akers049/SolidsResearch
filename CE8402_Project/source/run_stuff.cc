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

  // output inital mesh
  ep.output_results (0);



  ep.set_load_val(0.0);
  double load_step = ep.get_final_load()/ep.get_load_steps();
    std::cout << load_step << "\n";
  for(unsigned int i = 0; i < ep.get_load_steps(); i++)
  {

    // update lambda
    ep.increment_load_val(load_step);
    std::cout << "Load step "<< i + 1 << " With loading parameter lambda = " << ep.get_load_val() << std::endl;

    ep.newton_iterate();


    // output data if we're on the right step.
    if((i+1)%ep.get_output_every() == 0)
      ep.output_results ((i+1)/ep.get_output_every());
  }

  return(0);
}
