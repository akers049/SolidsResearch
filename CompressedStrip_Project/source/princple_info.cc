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



  std::vector<double> lambda_values;
  std::vector<double> congugate_lambda_values;
  std::vector<double> energy_values;
  std::vector<double> displacement_magnitude;
  lambda_values.push_back(0.0);
  congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
  energy_values.push_back( ep.system_energy/ep.get_number_unit_cells());
  displacement_magnitude.push_back(ep.present_solution.l2_norm()/ep.get_number_unit_cells());



  for(unsigned int i = 0; i < 1000; i++)
  {
    double lambda_eval = ((ep.critical_lambda_analytical + 0.1)*i)/1000.0;
    ep.set_present_lambda( lambda_eval );

    // get energy and congugate lambda value and save them.
    ep.assemble_system_energy_and_congugate_lambda();
    lambda_values.push_back(ep.get_present_lambda());
    congugate_lambda_values.push_back(ep.congugate_lambda/ep.get_number_unit_cells());
    energy_values.push_back(ep.system_energy/ep.get_number_unit_cells());
    displacement_magnitude.push_back(ep.present_solution.l2_norm()/(sqrt(1.0*ep.get_number_unit_cells())));
    ep.output_results(i+5000);


  }
  ep.output_load_info(lambda_values, energy_values, congugate_lambda_values, displacement_magnitude, 4);


  return 0;
}
