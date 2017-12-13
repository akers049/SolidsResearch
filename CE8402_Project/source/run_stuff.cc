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


  double nextU;
  double nextV;
  double next_load;
  double next_lambda;

  std::vector<double> u(ep.get_load_steps());
  std::vector<double> v(ep.get_load_steps());
  std::vector<double> lambda(ep.get_load_steps());
  std::vector<double> lambda_scaled(ep.get_load_steps());

  std::vector<double> next_sigma(3);
  std::vector<double> next_secondPiola(3);

  std::vector<std::vector<double>> sigma(ep.get_load_steps());
  std::vector<std::vector<double>> secondPiola(ep.get_load_steps());

  for(unsigned int i = 0; i < ep.get_load_steps(); i++)
  {

    // update lambda
    ep.increment_load_val(load_step);
    std::cout << "Load step "<< i + 1 << " With loading parameter lambda = " << ep.get_load_val() << std::endl;

    ep.newton_iterate();

    ep.get_characteristic_displacements_and_load(&next_lambda, &nextU, &nextV);
    ep.get_stress_components(&next_load, &next_sigma, &next_secondPiola);

    u[i] = nextU;
    v[i] = nextV;
    sigma[i] = next_sigma;
    secondPiola[i] = next_secondPiola;
    lambda_scaled[i] = next_lambda;
    lambda[i] = next_load;

    // output data if we're on the right step.
    if((i+1)%ep.get_output_every() == 0)
      ep.output_results ((i+1)/ep.get_output_every());
  }

  ep.output_load_info(lambda_scaled, u, v, 1);
  ep.output_stresses(lambda, &sigma, &secondPiola, 1);

  return(0);
}
