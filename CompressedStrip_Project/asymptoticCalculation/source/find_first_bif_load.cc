/*
 * find_first_bif_load.cc
 *
 *  Created on: Nov 1, 2017
 *      Author: Andrew Akerson
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

  return(0);
}