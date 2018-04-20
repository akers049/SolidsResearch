/*
 * find_first_bif_load.cc
 *
 *  Created on: Nov 1, 2017
 *      Author: Andrew Akerson
 */

#include "CompressedStrip_asymptotics.h"


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


  ep.assemble_asymptotic_integrals();

  double bif_val_1 = -0.5*ep.E_u1u1u1/(ep.dEdlambda_u1u1);
  double bif_val_2 = -(1.0/3.0)*(ep.E_u1u1u1u1 + 3.0*ep.E_u2u1u1)/(ep.dEdlambda_u1u1);

  std::cout << "\n\n LAMBDA_1 : " << bif_val_1 << std::endl;
  std::cout << " LAMBDA_2 : " <<  bif_val_2 << std::endl;

  return(0);
}
