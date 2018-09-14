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

  ep.set_present_lambda(ep.critical_lambda_analytical);
  ep.update_F0(ep.critical_lambda_analytical);

  ep.assemble_asymptotic_integrals();

//  ep.check_W_derivs();

  ep.assemble_vexex_eq();
  ep.output_results(0);
//  for(unsigned int i = 0; i < ep.vexex_eq.size(); i++)
//  {
//    std::cout << ep.vexex_eq[i] << std::endl;
//  }

  std::cout << "\n\n LAMBDA_1 : " << ep.get_deltaLambda1() << std::endl;
  std::cout << " LAMBDA_2 : " <<  ep.get_deltaLambda2() << std::endl;
  std::cout << " BIG_LAMBDA_2 : " << ep.get_deltaConguateLambda2() << std::endl;

  FILE* fid;

  char newFileName[MAXLINE];
  char *token;

  token = strtok(fileName, "/");
  token = strtok(NULL, ".");

  sprintf(newFileName, "%s", ep.output_directory);
  strcat(newFileName, "/");
  strcat(newFileName, token);
  strcat(newFileName, ".out");

  std::cout << "Printing File to : " << newFileName << std::endl;

  fid = std::fopen(newFileName, "w");
  if (fid == NULL)
  {
    std::cout << "Unable to open file \"" << fileName  << "\"" <<  std::endl;
  }
  else
  {
    fprintf(fid, "%.16f\n%.16f\n%.16f", ep.get_kappa() ,ep.get_deltaLambda2(),ep.get_deltaConguateLambda2());
  }



//  std::cout << ep.E_u1u1u1u1 << std::endl << ep.E_u2u1u1 << std::endl << ep.dEdlambda_u1u1 << std::endl << ep.E_u2u2 << std::endl;

  return(0);
}
