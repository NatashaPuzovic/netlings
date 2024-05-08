#include "Params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>

Params::Params(
  std::string paramFileName)
  : S_MAX(0),
    NUM_TIMESTEPS(0),
    NUM_NODES(0),
    POP_SIZE(0),
    NUM_GENERATIONS(0),
    NOISE_EXT(0),
    MU_NOISE_PER_GENE(0),
    NOISE_MUT_DISTR_MEAN(0),
    NOISE_MUT_DISTR_SD(0),
    REC_RATE(0),
    S_BASAL_ALL(0),
    STARTING_NOISE_INT_ALL(0),
    SELPRESS_NOISE_EVOL_ALL(0),
    NUM_OUT_GENS(0),
    PRNG_SEED(0)
{
  // read parameters
  std::ifstream paramFile(paramFileName.c_str());
  std::string paramName;
  double paramValue;
  std::map <std::string, double> paramMap;
  while(paramFile >> paramName >> paramValue)
  {
      paramMap[paramName] = paramValue;
  }
  // assign read values
  S_MAX = static_cast<int>(paramMap["S_MAX"]);
  NUM_TIMESTEPS = static_cast<size_t>(paramMap["NUM_TIMESTEPS"]);
  NUM_NODES = static_cast<size_t>(paramMap["NUM_NODES"]);
  POP_SIZE = static_cast<size_t>(paramMap["POP_SIZE"]);
  NUM_GENERATIONS = static_cast<size_t>(paramMap["NUM_GENERATIONS"]);
  NOISE_EXT = paramMap["NOISE_EXT"];
  MU_NOISE_PER_GENE = paramMap["MU_NOISE_PER_GENE"];
  NOISE_MUT_DISTR_MEAN = paramMap["NOISE_MUT_DISTR_MEAN"];
  NOISE_MUT_DISTR_SD = paramMap["NOISE_MUT_DISTR_SD"];
  REC_RATE = paramMap["REC_RATE"];
  S_BASAL_ALL = paramMap["S_BASAL_ALL"];
  STARTING_NOISE_INT_ALL = paramMap["STARTING_NOISE_INT_ALL"];
  SELPRESS_NOISE_EVOL_ALL = paramMap["SELPRESS_NOISE_EVOL_ALL"];
  NUM_OUT_GENS = static_cast<int>(paramMap["NUM_OUT_GENS"]);
  PRNG_SEED = static_cast<unsigned int>(paramMap["PRNG_SEED"]);
}

void Params::printAllParams()
{
  std::cout << "S_MAX: " << S_MAX << std::endl;
  std::cout << "NUM_TIMESTEPS: " << NUM_TIMESTEPS << std::endl;
  std::cout << "NUM_NODES: " << NUM_NODES << std::endl;
  std::cout << "POP_SIZE: " << POP_SIZE << std::endl;
  std::cout << "NUM_GENERATIONS: " << NUM_GENERATIONS << std::endl;
  std::cout << "NOISE_EXT: " << NOISE_EXT << std::endl;
  std::cout << "MU_NOISE_PER_GENE: " << MU_NOISE_PER_GENE << std::endl;
  std::cout << "NOISE_MUT_DISTR_MEAN: " << NOISE_MUT_DISTR_MEAN << std::endl;
  std::cout << "NOISE_MUT_DISTR_SD: " << NOISE_MUT_DISTR_SD << std::endl;
  std::cout << "REC_RATE: " << REC_RATE << std::endl;
  std::cout << "S_BASAL_ALL: " << S_BASAL_ALL << std::endl;
  std::cout << "STARTING_NOISE_INT_ALL: " << STARTING_NOISE_INT_ALL << std::endl;
  std::cout << "SELPRESS_NOISE_EVOL_ALL: " << SELPRESS_NOISE_EVOL_ALL << std::endl;
  std::cout << "NUM_OUT_GENS: " << NUM_OUT_GENS << std::endl;
  std::cout << "PRNG_SEED: " << PRNG_SEED << std::endl;
}
