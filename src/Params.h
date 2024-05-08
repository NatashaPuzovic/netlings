#pragma once
#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

class Params
{
  public:
    int S_MAX;
    size_t NUM_TIMESTEPS;
    size_t NUM_NODES;
    size_t POP_SIZE;
    size_t NUM_GENERATIONS;
    double NOISE_EXT;
    double MU_NOISE_PER_GENE;
    double NOISE_MUT_DISTR_MEAN;
    double NOISE_MUT_DISTR_SD;
    double REC_RATE;
    double S_BASAL_ALL;
    double STARTING_NOISE_INT_ALL;
    double SELPRESS_NOISE_EVOL_ALL;
    int NUM_OUT_GENS;
    unsigned int PRNG_SEED;

  public:
    Params(std::string paramFileName);
    void printAllParams();
};
