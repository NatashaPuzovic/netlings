#pragma once
#include <cstddef>
#include <vector>
#include <random>

#include "Netling.h"

class Population
{
  private:
    size_t popSize_;
    size_t numNodes_;
  public:
    std::vector<Netling> netlings;
  private:
    size_t pairedCell_;
    size_t recBreakpoint_;
    double sumF_;
    std::vector<double> popfitness_;
    std::vector<double> popfitnessCDF_;
    std::vector<double> randvals0to1_;
    std::vector<int> reprGenotypesIndices_;
    std::vector<Netling> parentNetlings_;
    bool mutEvent_;
    bool recEvent_;
    double recombCellAllele_;
    std::bernoulli_distribution mutEventBernDistr;
    std::normal_distribution<double> mutValNormDistr;
    std::uniform_real_distribution<float> unifDistr0to1;
    std::bernoulli_distribution recEventBernDistr;
    std::uniform_int_distribution<int> recPairUnifDistr;
    std::uniform_int_distribution<int> recBreakpointUnifDistr;

  public:
    Population(
      const size_t& popSize,
      const size_t& numNodes,
      const Netling& founderNetling,
      const double& muNoisePerGene,
      const double& noiseMutDistrMean,
      const double& noiseMutDistrSD,
      const double& recRate);
    void reproduceNoiseGenotypes(std::mt19937 &PRNG_ENG);
    void mutateNoiseGenotypes(std::mt19937 &PRNG_ENG);
    void recombineNoiseGenotypes(std::mt19937 &PRNG_ENG);
    void writeDistrStates();
    void contRun_readDistrStates();
};
