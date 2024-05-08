#include "Population.h"
#include <cstddef>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>

Population::Population(
  const size_t& popSize,
  const size_t& numNodes,
  const Netling& founderNetling,
  const double& muNoisePerGene,
  const double& noiseMutDistrMean,
  const double& noiseMutDistrSD,
  const double& recRate)
  : popSize_(popSize),
    numNodes_(numNodes),
    netlings(std::vector<Netling>(popSize, founderNetling)),
    pairedCell_(0),
    recBreakpoint_(0),
    sumF_(0),
    popfitness_(std::vector<double>(popSize, 0)),
    popfitnessCDF_(std::vector<double>(popSize, 0)),
    randvals0to1_(std::vector<double>(popSize, 0)),
    reprGenotypesIndices_(std::vector<int>(popSize, 0)),
    parentNetlings_(std::vector<Netling>(popSize, founderNetling)),
    mutEvent_(0),
    recEvent_(0),
    recombCellAllele_(0),
    mutEventBernDistr(std::bernoulli_distribution (muNoisePerGene)),
    mutValNormDistr(std::normal_distribution<double> (noiseMutDistrMean, noiseMutDistrSD)),
    unifDistr0to1(std::uniform_real_distribution<float> (0.0, 1.0)),
    recEventBernDistr(std::bernoulli_distribution (recRate)),
    recPairUnifDistr(std::uniform_int_distribution<int> (0, static_cast<int>(popSize)-1)),
    recBreakpointUnifDistr(std::uniform_int_distribution<int> (1, static_cast<int>(numNodes)-2))
  {
}

void Population::reproduceNoiseGenotypes(std::mt19937 &PRNG_ENG)
{
  parentNetlings_ = netlings;
  // extract fitnesses from the population into a vector
  for(size_t cell = 0; cell < popSize_; cell++)
  {
    popfitness_[cell] = netlings[cell].fitness;
  }
  // make CDF for population fitness
  popfitnessCDF_[0] = popfitness_[0];
  sumF_ = popfitness_[0];
  for(size_t cell = 1; cell < popSize_; cell++)
  {
    popfitnessCDF_[cell] = popfitnessCDF_[cell-1] + popfitness_[cell];
    sumF_ = sumF_ + popfitness_[cell];
  }

  for(size_t cell = 0; cell < popSize_; cell++)
  {
    randvals0to1_[cell] = unifDistr0to1(PRNG_ENG);
    randvals0to1_[cell] = randvals0to1_[cell]*sumF_;
  }

  // FASTER
  /*
  for(size_t cell = 0; cell < popSize_; cell++){
    auto samp_ = std::upper_bound(std::begin(popfitnessCDF_), std::end(popfitnessCDF_), randvals0to1_[cell]);
    auto idx_ = samp_ - std::begin(popfitnessCDF_);
    netlings[cell] = parentNetlings_[idx_];
  }
  */
  /*
  // new
  for(size_t cell = 0; cell < popSize_; cell++){
    for(size_t l = 0; l < popSize_; l++){
      if (randvals0to1_[cell] < popfitnessCDF_[l]){
        netlings[cell] = parentNetlings_[static_cast<int>(l)];
        break;
      }
    }
  }
  */

  // verified version
  for(size_t cell = 0; cell < popSize_; cell++)
  {
    for(size_t l = 0; l < popSize_; l++)
    {
      if(randvals0to1_[cell] < popfitnessCDF_[l])
      {
        reprGenotypesIndices_[cell] = static_cast<int>(l);
        break;
      }
    }
  }
  for(size_t cell = 0; cell < popSize_; cell++)
  {
    netlings[cell] = parentNetlings_[reprGenotypesIndices_[cell]];
  }
}


void Population::mutateNoiseGenotypes(std::mt19937 &PRNG_ENG)
{
  for(size_t cell = 0; cell < popSize_; cell++)
  {
    for(size_t j = 0; j < numNodes_; j++)
    {
      // draw bool whether there's a mutation event
      mutEvent_ = mutEventBernDistr(PRNG_ENG);
      // if there is a mutation event
      if(mutEvent_)
      {
        // draw mutation (new noise param value)
        netlings[cell].genotNoise[j] = mutValNormDistr(PRNG_ENG);
        // in case drawn value is negative, redraw
        while(netlings[cell].genotNoise[j] < 0)
        {
          netlings[cell].genotNoise[j] = mutValNormDistr(PRNG_ENG);
        }
      }
      // if there is no mutation event, do nothing
      else{}
    }
  }
}


void Population::recombineNoiseGenotypes(std::mt19937 &PRNG_ENG)
{
  for(size_t recombCell = 0; recombCell < popSize_; recombCell++)
  {
    // draw bool whether there's a recombination event
    recEvent_ = recEventBernDistr(PRNG_ENG);
    // if there is a recombination event
    if(recEvent_)
    {
      // draw other cell in the recombining pair
      pairedCell_ = static_cast<size_t>(recPairUnifDistr(PRNG_ENG));
      // in case the drawn cell is the same as the first cell, redraw
      while(pairedCell_ == recombCell)
      {
        pairedCell_ = static_cast<size_t>(recPairUnifDistr(PRNG_ENG));
      }
      // draw genome breakpoint
      recBreakpoint_ = static_cast<size_t>(recBreakpointUnifDistr(PRNG_ENG));
      // recombine genomes
      for(size_t j = 0; j < recBreakpoint_; j++)
      {
        recombCellAllele_ = netlings[recombCell].genotNoise[j];
        netlings[recombCell].genotNoise[j] = netlings[pairedCell_].genotNoise[j];
        netlings[pairedCell_].genotNoise[j] = recombCellAllele_;
      }
    }
    // if there is no recombination event, do nothing
    else{}
  }
}

// write out distribution states
void Population::writeDistrStates()
{
  std::ofstream distrStateFile_mutEvent;
  distrStateFile_mutEvent.open("prng.mutEventDistr.state", std::ios::binary);
  distrStateFile_mutEvent << mutEventBernDistr;
  distrStateFile_mutEvent.close();

  std::ofstream distrStateFile_mutVal;
  distrStateFile_mutVal.open("prng.mutValDistr.state", std::ios::binary);
  distrStateFile_mutVal << mutValNormDistr;
  distrStateFile_mutVal.close();

  std::ofstream distrStateFile_unifDistr;
  distrStateFile_unifDistr.open("prng.unifDistr.state", std::ios::binary);
  distrStateFile_unifDistr << unifDistr0to1;
  distrStateFile_unifDistr.close();

  std::ofstream distrStateFile_recEvent;
  distrStateFile_recEvent.open("prng.recEvent.state", std::ios::binary);
  distrStateFile_recEvent << recEventBernDistr;
  distrStateFile_recEvent.close();

  std::ofstream distrStateFile_recPair;
  distrStateFile_recPair.open("prng.recPair.state", std::ios::binary);
  distrStateFile_recPair << recPairUnifDistr;
  distrStateFile_recPair.close();

  std::ofstream distrStateFile_recBreakpoint;
  distrStateFile_recBreakpoint.open("prng.recBreakpoint.state", std::ios::binary);
  distrStateFile_recBreakpoint << recBreakpointUnifDistr;
  distrStateFile_recBreakpoint.close();
}

// read and assign distribution states
void Population::contRun_readDistrStates()
{
  std::ifstream prevDistrStateFile_mutEvent("prng.mutEventDistr.state", std::ios::binary);
  prevDistrStateFile_mutEvent >> mutEventBernDistr;
  prevDistrStateFile_mutEvent.close();

  std::ifstream prevDistrStateFile_mutVal("prng.mutValDistr.state", std::ios::binary);
  prevDistrStateFile_mutVal >> mutValNormDistr;
  prevDistrStateFile_mutVal.close();

  std::ifstream distrStateFile_unifDistr("prng.unifDistr.state", std::ios::binary);
  distrStateFile_unifDistr >> unifDistr0to1;
  distrStateFile_unifDistr.close();

  std::ifstream distrStateFile_recEvent("prng.recEvent.state", std::ios::binary);
  distrStateFile_recEvent >> recEventBernDistr;
  distrStateFile_recEvent.close();

  std::ifstream distrStateFile_recPair("prng.recPair.state", std::ios::binary);
  distrStateFile_recPair >> recPairUnifDistr;
  distrStateFile_recPair.close();

  std::ifstream distrStateFile_recBreakpoint("prng.recBreakpoint.state", std::ios::binary);
  distrStateFile_recBreakpoint >> recBreakpointUnifDistr;
  distrStateFile_recBreakpoint.close();
}
