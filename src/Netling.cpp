#include "Netling.h"
#include <cstddef>
#include <vector>
#include <math.h>
#include <random>

Netling::Netling(
  const int& sMax,
  const size_t& numTimesteps,
  const size_t& numNodes,
  const size_t& popSize,
  const double& noiseExt,
  const std::vector<double>& sBasal,
  const std::vector<std::vector<double>>& constrRmtrx,
  const std::vector<double>& constrNoiseInt)
  : sMax_(sMax),
    numTimesteps_(numTimesteps),
    numNodes_(numNodes),
    popSize_(popSize),
    noiseExt_(noiseExt),
    sBasal_(sBasal),
    genotRmtrx(constrRmtrx),
    genotNoise(constrNoiseInt),
    phenot(std::vector<double>(numNodes, 0)),
    fitness(0),
    stdDevs_(std::vector<double>(numNodes, 0)),
    sT_(std::vector<double>(numNodes, 0)),
    sTnext_(std::vector<double>(numNodes, 0)),
    actRateJ_(0),
    weightedDistances_(std::vector<double>(numNodes, 0)),
    sumDist_(0)
  {
};

void Netling::realizePhenotype(std::mt19937 &PRNG_ENG){
  // calculate stand. dev. of expression level distributions for each node
  for (size_t i = 0; i < numNodes_; i++)
  {
    stdDevs_[i] = sqrt(genotNoise[i] + noiseExt_);
  }
  // initialize sT_ vector with basal expression levels
  sT_ = sBasal_;
  //sTnext_.reserve(numNodes_);

  // loop over all numTimesteps_
    for(size_t t = 0; t < numTimesteps_; t++)
    {
      // fill matrix with current sT_
      //for (size_t j = 0; j < numNodes_; j++) {
      //  realization_[t][j] = sT_[j];
      //}
      // draw expression levels in next timestep
      for(size_t j = 0; j < numNodes_; j++)
      {
        // calculate activation rates in this num_timestep
        actRateJ_ = 0;
        for(size_t i = 0; i < numNodes_; i++)
        {
          actRateJ_ = actRateJ_ + sT_[i]*genotRmtrx[j][i];
        }
        // instance of class std::normal_distribution with specific mean and stddev
        std::normal_distribution<double> exprLevelNormDistr(actRateJ_ + sBasal_[j], stdDevs_[j]);
        //
        sTnext_[j] = exprLevelNormDistr(PRNG_ENG);
        // if expr. sT_ is below 0 or above sMax_, change to 0 or sMax_
        if(sTnext_[j] < 0)
        {
          sTnext_[j] = 0;
        }
        else if(sTnext_[j] > sMax_)
        {
          sTnext_[j] = sMax_;
        }
      }
      // assign new sT_ values
      sT_ = sTnext_;
    }
  for(size_t j = 0; j < numNodes_; j++)
  {
    //phenot[j] = realization_[numTimesteps_ - 1][j];
    phenot[j] = sT_[j];
  }
};

void Netling::calculateFitness(const std::vector<double>& sOpt, const std::vector<double>& selPress){
  sumDist_ = 0;
  for(size_t j = 0; j < numNodes_; j++)
  {
    // as in Laarits:
    // numerator = absolute(targetExpressionLevel_i - phenotype_i)
    // denominator = number of nodes * selectionStrength_i
    weightedDistances_[j] = fabs(sOpt[j] - phenot[j])/(static_cast<double>(numNodes_) * selPress[j]);
    // weightedDistances summed up over all genes
    sumDist_ = sumDist_ + weightedDistances_[j];
  }
  // fitness of a given phenotype = e^(-sumDist_)
  fitness = exp(-sumDist_);
};
