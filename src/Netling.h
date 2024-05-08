#pragma once
#include <cstddef>
#include <vector>
#include <random>


class Netling
{
  private:
    int sMax_;
    size_t numTimesteps_;
    size_t numNodes_;
    size_t popSize_;
    double noiseExt_;
    std::vector<double> sBasal_;
  public:
    std::vector<std::vector<double>> genotRmtrx;
    std::vector<double> genotNoise;
    std::vector<double> phenot;
    double fitness;
  private:
    std::vector<double> stdDevs_;
    std::vector<double> sT_;
    std::vector<double> sTnext_;
    double actRateJ_;
    std::vector<double> weightedDistances_;
    double sumDist_;

  public:
    Netling(
      const int& sMax,
      const size_t& numTimesteps,
      const size_t& numNodes,
      const size_t& popSize,
      const double& noiseExt,
      const std::vector<double>& sBasal,
      const std::vector<std::vector<double>>& constrRmtrx,
      const std::vector<double>& constrNoiseInt);
    void realizePhenotype(std::mt19937 &PRNG_ENG);
    void calculateFitness(const std::vector<double>& sOpt, const std::vector<double>& selPress);
};
