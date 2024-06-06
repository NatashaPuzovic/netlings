#pragma once
#include <cstddef>
#include <vector>
#include <random>


class Netling
{
public:
    std::vector<double> genotNoise;
    std::vector<double> genotSBasal;
    std::vector<std::vector<double>> genotRmtrx;
    int IDGenotype;
    int IDNet;
    int IDNoise;
    int IDSBasal;
    static int sMax;
    static size_t numTimesteps;
    static size_t numNodes;
    static double noiseExt;
    std::vector<double> phenot;
    double fitness;

private:
    std::vector<double> stdDevs_;
    std::vector<double> sT_;
    std::vector<double> sTnext_;
    double actRateJ_;
    std::vector<double> weightedDistances_;
    double sumDist_;
    std::normal_distribution<double> realizationNormDistribution_;

public:
    Netling(
        const std::vector<double>& constrNoiseInt,
        const std::vector<double>& constrSBasal,
        const std::vector<std::vector<double>>& constrRmtrx);
    void realizePhenotype(std::mt19937& PRNG_ENG);
    void calculateFitness(const std::vector<double>& sOpt, const std::vector<double>& selPress);
};
