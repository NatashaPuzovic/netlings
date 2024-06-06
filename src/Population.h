#pragma once
#include <cstddef>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>

#include "Netling.h"

class Population
{
private:
    size_t popSize_;
    size_t numNodes_;
    size_t pairedCell_;
    size_t recBreakpoint1_;
    size_t recBreakpoint2_;
    double sumF_;
    std::vector<double> popfitness_;
    std::vector<double> popfitnessCDF_;
    std::vector<double> randvals0to1_;
    std::vector<int> reprGenotypesIndices_;
    std::vector<Netling> parentNetlings_;
    bool mutEvent_;
    bool recEvent_;
    double recombCellAllele_;
    std::bernoulli_distribution mutationEventDistributionNoise_;
    std::uniform_real_distribution<double> mutationValueDistributionNoise_;
    std::bernoulli_distribution mutationEventDistributionSBasal_;
    std::uniform_real_distribution<double> mutationValueDistributionSBasal_;
    std::bernoulli_distribution mutationEventDistributionNet_;
    std::uniform_real_distribution<double> mutationValueDistributionNet_;
    std::uniform_real_distribution<float> unifDistr0to1_;
    std::bernoulli_distribution recEventBernDistr_;
    std::uniform_int_distribution<int> recPairUnifDistr_;
    std::uniform_int_distribution<int> recBreakpointUnifDistr_;
    std::ofstream outSummarizedPVFFile_;
    std::ofstream outNoisesFile_;
    std::ofstream outSBasalsFile_;
    std::ofstream outNetsFile_;
    //std::ofstream outGenealogyFile_;
    std::ofstream outEntirePopFile_;
    double sumEta_j_;
    double meanEta_j_;
    double sumSBasal_j_;
    double meanSBasal_j_;
public:
    std::vector<Netling> netlings;
    int lastAssignedIDGenotype;
    int lastAssignedIDNoise;
    int lastAssignedIDSBasal;
    int lastAssignedIDNet;

public:
    Population(
        const size_t& popSize,
        const size_t& numNodes,
        const Netling& founderNetling,
        const double& mutationRateNoisePerGene,
        const double& mutationDistrNoiseLowerBound,
        const double& mutationDistrNoiseUpperBound,
        const double& mutationRateSBasalPerGene,
        const double& mutationDistrSBasalLowerBound,
        const double& mutationDistrSBasalUpperBound,
        const double& mutationRateNetPerGene,
        const double& mutationDistrNetLowerBound,
        const double& mutationDistrNetUpperBound,
        const double& recRate,
        const std::string& outFileNameTemplate);

    void readFirstGenEntirePopulation(const std::string& passedFirstGenEntirePopulationFileName);
    void reproduceNetlings(std::mt19937& PRNG_ENG);

    void mutateNoiseGenotypes(std::mt19937& PRNG_ENG);
    void recombineNoiseGenotypes(std::mt19937& PRNG_ENG);
    void outputSummarizedNoiseGenotypes(size_t gen);

    void mutateSBasalGenotypes(std::mt19937& PRNG_ENG);
    void recombineSBasalGenotypes(std::mt19937& PRNG_ENG);
    void outputSummarizedSBasalGenotypes(size_t gen);

    void mutateNetExistingRegLinks(std::mt19937& PRNG_ENG);
    //void outputNetworkTopo(size_t gen);

    void outputSummarizedPVF(size_t gen);
    void outputEntirePopulation(size_t gen);

    void outputDistrStates();
    void printFinishMessage(std::string outFileNameTemplate);
    void closeOutputFiles();
};
