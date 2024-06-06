#include "Netling.h"
#include <cstddef>
#include <vector>
#include <math.h>
#include <random>

//#include <chrono>
//#include <iostream>

int Netling::sMax;
size_t Netling::numTimesteps;
size_t Netling::numNodes;
double Netling::noiseExt;

Netling::Netling(
    const std::vector<double>& constrNoiseInt,
    const std::vector<double>& constrSBasal,
    const std::vector<std::vector<double>>& constrRmtrx)
    : genotNoise(constrNoiseInt),
    genotSBasal(constrSBasal),
    genotRmtrx(constrRmtrx),
    IDGenotype(1),
    IDNet(1),
    IDNoise(1),
    IDSBasal(1),
    phenot(std::vector<double>(numNodes, 0)),
    fitness(0),
    stdDevs_(std::vector<double>(numNodes, 0)),
    sT_(std::vector<double>(numNodes, 0)),
    sTnext_(std::vector<double>(numNodes, 0)),
    actRateJ_(0),
    weightedDistances_(std::vector<double>(numNodes, 0)),
    sumDist_(0),
    realizationNormDistribution_(std::normal_distribution<double>(0, 1))
{
}

void Netling::realizePhenotype(std::mt19937& PRNG_ENG)
{
    // calculate stand. dev. of expression level distributions for each node
    for (size_t i = 0; i < numNodes; ++i)
    {
        stdDevs_[i] = sqrt(genotNoise[i] + noiseExt);
    }
    // initialize sT_ vector with basal expression levels
    sT_ = genotSBasal;
    //sTnext_.reserve(numNodes);

    // loop over all numTimesteps
    for (size_t t = 0; t < numTimesteps; ++t)
    {
        // fill matrix with current sT_
        //for (size_t j = 0; j < numNodes; j++) {
        //  realization_[t][j] = sT_[j];
        //}
        // draw expression levels in next timestep
        for (size_t j = 0; j < numNodes; ++j)
        {
            // calculate activation rates in this num_timestep
            actRateJ_ = 0;
            auto& tmp = genotRmtrx[j];
            for (size_t i = 0; i < numNodes; ++i)
            {
                //actRateJ_ += sT_[i]*genotRmtrx[j][i];
                actRateJ_ += sT_[i] * tmp[i];
            }
            // instance of class std::normal_distribution with specific mean and stddev
            //std::normal_distribution<double> exprLevelNormDistr(actRateJ_ + genotSBasal[j], stdDevs_[j]);
            //sTnext_[j] = exprLevelNormDistr(PRNG_ENG);

            sTnext_[j] = realizationNormDistribution_(PRNG_ENG);
            sTnext_[j] = sTnext_[j] * stdDevs_[j] + actRateJ_ + genotSBasal[j];

            // if expr. sT_ is below 0 or above sMax, change to 0 or sMax
            if (sTnext_[j] < 0)
            {
                sTnext_[j] = 0;
            }
            else if (sTnext_[j] > sMax)
            {
                sTnext_[j] = sMax;
            }
        }
        // assign new sT_ values
        sT_ = sTnext_;
    }
    for (size_t j = 0; j < numNodes; ++j)
    {
        //phenot[j] = realization_[numTimesteps - 1][j];
        phenot[j] = sT_[j];
    }
}

void Netling::calculateFitness(const std::vector<double>& sOpt, const std::vector<double>& selPress)
{
    sumDist_ = 0;
    for (size_t j = 0; j < numNodes; j++)
    {
        // as in Laarits:
        // numerator = absolute(targetExpressionLevel_i - phenotype_i)
        // denominator = number of nodes * selectionStrength_i

        // if selective pressure is 0, the genes is not selected and its phenotype is not important for fitness
        if (selPress[j] == 0) {
            continue;
        }
        else {
            weightedDistances_[j] = fabs(sOpt[j] - phenot[j]) / (static_cast<double>(numNodes) * selPress[j]);
            // weightedDistances summed up over all genes
            sumDist_ = sumDist_ + weightedDistances_[j];
        }
    }
    // fitness of a given phenotype = e^(-sumDist_)
    fitness = exp(-sumDist_);
}
