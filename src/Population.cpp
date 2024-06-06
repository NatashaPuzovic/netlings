#include "Population.h"
#include <cstddef>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>

//#include <chrono>
//#include <iostream>

Population::Population(
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
    const std::string& outFileNameTemplate)
    : popSize_(popSize),
    numNodes_(numNodes),
    pairedCell_(0),
    recBreakpoint1_(0),
    recBreakpoint2_(0),
    sumF_(0),
    popfitness_(std::vector<double>(popSize, 0)),
    popfitnessCDF_(std::vector<double>(popSize, 0)),
    randvals0to1_(std::vector<double>(popSize, 0)),
    reprGenotypesIndices_(std::vector<int>(popSize, 0)),
    parentNetlings_(std::vector<Netling>(popSize, founderNetling)),
    mutEvent_(0),
    recEvent_(0),
    recombCellAllele_(0),
    mutationEventDistributionNoise_(std::bernoulli_distribution(mutationRateNoisePerGene)),
    mutationValueDistributionNoise_(std::uniform_real_distribution<double>(mutationDistrNoiseLowerBound, mutationDistrNoiseUpperBound)),
    mutationEventDistributionSBasal_(std::bernoulli_distribution(mutationRateSBasalPerGene)),
    mutationValueDistributionSBasal_(std::uniform_real_distribution<double>(mutationDistrSBasalLowerBound, mutationDistrSBasalUpperBound)),
    mutationEventDistributionNet_(std::bernoulli_distribution(mutationRateNetPerGene)),
    mutationValueDistributionNet_(std::uniform_real_distribution<double>(mutationDistrNetLowerBound, mutationDistrNetUpperBound)),
    unifDistr0to1_(std::uniform_real_distribution<float>(0.0, 1.0)),
    recEventBernDistr_(std::bernoulli_distribution(recRate)),
    recPairUnifDistr_(std::uniform_int_distribution<int>(0, static_cast<int>(popSize) - 1)),
    recBreakpointUnifDistr_(std::uniform_int_distribution<int>(1, static_cast<int>(numNodes) - 2)),
    outSummarizedPVFFile_((outFileNameTemplate + ".hist").c_str(), std::fstream::out),
    outNoisesFile_((outFileNameTemplate + ".noises").c_str(), std::fstream::out),
    outSBasalsFile_((outFileNameTemplate + ".sbasals").c_str(), std::fstream::out),
    outNetsFile_((outFileNameTemplate + ".nets").c_str(), std::fstream::out),
    //outGenealogyFile_((outFileNameTemplate + ".genealogy").c_str(), std::fstream::out),
    outEntirePopFile_((outFileNameTemplate + ".pop").c_str(), std::fstream::out),
    sumEta_j_(0),
    meanEta_j_(0),
    sumSBasal_j_(0),
    meanSBasal_j_(0),
    netlings(std::vector<Netling>(popSize, founderNetling)),
    lastAssignedIDGenotype(1),
    lastAssignedIDNoise(1),
    lastAssignedIDSBasal(1),
    lastAssignedIDNet(1)
{
}

void Population::readFirstGenEntirePopulation(const std::string& passedFirstGenEntirePopulationFileName)
{
    // open connection to file
    std::ifstream firstGenEntirePopulationFile;
    firstGenEntirePopulationFile.open(passedFirstGenEntirePopulationFileName.c_str());

    if (!firstGenEntirePopulationFile.is_open())
    {
        std::cout << "Error: Failed to open first generation file: " << passedFirstGenEntirePopulationFileName << std::endl;
    }
    else
    {
        int generationNumber = 0;
        int cellNumber = 0;
        for (size_t cell = 0; cell < popSize_; ++cell)
        {
            firstGenEntirePopulationFile >> generationNumber;
            firstGenEntirePopulationFile >> cellNumber;
            for (size_t i = 0; i < numNodes_; ++i)
            {
                firstGenEntirePopulationFile >> netlings[cell].genotNoise[i];
            }
            for (size_t i = 0; i < numNodes_; ++i)
            {
                firstGenEntirePopulationFile >> netlings[cell].genotSBasal[i];
            }
            for (size_t i = 0; i < numNodes_; ++i)
            {
                for (size_t j = 0; j < numNodes_; ++j)
                {
                    firstGenEntirePopulationFile >> netlings[cell].genotRmtrx[i][j];
                }
            }
            for (size_t i = 0; i < numNodes_; ++i)
            {
                firstGenEntirePopulationFile >> netlings[cell].phenot[i];
            }
            firstGenEntirePopulationFile >> netlings[cell].fitness;
        }
        firstGenEntirePopulationFile.close();
    }
}

void Population::reproduceNetlings(std::mt19937& PRNG_ENG)
{

    parentNetlings_ = netlings;

    // extract fitnesses from the population into a vector
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        popfitness_[cell] = netlings[cell].fitness;
    }

    // make CDF for population fitness
    popfitnessCDF_[0] = netlings[0].fitness;
    sumF_ = netlings[0].fitness;
    for (size_t cell = 1; cell < popSize_; cell++)
    {
        popfitnessCDF_[cell] = popfitnessCDF_[cell - 1] + popfitness_[cell];
        sumF_ = sumF_ + popfitness_[cell];
    }

    for (size_t cell = 0; cell < popSize_; cell++)
    {
        randvals0to1_[cell] = unifDistr0to1_(PRNG_ENG);
        randvals0to1_[cell] = randvals0to1_[cell] * sumF_;
    }
    // FASTER
    /*
    for(size_t cell = 0; cell < popSize_; cell++){
      auto samp_ = std::upper_bound(std::begin(popfitnessCDF_), std::end(popfitnessCDF_), randvals0to1_[cell]);
      auto idx_ = samp_ - std::begin(popfitnessCDF_);
      netlings[cell] = parentNetlings_[idx_];
    }
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
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        for (size_t l = 0; l < popSize_; l++)
        {
            if (randvals0to1_[cell] < popfitnessCDF_[l])
            {
                reprGenotypesIndices_[cell] = static_cast<int>(l);
                break;
            }
        }
    }
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        netlings[cell] = parentNetlings_[reprGenotypesIndices_[cell]];
    }
}

void Population::mutateNoiseGenotypes(std::mt19937& PRNG_ENG)
{
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        for (size_t i = 0; i < numNodes_; i++)
        {
            // draw bool whether there's a mutation event
            mutEvent_ = mutationEventDistributionNoise_(PRNG_ENG);
            // if there is a mutation event
            if (mutEvent_)
            {
                // draw mutation (new noise param value)
                netlings[cell].genotNoise[i] = mutationValueDistributionNoise_(PRNG_ENG);
                // in case drawn value is negative, redraw
                while (netlings[cell].genotNoise[i] < 0)
                {
                    netlings[cell].genotNoise[i] = mutationValueDistributionNoise_(PRNG_ENG);
                }
                // update the unique ID identifiers
                /*
                lastAssignedIDGenotype++;
                lastAssignedIDNoise++;
                // write genealogy (genotype IDs of the parent and the new mutant child)
                outGenealogyFile_ << gen + 1 << " ";
                outGenealogyFile_ << netlings[cell].IDGenotype << " ";
                netlings[cell].IDGenotype = lastAssignedIDGenotype;
                outGenealogyFile_ << netlings[cell].IDGenotype << " ";
                netlings[cell].IDNoise = lastAssignedIDNoise;
                outGenealogyFile_ << std::endl;
                */
            }
        }
    }
}

void Population::recombineNoiseGenotypes(std::mt19937& PRNG_ENG)
{
    for (size_t recombCell_ = 0; recombCell_ < popSize_; recombCell_++)
    {
        // draw bool whether there's a recombination event
        recEvent_ = recEventBernDistr_(PRNG_ENG);
        // if there is a recombination event
        if (recEvent_)
        {
            // draw other cell in the recombining pair
            pairedCell_ = static_cast<size_t>(recPairUnifDistr_(PRNG_ENG));
            // in case the drawn cell is the same as the first cell, redraw
            while (pairedCell_ == recombCell_)
            {
                pairedCell_ = static_cast<size_t>(recPairUnifDistr_(PRNG_ENG));
            }
            // draw genome recombination breakpoints
            recBreakpoint1_ = static_cast<size_t>(recBreakpointUnifDistr_(PRNG_ENG));
            recBreakpoint2_ = static_cast<size_t>(recBreakpointUnifDistr_(PRNG_ENG));
            // in case the 2nd recombination breakpoint is the same as 1st, redraw
            while (recBreakpoint2_ == recBreakpoint1_)
            {
                recBreakpoint2_ = static_cast<size_t>(recBreakpointUnifDistr_(PRNG_ENG));
            }

            if (recBreakpoint1_ < recBreakpoint2_)
            {
                // recombine genomes
                for (size_t j = recBreakpoint1_; j <= recBreakpoint2_; j++)
                {
                    recombCellAllele_ = netlings[recombCell_].genotNoise[j];
                    netlings[recombCell_].genotNoise[j] = netlings[pairedCell_].genotNoise[j];
                    netlings[pairedCell_].genotNoise[j] = recombCellAllele_;
                }
            }
            else
            {
                // recombine genomes
                for (size_t j = recBreakpoint2_; j <= recBreakpoint1_; j++)
                {
                    recombCellAllele_ = netlings[recombCell_].genotNoise[j];
                    netlings[recombCell_].genotNoise[j] = netlings[pairedCell_].genotNoise[j];
                    netlings[pairedCell_].genotNoise[j] = recombCellAllele_;
                }
            }
        }
    }
}

void Population::outputSummarizedNoiseGenotypes(size_t gen)
{
    // write the mean values of noise genotypes population-wide
    outNoisesFile_ << gen + 1 << " ";
    for (size_t j = 0; j < numNodes_; j++)
    {
        sumEta_j_ = 0;
        for (size_t cell = 0; cell < popSize_; cell++)
        {
            sumEta_j_ += netlings[cell].genotNoise[j];
        }
        meanEta_j_ = sumEta_j_ / static_cast<double>(popSize_);
        outNoisesFile_ << meanEta_j_ << " ";
    }
    outNoisesFile_ << std::endl;
}

void Population::mutateSBasalGenotypes(std::mt19937& PRNG_ENG)
{
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        for (size_t i = 0; i < numNodes_; i++)
        {
            // draw bool whether there's a mutation event
            mutEvent_ = mutationEventDistributionSBasal_(PRNG_ENG);
            // if there is a mutation event
            if (mutEvent_)
            {
                // draw mutation (new noise param value)
                netlings[cell].genotSBasal[i] = mutationValueDistributionSBasal_(PRNG_ENG);
                // in case drawn value is negative, redraw
                while (netlings[cell].genotSBasal[i] < 0)
                {
                    netlings[cell].genotSBasal[i] = mutationValueDistributionSBasal_(PRNG_ENG);
                }
                // update the unique ID identifiers
                /*
                lastAssignedIDGenotype++;
                lastAssignedIDSBasal++;
                // write genealogy (genotype IDs of the parent and the new mutant child)
                outGenealogyFile_ << gen + 1 << " ";
                outGenealogyFile_ << netlings[cell].IDGenotype << " ";
                netlings[cell].IDGenotype = lastAssignedIDGenotype;
                outGenealogyFile_ << netlings[cell].IDSBasal << " ";
                netlings[cell].IDSBasal = lastAssignedIDSBasal;
                outGenealogyFile_ << std::endl;
                */
            }
        }
    }
}

void Population::recombineSBasalGenotypes(std::mt19937& PRNG_ENG)
{
    for (size_t recombCell_ = 0; recombCell_ < popSize_; recombCell_++)
    {
        // draw bool whether there's a recombination event
        recEvent_ = recEventBernDistr_(PRNG_ENG);
        // if there is a recombination event
        if (recEvent_)
        {
            // draw other cell in the recombining pair
            pairedCell_ = static_cast<size_t>(recPairUnifDistr_(PRNG_ENG));
            // in case the drawn cell is the same as the first cell, redraw
            while (pairedCell_ == recombCell_)
            {
                pairedCell_ = static_cast<size_t>(recPairUnifDistr_(PRNG_ENG));
            }
            // draw genome recombination breakpoints
            recBreakpoint1_ = static_cast<size_t>(recBreakpointUnifDistr_(PRNG_ENG));
            recBreakpoint2_ = static_cast<size_t>(recBreakpointUnifDistr_(PRNG_ENG));
            // in case the 2nd recombination breakpoint is the same as 1st, redraw
            while (recBreakpoint2_ == recBreakpoint1_)
            {
                recBreakpoint2_ = static_cast<size_t>(recBreakpointUnifDistr_(PRNG_ENG));
            }
            if (recBreakpoint1_ < recBreakpoint2_)
            {
                // recombine genomes
                for (size_t j = recBreakpoint1_; j <= recBreakpoint2_; j++)
                {
                    recombCellAllele_ = netlings[recombCell_].genotSBasal[j];
                    netlings[recombCell_].genotSBasal[j] = netlings[pairedCell_].genotSBasal[j];
                    netlings[pairedCell_].genotSBasal[j] = recombCellAllele_;
                }
            }
            else
            {
                // recombine genomes
                for (size_t j = recBreakpoint2_; j <= recBreakpoint1_; j++)
                {
                    recombCellAllele_ = netlings[recombCell_].genotSBasal[j];
                    netlings[recombCell_].genotSBasal[j] = netlings[pairedCell_].genotSBasal[j];
                    netlings[pairedCell_].genotSBasal[j] = recombCellAllele_;
                }
            }
        }
    }
}

void Population::outputSummarizedSBasalGenotypes(size_t gen)
{
    // write the mean values of noise genotypes population-wide
    outSBasalsFile_ << gen + 1 << " ";
    for (size_t j = 0; j < numNodes_; j++)
    {
        sumSBasal_j_ = 0;
        for (size_t cell = 0; cell < popSize_; cell++)
        {
            sumSBasal_j_ += netlings[cell].genotSBasal[j];
        }
        meanSBasal_j_ = sumSBasal_j_ / static_cast<double>(popSize_);
        outSBasalsFile_ << meanSBasal_j_ << " ";
    }
    outSBasalsFile_ << std::endl;
}

void Population::mutateNetExistingRegLinks(std::mt19937& PRNG_ENG)
{
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        for (size_t i = 0; i < numNodes_; i++)
        {
            for (size_t j = 0; j < numNodes_; j++)
            {
                if (netlings[cell].genotRmtrx[i][j] != 0)
                {
                    // draw bool whether there's a mutation event
                    mutEvent_ = mutationEventDistributionNet_(PRNG_ENG);
                    // if there is a mutation event
                    if (mutEvent_)
                    {
                        // draw mutation (new regulatory link value)
                        netlings[cell].genotRmtrx[i][j] = mutationValueDistributionNet_(PRNG_ENG);

                        // update the unique ID identifiers
                        //lastAssignedIDGenotype++;
                        //lastAssignedIDNet++;
                        // write genealogy (genotype IDs of the parent and the new mutant child)
                        //outGenealogyFile_ << gen + 1 << " ";
                        //outGenealogyFile_ << netlings[cell].IDGenotype << " ";
                        //netlings[cell].IDGenotype = lastAssignedIDGenotype;
                        //outGenealogyFile_ << netlings[cell].IDNet << " ";
                        //netlings[cell].IDNet = lastAssignedIDNet;
                        //outGenealogyFile_ << std::endl;
                        // write net part of genotypes
                        //outNetsFile_ << netlings[cell].IDNet << " ";
                        //for (size_t k = 0; k < numNodes_; k++)
                        //{
                        //  for (size_t l = 0; l < numNodes_; l++)
                        //  {
                        //    outNetsFile_ << netlings[cell].genotRmtrx[k][l] << " ";
                        //  }
                        //}
                        //outNetsFile_ << std::endl;
                    }
                }
            }
        }
    }
}

/*
void Population::outputNetworkTopo(size_t gen)
{
  // write the net topo configuration part of genotype in first gen
  outNetsFile_ << gen + 1 << " ";
  for (size_t i = 0; i < numNodes_; i++)
  {
    for (size_t j = 0; j < numNodes_; j++)
    {
      outNetsFile_ << netlings[0].genotRmtrx[i][j] << " ";
    }
  }
  outNetsFile_ << std::endl;
}
*/

void Population::outputSummarizedPVF(size_t gen)
{
    //double sumG_j = 0;
    double sumP_j = 0;
    //double meanG_j = 0;
    double varP_j = 0;
    double meanF = 0;
    std::vector<double> meanP(numNodes_, 0);

    outSummarizedPVFFile_ << gen + 1 << " ";
    /*
    // mean genotypes
    for(size_t j = 0; j < prms.NUM_NODES; j++)
    {
      sumG_j = 0;
      for(size_t cell = 0; cell < popSize_; cell++)
      {
        sumG_j += pop.netlings[cell].genotNoise[j];
      }
      meanG_j = sumG_j/static_cast<double>(popSize_);
      outSummarizedPVFfile << meanG_j << " ";
    }
    */
    // mean phenotypes
    for (size_t j = 0; j < numNodes_; j++)
    {
        sumP_j = 0;
        for (size_t cell = 0; cell < popSize_; cell++)
        {
            sumP_j += netlings[cell].phenot[j];
        }
        meanP[j] = sumP_j / static_cast<double>(popSize_);
        outSummarizedPVFFile_ << meanP[j] << " ";
    }
    // variances of phenotypes
    for (size_t j = 0; j < numNodes_; j++)
    {
        varP_j = 0;
        for (size_t cell = 0; cell < popSize_; cell++)
        {
            varP_j += (netlings[cell].phenot[j] - meanP[j]) * (netlings[cell].phenot[j] - meanP[j]);
        }
        varP_j = varP_j / static_cast<double>(popSize_);
        outSummarizedPVFFile_ << varP_j << " ";
    }
    // mean fitness
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        meanF += netlings[cell].fitness;
    }
    meanF = meanF / static_cast<double>(popSize_);
    outSummarizedPVFFile_ << meanF << " ";

    // finish line
    outSummarizedPVFFile_ << std::endl;
}

void Population::outputEntirePopulation(size_t gen)
{
    for (size_t cell = 0; cell < popSize_; cell++)
    {
        outEntirePopFile_ << gen + 1 << " " << cell + 1 << " ";
        //outEntirePopFile_ << netlings[cell].IDGenotype << " ";
        //outEntirePopFile_ << netlings[cell].IDNoise << " ";
        //outEntirePopFile_ << netlings[cell].IDSBasal << " ";
        //outEntirePopFile_ << netlings[cell].IDNet << " ";

        for (size_t j = 0; j < numNodes_; j++)
        {
            outEntirePopFile_ << netlings[cell].genotNoise[j] << " ";
        }
        for (size_t j = 0; j < numNodes_; j++)
        {
            outEntirePopFile_ << netlings[cell].genotSBasal[j] << " ";
        }
        //outNetsFile_ << netlings[cell].IDNet << " ";
        for (size_t k = 0; k < numNodes_; k++)
        {
            for (size_t l = 0; l < numNodes_; l++)
            {
                outEntirePopFile_ << netlings[cell].genotRmtrx[k][l] << " ";
            }
        }
        for (size_t j = 0; j < numNodes_; j++)
        {
            outEntirePopFile_ << netlings[cell].phenot[j] << " ";
        }
        //outNetsFile_ << std::endl;
        outEntirePopFile_ << netlings[cell].fitness;
        outEntirePopFile_ << std::endl;
    }
}

void Population::outputDistrStates()
{
    // write out distribution states
    std::ofstream distrStateFile_mutEvent;
    distrStateFile_mutEvent.open("prng.mutEventDistr.state", std::ios::binary);
    distrStateFile_mutEvent << mutationEventDistributionNoise_;
    distrStateFile_mutEvent.close();

    std::ofstream distrStateFile_mutVal;
    distrStateFile_mutVal.open("prng.mutValDistr.state", std::ios::binary);
    distrStateFile_mutVal << mutationValueDistributionNoise_;
    distrStateFile_mutVal.close();

    std::ofstream distrStateFile_unifDistr;
    distrStateFile_unifDistr.open("prng.unifDistr.state", std::ios::binary);
    distrStateFile_unifDistr << unifDistr0to1_;
    distrStateFile_unifDistr.close();

    std::ofstream distrStateFile_recEvent;
    distrStateFile_recEvent.open("prng.recEvent.state", std::ios::binary);
    distrStateFile_recEvent << recEventBernDistr_;
    distrStateFile_recEvent.close();

    std::ofstream distrStateFile_recPair;
    distrStateFile_recPair.open("prng.recPair.state", std::ios::binary);
    distrStateFile_recPair << recPairUnifDistr_;
    distrStateFile_recPair.close();

    std::ofstream distrStateFile_recBreakpoint;
    distrStateFile_recBreakpoint.open("prng.recBreakpoint.state", std::ios::binary);
    distrStateFile_recBreakpoint << recBreakpointUnifDistr_;
    distrStateFile_recBreakpoint.close();
}

void Population::printFinishMessage(std::string outFileNameTemplate)
{
    std::cout << std::endl << "Written summarized phenotypes, variance of phenotypes and fitness of chosen generations to " <<
        outFileNameTemplate << ".hist." << std::endl;
    std::cout << "Written summarized noise genotypes of chosen generations to " <<
        outFileNameTemplate << ".noises." << std::endl;
    std::cout << "Written summarized sbasal genotypes of chosen generations to " <<
        outFileNameTemplate << ".sbasals." << std::endl;
    std::cout << "Written net genotypes of chosen generations to " <<
        outFileNameTemplate << ".nets." << std::endl;
    //std::cout << "Written genealogy to " <<
    //outFileNameTemplate << ".genealogy." << std::endl;
}

void Population::closeOutputFiles()
{
    if (outSummarizedPVFFile_.is_open())
    {
        outSummarizedPVFFile_.close();
    }
    if (outNetsFile_.is_open())
    {
        outNetsFile_.close();
    }
    if (outNoisesFile_.is_open())
    {
        outNoisesFile_.close();
    }
    if (outSBasalsFile_.is_open())
    {
        outSBasalsFile_.close();
    }
    if (outEntirePopFile_.is_open())
    {
        outEntirePopFile_.close();
    }
}

/*
// read and assign distribution states
void Population::contRun_readDistrStates()
{
  std::ifstream prevDistrStateFile_mutEvent("prng.mutEventDistr.state", std::ios::binary);
  prevDistrStateFile_mutEvent >> mutEventBernDistr_;
  prevDistrStateFile_mutEvent.close();

  std::ifstream prevDistrStateFile_mutVal("prng.mutValDistr.state", std::ios::binary);
  prevDistrStateFile_mutVal >> mutValNormDistr_;
  prevDistrStateFile_mutVal.close();

  std::ifstream distrStateFile_unifDistr("prng.unifDistr.state", std::ios::binary);
  distrStateFile_unifDistr >> unifDistr0to1_;
  distrStateFile_unifDistr.close();

  std::ifstream distrStateFile_recEvent("prng.recEvent.state", std::ios::binary);
  distrStateFile_recEvent >> recEventBernDistr_;
  distrStateFile_recEvent.close();

  std::ifstream distrStateFile_recPair("prng.recPair.state", std::ios::binary);
  distrStateFile_recPair >> recPairUnifDistr_;
  distrStateFile_recPair.close();

  std::ifstream distrStateFile_recBreakpoint("prng.recBreakpoint.state", std::ios::binary);
  distrStateFile_recBreakpoint >> recBreakpointUnifDistr_;
  distrStateFile_recBreakpoint.close();
}
*/
