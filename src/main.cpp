#include <iostream>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <random>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <chrono>
#include <bitset>

#include "Netling.h"
#include "Params.h"
#include "Population.h"

////// FUNCTIONS

void outputPRNGState(std::mt19937& PRNG_ENG)
{
    std::ofstream PRNGStateFile;
    PRNGStateFile.open("prng.engine.state", std::ios::binary);
    PRNGStateFile << PRNG_ENG;
    PRNGStateFile.close();
}

void contRun_readPRNGState(std::mt19937& PRNG_ENG)
{
    std::ifstream prevPRNGStateFile("prng.engine.state", std::ios::binary);
    prevPRNGStateFile >> PRNG_ENG;
    prevPRNGStateFile.close();
}

void printUsage()
{
    std::cout << "USAGE" << std::endl <<

        "Usage for evolutionary simulations:" << std::endl <<
        "$ netlings paramScalarFile paramVectorFile populationFile outFile" << std::endl <<
        " Input (3):" << std::endl <<
        "   paramScalarFile    parameters (scalar) for the simulations, including the PRNG seed." << std::endl <<
        "   paramVectorFile    parameters (vector) for the simulations." << std::endl <<
        "   populationFile     all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)" << std::endl <<
        "                      of the first generation. The genotypes specified in the firstGenPopFile may be realized into the phenotype or not, " << std::endl <<
        "                      it does not matter because they will be realized again in any case." << std::endl <<
        " Output (5):" << std::endl <<
        "   outFile.hist       summarized phenotypes (mean phenotype, variance of phenotypes and mean fitness) per generation." << std::endl <<
        "   outFile.noises     genotypes 1: mean noise genotypes per gene per generation." << std::endl <<
        "   outFile.sbasals    genotypes 2: mean sbasal genotypes per gene per generation." << std::endl <<
        "   outFile.nets       genotypes 3: regulatory networks of all individuals in a generation." << std::endl <<
        "   outFile.pop        all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)" << std::endl <<
        "                      of the last generation." << std::endl <<
        "   prng.*             PRNG engine & distribution states, in case the run is interrupted." << std::endl << std::endl <<

        "Usage for realization of genotypes:" << std::endl <<
        "$ netlings paramScalarFile paramVectorFile populationFile -realize" << std::endl <<
        " Input (3):" << std::endl <<
        "   paramScalarFile            parameters (scalar) for the simulations, including the PRNG seed." << std::endl <<
        "   paramVectorFile            parameters (vector) for the simulations." << std::endl <<
        "   populationFile             all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)" << std::endl <<
        "                              of one generation. The genotypes specified in the populationFile may be realized into the phenotype or not, " << std::endl <<
        "                              it does not matter because they will be realized again in any case." << std::endl <<
        " Output (1):" << std::endl <<
        "   populationFile.realized    all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)" << std::endl <<
        "                              of one generation, with the phenotypes realized." << std::endl << std::endl <<

        "OPTIONS" << std::endl <<
        "   --help       print help" << std::endl;
}

//////  MAIN
int main(int argc, char* argv[])
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~                                             ~" << std::endl;
    std::cout << "~           netlings version 1.1.1            ~" << std::endl;
    std::cout << "~                                             ~" << std::endl;
    std::cout << "~     Evolve networks and gene expression     ~" << std::endl;
    std::cout << "~                                             ~" << std::endl;
    std::cout << "~  Author: Nataša Puzović, Julien Y. Dutheil  ~" << std::endl;
    std::cout << "~  Email:  npuzovic@protonmail.com            ~" << std::endl;
    std::cout << "~                                             ~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    // if insufficient input is provided
    if (argc != 5)
    {
        std::cout << std::endl << "Invalid argument input. Printing usage." << std::endl << std::endl;
        printUsage();
        return 1;
    }

    // sort commandline arguments
    std::string paramsScalarFileName = argv[1];
    std::string paramsVectorFileName = argv[2];
    std::string populationFileName = argv[3];
    std::string outFileNameTemplate = argv[4];

    // REALIZATION OF GENOTYPES
    if (outFileNameTemplate == "-realize")
    {
        // rename outfile so the realized population can be written to outfile.realized.
        outFileNameTemplate = populationFileName;

        // read parameters and print read parameters to terminal
        Params prms(paramsScalarFileName, paramsVectorFileName);
        prms.printAllParams();

        // construct PRNG initialized with the seed passed as parameter (or read from a previous run).
        // If continuation of a previous run is chosen, the state will be reassigned ~50 lines below.
        std::mt19937 PRNG_ENG(prms.PRNG_SEED);

        // get steady state expression levels from a realization of a non-noisy network
        Netling::sMax = prms.S_MAX;
        Netling::numTimesteps = prms.NUM_TIMESTEPS;
        Netling::numNodes = prms.NUM_NODES;
        Netling::noiseExt = prms.NOISE_EXT;
        std::vector<std::vector<double>> emptyRmtrx(prms.NUM_NODES, std::vector<double>(prms.NUM_NODES));

        // set optimal expression levels from the read params
        std::vector<double> sOpt(prms.NUM_NODES, 0);
        for (size_t i = 0; i < prms.NUM_NODES; i++)
        {
            sOpt[i] = prms.S_OPT[i] + prms.S_OPT_DELTA[i];
            if (sOpt[i] < 0)
            {
                std::cout << "WARNING: Optimal expression level of a gene is below 0." << std::endl;
            }
            if (sOpt[i] > prms.S_MAX)
            {
                std::cout << "WARNING: Optimal expression level of a gene is above the maximum expression level." << std::endl;
            }
        }
        std::cout << std::endl << "Optimal expression levels (sOpt):" << std::endl;
        for (auto& s : sOpt) { std::cout << s << " "; }
        std::cout << std::endl;

        // generate population with default netling genotypes and configurations
        Netling netling0(std::vector<double>(prms.NUM_NODES, 0),
            std::vector<double>(prms.NUM_NODES, 0),
            emptyRmtrx);
        Population pop(prms.POP_SIZE, prms.NUM_NODES, netling0,
            prms.MUTATION_RATE_NOISE_PER_GENE, prms.MUTATION_DISTR_NOISE_LOWERBOUND, prms.MUTATION_DISTR_NOISE_UPPERBOUND,
            prms.MUTATION_RATE_SBASAL_PER_GENE, prms.MUTATION_DISTR_SBASAL_LOWERBOUND, prms.MUTATION_DISTR_SBASAL_UPPERBOUND,
            prms.MUTATION_RATE_NET_PER_LINK, prms.MUTATION_DISTR_NET_LOWERBOUND, prms.MUTATION_DISTR_NET_UPPERBOUND,
            prms.REC_RATE, outFileNameTemplate);

        // change the genotypes in default population to the input values
        pop.readFirstGenEntirePopulation(populationFileName);

        // realize and calculate fitness
        for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
        {
            pop.netlings[cell].realizePhenotype(PRNG_ENG);
            pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
        }

        // output the population
        pop.outputEntirePopulation(0);
        pop.closeOutputFiles();

        // rename the .pop output file into .realized file
        std::rename((outFileNameTemplate + ".pop").c_str(), (outFileNameTemplate + ".realized").c_str());
        std::cout << "Realized genotypes written to " << outFileNameTemplate << ".realized." << std::endl;
        std::remove((outFileNameTemplate + ".hist").c_str());
        std::remove((outFileNameTemplate + ".noises").c_str());
        std::remove((outFileNameTemplate + ".sbasals").c_str());
        std::remove((outFileNameTemplate + ".nets").c_str());
        return 0;

    }
    else

        // EVOLUTIONARY SIMULATIONS
    {
        // start clock to time execution
        auto simStart = std::chrono::steady_clock::now();

        // read parameters and print read parameters to terminal
        Params prms(paramsScalarFileName, paramsVectorFileName);
        prms.printAllParams();

        // bit code to identify the selection scenario to simulate
        std::bitset<5> scenarioCase; // 5-bit code: |mut.noise|mut.S|mut.net|sel|constSel|
        if (prms.CONST_SEL) { scenarioCase.set(0, 1); }
        if (prms.SEL) { scenarioCase.set(1, 1); }  // selection or not?
        if (prms.MUTABLE_NET) { scenarioCase.set(2, 1); }  // mutable net topology or not?
        if (prms.MUTABLE_S_BASAL) { scenarioCase.set(3, 1); }  // mutable basal expr. level or not?
        if (prms.MUTABLE_NOISE) { scenarioCase.set(4, 1); }  // mutable intr. expr. noise or not?

        // all scenario codes:
        std::bitset<5> constSelection_mutNoise("10011");
        std::bitset<5> constSelection_mutNoiseSBasal("11011");
        std::bitset<5> constSelection_mutNoiseNet("10111");
        std::bitset<5> constSelection_mutSBasal("01011");
        std::bitset<5> constSelection_mutNet("00111");

        std::bitset<5> fluctSelection_mutNoise("10010");
        std::bitset<5> fluctSelection_mutNoiseSBasal("11010");
        std::bitset<5> fluctSelection_mutNoiseNet("10110");
        std::bitset<5> fluctSelection_mutSBasal("01010");
        std::bitset<5> fluctSelection_mutNet("00110");

        std::bitset<5> neutrality_mutNoise("10000");
        std::bitset<5> neutrality_mutNoiseSBasal("11000");
        std::bitset<5> neutrality_mutNoiseNet("10100");
        std::bitset<5> neutrality_mutSBasal("01000");
        std::bitset<5> neutrality_mutNet("00100");


        std::cout << std::endl << "scenarioCase " << scenarioCase << ":" << std::endl;

        // construct PRNG initialized with the seed passed as parameter (or read from a previous run).
        // If continuation of a previous run is chosen, the state will be reassigned ~50 lines below.
        std::mt19937 PRNG_ENG(prms.PRNG_SEED);

        // get steady state expression levels from a realization of a non-noisy network
        Netling::sMax = prms.S_MAX;
        Netling::numTimesteps = prms.NUM_TIMESTEPS;
        Netling::numNodes = prms.NUM_NODES;
        Netling::noiseExt = prms.NOISE_EXT;
        std::vector<std::vector<double>> emptyRmtrx(prms.NUM_NODES, std::vector<double>(prms.NUM_NODES));

        // set optimal expression levels from the read params
        std::vector<double> sOpt(prms.NUM_NODES, 0);
        for (size_t i = 0; i < prms.NUM_NODES; i++)
        {
            sOpt[i] = prms.S_OPT[i] + prms.S_OPT_DELTA[i];
            if (sOpt[i] < 0)
            {
                std::cout << "WARNING: Optimal expression level of a gene is below 0." << std::endl;
            }
            if (sOpt[i] > prms.S_MAX)
            {
                std::cout << "WARNING: Optimal expression level of a gene is above the maximum expression level." << std::endl;
            }
        }
        std::cout << std::endl << "Optimal expression levels (sOpt):" << std::endl;
        for (auto& s : sOpt) { std::cout << s << " "; }
        std::cout << std::endl;

        // generate population with default netling genotypes and configurations
        Netling netling0(std::vector<double>(prms.NUM_NODES, 0),
            std::vector<double>(prms.NUM_NODES, 0),
            emptyRmtrx);
        Population pop(prms.POP_SIZE, prms.NUM_NODES, netling0,
            prms.MUTATION_RATE_NOISE_PER_GENE, prms.MUTATION_DISTR_NOISE_LOWERBOUND, prms.MUTATION_DISTR_NOISE_UPPERBOUND,
            prms.MUTATION_RATE_SBASAL_PER_GENE, prms.MUTATION_DISTR_SBASAL_LOWERBOUND, prms.MUTATION_DISTR_SBASAL_UPPERBOUND,
            prms.MUTATION_RATE_NET_PER_LINK, prms.MUTATION_DISTR_NET_LOWERBOUND, prms.MUTATION_DISTR_NET_UPPERBOUND,
            prms.REC_RATE, outFileNameTemplate);

        // change the genotypes in default population to the input values
        pop.readFirstGenEntirePopulation(populationFileName);

        // calculate generation intervals to output
        size_t genForOut = 0;
        size_t genForOutInterval = static_cast<size_t>(prms.NUM_GENERATIONS / prms.NUM_OUT_GENS);

        ///// SCENARIO (NOISE): const selection, mutable noise
        if (scenarioCase == constSelection_mutNoise)
        {
            std::cout << std::endl << "Evolutionary scenario: constant selection with mutable noise." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // realize phenotypes & calc fitness for the entire population
                for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                {
                    pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                }
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
            }

            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE): fluctuating selection, mutable noise
        if (scenarioCase == fluctSelection_mutNoise)
        {
            std::cout << std::endl << "Evolutionary scenario: fluctuating selection with mutable noise." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // set up sOpt for the second environment
            std::vector<double> sOpt2(prms.NUM_NODES, 0);
            for (size_t i = 0; i < prms.NUM_NODES; i++)
            {
                sOpt2[i] = prms.S_OPT[i] - prms.S_OPT_DELTA[i];
                if (sOpt2[i] < 0)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is below 0." << std::endl;
                }
                if (sOpt2[i] > prms.S_MAX)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is above the maximum expression level." << std::endl;
                }
            }
            std::cout << std::endl << "Optimal expression levels, second environmment (sOpt2):" << std::endl;
            for (auto& s : sOpt2) { std::cout << s << " "; }
            std::cout << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                if (gen % 2 == 0) {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                    }
                }
                else {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt2, prms.SELPRESS);
                    }
                }
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE): neutrality, mutable noise
        if (scenarioCase == neutrality_mutNoise)
        {
            std::cout << std::endl << "Evolutionary scenario: neutral evolution with mutable noise." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under neutrality..." << std::endl;

            // start sim from beginning
              // set fitness of all individuals to 0.5, it will not change throughout evol
            for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
            {
                pop.netlings[cell].fitness = 0.5;
            }
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // realize phenotypes for the entire population only in output generations
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    }
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE + SBASAL): const selection, mutable noise + SBasal
        if (scenarioCase == constSelection_mutNoiseSBasal)
        {
            std::cout << std::endl << "Evolutionary scenario: constant selection with mutable noise and basal expression level." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // realize phenotypes & calc fitness for the entire population
                for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                {
                    pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                }
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);
                    pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }

                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                pop.mutateSBasalGenotypes(PRNG_ENG);
                pop.recombineSBasalGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE + SBASAL): fluctuating selection, mutable noise + SBasal
        if (scenarioCase == fluctSelection_mutNoiseSBasal)
        {
            std::cout << std::endl << "Evolutionary scenario: fluctuating selection with mutable noise and basal expression level." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // set up sOpt for the second environment
            std::vector<double> sOpt2(prms.NUM_NODES, 0);
            for (size_t i = 0; i < prms.NUM_NODES; i++)
            {
                sOpt2[i] = prms.S_OPT[i] - prms.S_OPT_DELTA[i];
                if (sOpt2[i] < 0)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is below 0." << std::endl;
                }
                if (sOpt2[i] > prms.S_MAX)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is above the maximum expression level." << std::endl;
                }
            }
            std::cout << std::endl << "Optimal expression levels, second environmment (sOpt2):" << std::endl;
            for (auto& s : sOpt2) { std::cout << s << " "; }
            std::cout << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {

                if (gen % 2 == 0) {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                    }
                }
                else {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt2, prms.SELPRESS);
                    }
                }

                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);
                    pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }

                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                pop.mutateSBasalGenotypes(PRNG_ENG);
                pop.recombineSBasalGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE + SBASAL): neutrality, mutable noise + SBasal
        if (scenarioCase == neutrality_mutNoiseSBasal)
        {
            std::cout << std::endl << "Evolutionary scenario: neutral evolution with mutable noise and basal expression level." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under neutrality..." << std::endl;

            // start sim from beginning
              // set fitness of all individuals to 0.5, it will not change throughout evol
            for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
            {
                pop.netlings[cell].fitness = 0.5;
            }
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // realize phenotypes for the entire population only in output generations
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    }
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);
                    pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                pop.mutateSBasalGenotypes(PRNG_ENG);
                pop.recombineSBasalGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE + NET): const selection, mutable noise + net topology
        if (scenarioCase == constSelection_mutNoiseNet)
        {
            std::cout << std::endl << "Evolutionary scenario: constant selection with mutable noise and network topology links." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // realize phenotypes & calc fitness for the entire population
                for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                {
                    pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                }
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                pop.mutateNetExistingRegLinks(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE + NET): fluctuating selection, mutable noise + net topology
        if (scenarioCase == fluctSelection_mutNoiseNet)
        {
            std::cout << std::endl << "Evolutionary scenario: fluctuating selection with mutable noise and network topology links." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // set up sOpt for the second environment
            std::vector<double> sOpt2(prms.NUM_NODES, 0);
            for (size_t i = 0; i < prms.NUM_NODES; i++)
            {
                sOpt2[i] = prms.S_OPT[i] - prms.S_OPT_DELTA[i];
                if (sOpt2[i] < 0)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is below 0." << std::endl;
                }
                if (sOpt2[i] > prms.S_MAX)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is above the maximum expression level." << std::endl;
                }
            }
            std::cout << std::endl << "Optimal expression levels, second environmment (sOpt2):" << std::endl;
            for (auto& s : sOpt2) { std::cout << s << " "; }
            std::cout << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                if (gen % 2 == 0) {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                    }
                }
                else {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt2, prms.SELPRESS);
                    }
                }

                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                pop.mutateNetExistingRegLinks(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (NOISE + NET): neutrality, mutable noise + net topology
        if (scenarioCase == neutrality_mutNoiseNet)
        {
            std::cout << std::endl << "Evolutionary scenario: neutral evolution with mutable noise and network topology links." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under neutrality..." << std::endl;

            // set fitness of all individuals to 0.5, it will not change throughout evol
            for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
            {
                pop.netlings[cell].fitness = 0.5;
            }
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // realize phenotypes for the entire population only in output generations
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    }
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNoiseGenotypes(PRNG_ENG);
                pop.recombineNoiseGenotypes(PRNG_ENG);
                pop.mutateNetExistingRegLinks(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        // here
          ///// SCENARIO (SBASAL): const selection, mutable SBasal
        if (scenarioCase == constSelection_mutSBasal)
        {
            std::cout << std::endl << "Evolutionary scenario: constant selection with mutable basal expression level." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // realize phenotypes & calc fitness for the entire population
                for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                {
                    pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                }
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }

                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateSBasalGenotypes(PRNG_ENG);
                pop.recombineSBasalGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (SBASAL): fluctuating selection, mutable SBasal
        if (scenarioCase == fluctSelection_mutSBasal)
        {
            std::cout << std::endl << "Evolutionary scenario: fluctuating selection with mutable basal expression level." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // set up sOpt for the second environment
            std::vector<double> sOpt2(prms.NUM_NODES, 0);
            for (size_t i = 0; i < prms.NUM_NODES; i++)
            {
                sOpt2[i] = prms.S_OPT[i] - prms.S_OPT_DELTA[i];
                if (sOpt2[i] < 0)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is below 0." << std::endl;
                }
                if (sOpt2[i] > prms.S_MAX)
                {
                    std::cout << "WARNING: Optimal expression level of a gene is above the maximum expression level." << std::endl;
                }
            }
            std::cout << std::endl << "Optimal expression levels, second environmment (sOpt2):" << std::endl;
            for (auto& s : sOpt2) { std::cout << s << " "; }
            std::cout << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {

                if (gen % 2 == 0) {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                    }
                }
                else {
                    // realize phenotypes & calc fitness for the entire population
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                        pop.netlings[cell].calculateFitness(sOpt2, prms.SELPRESS);
                    }
                }

                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedNoiseGenotypes(gen);
                    pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }

                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateSBasalGenotypes(PRNG_ENG);
                pop.recombineSBasalGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        ///// SCENARIO (SBASAL): neutrality, mutable SBasal
        if (scenarioCase == neutrality_mutSBasal)
        {
            std::cout << std::endl << "Evolutionary scenario: neutral evolution with mutable basal expression level." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under neutrality..." << std::endl;

            // start sim from beginning
              // set fitness of all individuals to 0.5, it will not change throughout evol
            for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
            {
                pop.netlings[cell].fitness = 0.5;
            }
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // realize phenotypes for the entire population only in output generations
                    for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                    {
                        pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    }
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }
                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateSBasalGenotypes(PRNG_ENG);
                pop.recombineSBasalGenotypes(PRNG_ENG);
            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }


        ///// SCENARIO (NET): const selection, mutable net
        if (scenarioCase == constSelection_mutNet)
        {
            std::cout << std::endl << "Evolutionary scenario: constant selection with mutable network topology." << std::endl;
            std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

            // start sim from beginning
            for (size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
            {
                // realize phenotypes & calc fitness for the entire population
                for (size_t cell = 0; cell < prms.POP_SIZE; cell++)
                {
                    pop.netlings[cell].realizePhenotype(PRNG_ENG);
                    pop.netlings[cell].calculateFitness(sOpt, prms.SELPRESS);
                }
                // write out if it's an intermediate generations
                if (gen == genForOut)
                {
                    // write PVF (phenotypes, variance of phenotypes, fitness)
                    pop.outputSummarizedPVF(gen);
                    //pop.outputNetworkTopo();
                    //pop.outputSummarizedNoiseGenotypes(gen);
                    //pop.outputSummarizedSBasalGenotypes(gen);

                    // write PRNG engine & distribution state (in case run is to be continued)
                    outputPRNGState(PRNG_ENG);
                    pop.outputDistrStates();

                    // update counter for generation output
                    if (gen == 0)
                    {
                        genForOut += genForOutInterval - 1;
                    }
                    else
                    {
                        genForOut += genForOutInterval;
                    }
                    //std::cout << "Written generation " << gen + 1 << " to output files." << std::endl;
                }
                // write out the entire population if it's the last generation
                if (gen == (prms.NUM_GENERATIONS - 1))
                {
                    pop.outputEntirePopulation(gen);
                }

                // reproduce, mutate, recombine
                pop.reproduceNetlings(PRNG_ENG);
                pop.mutateNetExistingRegLinks(PRNG_ENG);

            }
            // close connections to output files & print finish message
            pop.closeOutputFiles();
            pop.printFinishMessage(outFileNameTemplate);
        }

        // report elapsed execution time
        auto simStop = std::chrono::steady_clock::now();
        auto simElapsedTime = std::chrono::duration<double>(simStop - simStart);
        std::cout << std::endl << "Simulation completed in " << simElapsedTime.count() / 60 << " minutes." << std::endl;
    }
}
