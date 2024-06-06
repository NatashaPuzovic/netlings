#include "Params.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

Params::Params(
    const std::string& passedInParamsScalarFileName,
    const std::string& passedInParamsVectorFileName)
    : inParamsScalarFileName(passedInParamsScalarFileName),
    inParamsVectorFileName(passedInParamsVectorFileName),
    MUTABLE_NOISE(0),
    MUTABLE_S_BASAL(0),
    MUTABLE_NET(0),
    SEL(0),
    CONST_SEL(0),
    S_MAX(0),
    NUM_TIMESTEPS(0),
    NUM_NODES(0),
    POP_SIZE(0),
    NUM_GENERATIONS(0),
    NOISE_EXT(0),
    MUTATION_RATE_NOISE_PER_GENE(0),
    MUTATION_DISTR_NOISE_LOWERBOUND(0),
    MUTATION_DISTR_NOISE_UPPERBOUND(0),
    MUTATION_RATE_SBASAL_PER_GENE(0),
    MUTATION_DISTR_SBASAL_LOWERBOUND(0),
    MUTATION_DISTR_SBASAL_UPPERBOUND(0),
    MUTATION_RATE_NET_PER_LINK(0),
    MUTATION_DISTR_NET_LOWERBOUND(0),
    MUTATION_DISTR_NET_UPPERBOUND(0),
    REC_RATE(0),
    NUM_OUT_GENS(0),
    PRNG_SEED(0),
    SELPRESS(std::vector<double>{}),
    S_OPT(std::vector<double>{}),
    S_OPT_DELTA(std::vector<double>{}),
    NUM_SELECTEDGENES(0)
{
    // Scalar params: open connection to in file
    std::ifstream inParamsScalarFile;
    inParamsScalarFile.open(inParamsScalarFileName.c_str());

    if (!inParamsScalarFile.is_open()) {
        std::cout << "Error: Failed to open parameters (scalar) file: " << inParamsScalarFileName << std::endl;
    }
    else {

        // read params into a map
        std::string paramName;
        double paramValue;
        std::map <std::string, double> paramMap;
        while (inParamsScalarFile >> paramName >> paramValue)
        {
            paramMap[paramName] = paramValue;
        }

        // assign read values to members
        MUTABLE_NOISE = static_cast<bool>(paramMap["MUTABLE_NOISE"]);
        MUTABLE_S_BASAL = static_cast<bool>(paramMap["MUTABLE_S_BASAL"]);
        MUTABLE_NET = static_cast<bool>(paramMap["MUTABLE_NET"]);
        SEL = static_cast<bool>(paramMap["SEL"]);
        CONST_SEL = static_cast<bool>(paramMap["CONST_SEL"]);
        S_MAX = static_cast<int>(paramMap["S_MAX"]);
        NUM_TIMESTEPS = static_cast<size_t>(paramMap["NUM_TIMESTEPS"]);
        NUM_NODES = static_cast<size_t>(paramMap["NUM_NODES"]);
        POP_SIZE = static_cast<size_t>(paramMap["POP_SIZE"]);
        NUM_GENERATIONS = static_cast<size_t>(paramMap["NUM_GENERATIONS"]);
        NOISE_EXT = paramMap["NOISE_EXT"];
        MUTATION_RATE_NOISE_PER_GENE = paramMap["MUTATION_RATE_NOISE_PER_GENE"];
        MUTATION_DISTR_NOISE_LOWERBOUND = paramMap["MUTATION_DISTR_NOISE_LOWERBOUND"];
        MUTATION_DISTR_NOISE_UPPERBOUND = paramMap["MUTATION_DISTR_NOISE_UPPERBOUND"];
        MUTATION_RATE_SBASAL_PER_GENE = paramMap["MUTATION_RATE_SBASAL_PER_GENE"];
        MUTATION_DISTR_SBASAL_LOWERBOUND = paramMap["MUTATION_DISTR_SBASAL_LOWERBOUND"];
        MUTATION_DISTR_SBASAL_UPPERBOUND = paramMap["MUTATION_DISTR_SBASAL_UPPERBOUND"];
        MUTATION_RATE_NET_PER_LINK = paramMap["MUTATION_RATE_NET_PER_LINK"];
        MUTATION_DISTR_NET_LOWERBOUND = paramMap["MUTATION_DISTR_NET_LOWERBOUND"];
        MUTATION_DISTR_NET_UPPERBOUND = paramMap["MUTATION_DISTR_NET_UPPERBOUND"];
        REC_RATE = paramMap["REC_RATE"];
        NUM_OUT_GENS = static_cast<int>(paramMap["NUM_OUT_GENS"]);
        PRNG_SEED = static_cast<unsigned int>(paramMap["PRNG_SEED"]);
    }

    // Vector params: open connection to in file
    std::ifstream inParamsVectorFile;
    inParamsVectorFile.open(inParamsVectorFileName.c_str());

    if (!inParamsVectorFile.is_open()) {
        std::cout << "Error: Failed to open parameters (vector) file: " << inParamsVectorFileName << std::endl;
    }
    else {

        SELPRESS.resize(NUM_NODES);
        S_OPT.resize(NUM_NODES);
        S_OPT_DELTA.resize(NUM_NODES);

        // assign read values
        inParamsVectorFile.ignore(256, ' '); //ignore first element, i.e. the name of the param
        for (size_t j = 0; j < NUM_NODES; j++)
        {
            inParamsVectorFile >> SELPRESS[j];
        }
        inParamsVectorFile.ignore(256, ' '); //ignore first element, i.e. the name of the param
        for (size_t j = 0; j < NUM_NODES; j++)
        {
            inParamsVectorFile >> S_OPT[j];
        }
        inParamsVectorFile.ignore(256, ' '); //ignore first element, i.e. the name of the param
        for (size_t j = 0; j < NUM_NODES; j++)
        {
            inParamsVectorFile >> S_OPT_DELTA[j];
        }

        // close connection to in file
        inParamsVectorFile.close();
    }
    // determine number of selected genes
    for (size_t j = 0; j < NUM_NODES; j++)
    {
        if (SELPRESS[j] != 0)
        {
            NUM_SELECTEDGENES++;
        }
    }
}

void Params::printAllParams() const
{
    std::cout << std::endl << "Parameters (scalar) read from file " << inParamsScalarFileName << ":" << std::endl;
    std::cout << "MUTABLE_NOISE: " << MUTABLE_NOISE << std::endl;
    std::cout << "MUTABLE_S_BASAL: " << MUTABLE_S_BASAL << std::endl;
    std::cout << "MUTABLE_NET: " << MUTABLE_NET << std::endl;
    std::cout << "SEL: " << SEL << std::endl;
    std::cout << "CONST_SEL: " << CONST_SEL << std::endl;
    std::cout << "S_MAX: " << S_MAX << std::endl;
    std::cout << "NUM_TIMESTEPS: " << NUM_TIMESTEPS << std::endl;
    std::cout << "NUM_NODES: " << NUM_NODES << std::endl;
    std::cout << "POP_SIZE: " << POP_SIZE << std::endl;
    std::cout << "NUM_GENERATIONS: " << NUM_GENERATIONS << std::endl;
    std::cout << "NOISE_EXT: " << NOISE_EXT << std::endl;
    std::cout << "MUTATION_RATE_NOISE_PER_GENE: " << MUTATION_RATE_NOISE_PER_GENE << std::endl;
    std::cout << "MUTATION_DISTR_NOISE_LOWERBOUND: " << MUTATION_DISTR_NOISE_LOWERBOUND << std::endl;
    std::cout << "MUTATION_DISTR_NOISE_UPPERBOUND: " << MUTATION_DISTR_NOISE_UPPERBOUND << std::endl;
    std::cout << "MUTATION_RATE_SBASAL_PER_GENE: " << MUTATION_RATE_SBASAL_PER_GENE << std::endl;
    std::cout << "MUTATION_DISTR_SBASAL_LOWERBOUND: " << MUTATION_DISTR_SBASAL_LOWERBOUND << std::endl;
    std::cout << "MUTATION_DISTR_SBASAL_UPPERBOUND: " << MUTATION_DISTR_SBASAL_UPPERBOUND << std::endl;
    std::cout << "MUTATION_RATE_NET_PER_LINK: " << MUTATION_RATE_NET_PER_LINK << std::endl;
    std::cout << "MUTATION_DISTR_NET_LOWERBOUND: " << MUTATION_DISTR_NET_LOWERBOUND << std::endl;
    std::cout << "MUTATION_DISTR_NET_UPPERBOUND: " << MUTATION_DISTR_NET_UPPERBOUND << std::endl;
    std::cout << "REC_RATE: " << REC_RATE << std::endl;
    std::cout << "NUM_OUT_GENS: " << NUM_OUT_GENS << std::endl;
    std::cout << "PRNG_SEED: " << PRNG_SEED << std::endl;
    std::cout << "NUM_SELECTEDGENES: " << NUM_SELECTEDGENES << std::endl;

    std::cout << std::endl << "Parameters (vector) read from file " << inParamsVectorFileName << ":" << std::endl;

    std::cout << "SELPRESS: ";
    for (auto& rho : SELPRESS)
    {
        std::cout << rho << " ";
    }
    std::cout << std::endl;

    std::cout << "S_OPT: ";
    for (auto& s : S_OPT)
    {
        std::cout << s << " ";
    }
    std::cout << std::endl;

    std::cout << "S_OPT_DELTA: ";
    for (auto& s : S_OPT_DELTA)
    {
        std::cout << s << " ";
    }
    std::cout << std::endl;
}
