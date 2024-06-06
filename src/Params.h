#pragma once
#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

class Params
{
public:
    const std::string inParamsScalarFileName;
    const std::string inParamsVectorFileName;
    bool MUTABLE_NOISE;
    bool MUTABLE_S_BASAL;
    bool MUTABLE_NET;
    bool SEL;
    bool CONST_SEL;
    int S_MAX;
    size_t NUM_TIMESTEPS;
    size_t NUM_NODES;
    size_t POP_SIZE;
    size_t NUM_GENERATIONS;
    double NOISE_EXT;
    double MUTATION_RATE_NOISE_PER_GENE;
    double MUTATION_DISTR_NOISE_LOWERBOUND;
    double MUTATION_DISTR_NOISE_UPPERBOUND;
    double MUTATION_RATE_SBASAL_PER_GENE;
    double MUTATION_DISTR_SBASAL_LOWERBOUND;
    double MUTATION_DISTR_SBASAL_UPPERBOUND;
    double MUTATION_RATE_NET_PER_LINK;
    double MUTATION_DISTR_NET_LOWERBOUND;
    double MUTATION_DISTR_NET_UPPERBOUND;
    double REC_RATE;
    int NUM_OUT_GENS;
    unsigned int PRNG_SEED;
    std::vector<double> SELPRESS;
    std::vector<double> S_OPT;
    std::vector<double> S_OPT_DELTA;
    int NUM_SELECTEDGENES;

public:
    Params(const std::string& passedInParamsScalarFileName,
        const std::string& passedInParamsVectorFileName);
    void printAllParams() const;
};
