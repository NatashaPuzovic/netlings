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

#include "Params.h"
#include "Netling.h"
#include "Population.h"


////// FUNCTIONS
// read matrix
std::vector<std::vector<double>> readMatrix(size_t nrow, size_t ncol, std::string infile_name)
{
  std::vector<std::vector<double>> Mat( nrow , std::vector<double> (ncol));
  std::ifstream textfile(infile_name.c_str());
  for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            textfile >> Mat[i][j];
        }
  }
  return Mat;
}

// print matrix
void PrintMat(std::vector<std::vector<double>> mat, size_t rows, size_t cols)
{
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      std::cout << mat[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void writePRNGState(std::mt19937 &PRNG_ENG)
{
  std::ofstream PRNGStateFile;
  PRNGStateFile.open("prng.engine.state", std::ios::binary);
  PRNGStateFile << PRNG_ENG;
  PRNGStateFile.close();
}

void contRun_readPRNGState(std::mt19937 &PRNG_ENG)
{
  std::ifstream prevPRNGStateFile("prng.engine.state", std::ios::binary);
  prevPRNGStateFile >> PRNG_ENG;
  prevPRNGStateFile.close();
}

void writeGPF(size_t gen, Params prms, Population &pop, std::ofstream &outGPFfile)
{
  for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
  {
    outGPFfile << gen + 1 << " " << cell + 1 << " ";
    for(size_t j = 0; j < prms.NUM_NODES; j++)
    {
      outGPFfile << pop.netlings[cell].genotNoise[j] << " ";
    }
    for(size_t j = 0; j < prms.NUM_NODES; j++)
    {
      outGPFfile << pop.netlings[cell].phenot[j] << " ";
    }
    outGPFfile << pop.netlings[cell].fitness;
    outGPFfile << std::endl;
  }
}

void writeSummarizedGPF(size_t gen, Params prms, Population &pop, std::ofstream &outSummarizedGPFfile)
{
  double sumG_j = 0;
  double sumP_j = 0;
  double meanG_j = 0;
  double varP_j = 0;
  double meanF = 0;
  std::vector<double> meanP (prms.NUM_NODES, 0);

  outSummarizedGPFfile << gen + 1 << " ";
  // mean genotypes
  for(size_t j = 0; j < prms.NUM_NODES; j++)
  {
    sumG_j = 0;
    for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
    {
      sumG_j += pop.netlings[cell].genotNoise[j];
    }
    meanG_j = sumG_j/static_cast<double>(prms.POP_SIZE);
    outSummarizedGPFfile << meanG_j << " ";
  }
  // mean phenotypes
  for(size_t j = 0; j < prms.NUM_NODES; j++)
  {
    sumP_j = 0;
    for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
    {
      sumP_j += pop.netlings[cell].phenot[j];
    }
    meanP[j] = sumP_j/static_cast<double>(prms.POP_SIZE);
    outSummarizedGPFfile << meanP[j] << " ";
  }
  // variances of phenotypes
  for(size_t j = 0; j < prms.NUM_NODES; j++)
  {
    varP_j = 0;
    for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
    {
      varP_j += (pop.netlings[cell].phenot[j] - meanP[j])*(pop.netlings[cell].phenot[j] - meanP[j]);
    }
    varP_j = varP_j/static_cast<double>(prms.POP_SIZE);
    outSummarizedGPFfile << varP_j << " ";
  }
  // add mean fitness
  for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
  {
    meanF += pop.netlings[cell].fitness;
  }
  meanF = meanF/static_cast<double>(prms.POP_SIZE);
  outSummarizedGPFfile << meanF << " ";

  // finish line
  outSummarizedGPFfile << std::endl;
}

void printUsage()
{
  std::cout <<
  "enoise1.2" << std::endl <<
  "Usage: enoise networkFile paramFile outfile [-continue]" << std::endl <<
  std::endl <<
  " Input:" << std::endl <<
  "   networkFile: regulatory matrix that defines the strength of interactions" << std::endl <<
  "                of each pair of genes in the network." << std::endl <<
  "   paramFile:   parameters for the simulations, including the PRNG seed." << std::endl <<
  " Output:" << std::endl <<
  "   outfile:     genotype, phenotype and fitness values of individuals in" << std::endl <<
  "                specified generations to output." << std::endl <<
  "   prng* files: PRNG engine & distribution states, in case the run is interrupted." << std::endl <<
  std::endl <<
  "Options:" << std::endl <<
  "   -continue    continue an interrupted run" << std::endl <<
  "                this option reads the outfile and the PRNG states of a previous run" << std::endl <<
  "                and continues the simulation from the last recorded generation." << std::endl <<
  "   --help       print this help" << std::endl;
}

//////  MAIN
int main(int argc, char *argv[])
{
  auto simStart = std::chrono::steady_clock::now();

  // sort commandline arguments
  if((argc != 4) && (argc != 5))
  {
    printUsage();
    return 1;
  }

  std::string netFileName = argv[1];
  std::string paramFileName = argv[2];
  std::string outGPFfileName = argv[3];
  std::string continueMode = "no";
  std::string outSummarizedGPFfileName = argv[3];
  outSummarizedGPFfileName.append(".summary");

  if(netFileName == "--help")
  {
    printUsage();
  }

  if(argc == 5)
  {
    continueMode = argv[4];
  }

  // read parameters
  Params prms(paramFileName);
  std::cout << std::endl << "Parameters from " << paramFileName << ":" << std::endl;
  prms.printAllParams();

  std::vector<double> sBasal (prms.NUM_NODES, prms.S_BASAL_ALL);
  std::vector<double> initNoiseInt (prms.NUM_NODES, prms.STARTING_NOISE_INT_ALL);
  std::vector<double> selPressure (prms.NUM_NODES, prms.SELPRESS_NOISE_EVOL_ALL);

  // read in rmtrx
  const std::vector<std::vector<double>> rmtrx = readMatrix(prms.NUM_NODES, prms.NUM_NODES, netFileName);
  std::cout << std::endl << "Network configuration from " << netFileName << ":" << std::endl;
  PrintMat(rmtrx, prms.NUM_NODES, prms.NUM_NODES);

  // construct PRNG initialized with the seed passed as parameter (or read from a previous run).
  // If continuation of a previous run is chosen, the state will be reassigned ~50 lines below.
  std::mt19937 PRNG_ENG(prms.PRNG_SEED);

  // get sOpt from a realization of a non-noisy network
  std::vector<double> nullNoise (prms.NUM_NODES, 0);
  Netling netling0(prms.S_MAX, prms.NUM_TIMESTEPS, prms.NUM_NODES, prms.POP_SIZE, prms.NOISE_EXT,
    sBasal, rmtrx, nullNoise);
  netling0.realizePhenotype(PRNG_ENG);
  std::vector<double> sOpt (prms.NUM_NODES, 0);
  for(size_t i = 0; i < prms.NUM_NODES; i++)
  {
    sOpt[i] = netling0.phenot[i];
  }
  std::cout << std::endl << "Optimal expression levels (sOpt):" << std::endl;
  for(auto& s : sOpt)
  {
    std::cout << s << " ";
  }
  std::cout << std::endl;

  // founding genotype & population
  Netling founderNetling(prms.S_MAX, prms.NUM_TIMESTEPS, prms.NUM_NODES, prms.POP_SIZE, prms.NOISE_EXT,
    sBasal, rmtrx, initNoiseInt);
  Population pop(prms.POP_SIZE, prms.NUM_NODES, founderNetling,
    prms.MU_NOISE_PER_GENE, prms.NOISE_MUT_DISTR_MEAN, prms.NOISE_MUT_DISTR_SD, prms.REC_RATE);

  // calc. generation intervals to output
  size_t genForOut = 0;
  size_t genForOutInterval = static_cast<size_t>(prms.NUM_GENERATIONS/prms.NUM_OUT_GENS);
  //std::cout << "genfor out interval: " << genForOutInterval << std::endl;
  int lastGen = 0;

  // if continue
  if(continueMode == "-continue")
  {
    // read previous prng.engine.state & assign to current PRNG engine
    contRun_readPRNGState(PRNG_ENG);
    pop.contRun_readDistrStates();

    std::cout << "Info: Loaded PRNG engine and distribution states from prng.state files" << std::endl;

    // find number of lines in GPF & read last generation
    int lineCount = 0;
    std::string line;
    //outGPFfile.open(outGPFfileName.c_str(), std::ios_base::app)
    std::ifstream inGPFFile(outGPFfileName.c_str());
    if(inGPFFile.is_open())
    {
      while(!inGPFFile.eof())
      {
        getline(inGPFFile, line);
        lineCount++;
      }
      inGPFFile.close();
    }
    int numWrittenGens = ((lineCount - 1)/static_cast<int>(prms.POP_SIZE));
    lastGen = (numWrittenGens - 1)*static_cast<int>(genForOutInterval);
    std::cout << "Previous run stopped at generation: " << lastGen << std::endl;

    //inGPFFile.open("GPF");
    inGPFFile.open(outGPFfileName.c_str());
    int ignoreLines = (numWrittenGens - 1)*static_cast<int>(prms.POP_SIZE);
    size_t numCols = prms.NUM_NODES*2 + 3;
    std::vector<double> prevCellInfo (numCols, 0);
    // loop through unwanted info
    for(int lineNum = 0; lineNum < ignoreLines; lineNum++)
    {
      getline(inGPFFile, line);
    }
    // then read the last generation
    for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
    {
      getline(inGPFFile, line);
      std::stringstream prevCellInfoStream(line);
      for(size_t i = 0; i < numCols; i++)
      {
        prevCellInfoStream >> prevCellInfo[i];
      }
      // assign genotypes of netlings from last gen to current pop
      for(size_t j = 0; j < prms.NUM_NODES; j++)
      {
        pop.netlings[cell].genotNoise[j] = prevCellInfo[j + 2];
      }
      // assign fitnesses of netlings from last gen to current pop (phenotypes are not needed)
      pop.netlings[cell].fitness = prevCellInfo[numCols - 1];
    }
    inGPFFile.close();
  }

  // open connection to output GPF textfiles (new for new run, append for continued run)
  /* I commented the opening of the long GPF file because I'm not using it anymore.
  std::ofstream outGPFfile;
  if(continueMode == "-continue")
  {
    outGPFfile.open(outGPFfileName.c_str(), std::ios_base::app);
  }
  else
  {
    outGPFfile.open(outGPFfileName.c_str(), std::fstream::out);
  }
  */
  std::ofstream outSummarizedGPFfile;
  if(continueMode == "-continue")
  {
    outSummarizedGPFfile.open(outSummarizedGPFfileName.c_str(), std::ios_base::app);
  }
  else
  {
    outSummarizedGPFfile.open(outSummarizedGPFfileName.c_str(), std::fstream::out);
  }

  ///// SELECTION
  if(prms.SELPRESS_NOISE_EVOL_ALL == 1)
  {
    std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

    // start sim from beginning
    if(continueMode == "no")
    {
      for(size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
      {
        // realize phenotypes & calc fitness for the entire population
        for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
        {
          pop.netlings[cell].realizePhenotype(PRNG_ENG);
          pop.netlings[cell].calculateFitness(sOpt, selPressure);
        }
        // write out if it's an intermediate generations
        if(gen == genForOut)
        {
          //write GPF
          //writeGPF(gen, prms, pop, outGPFfile);
          writeSummarizedGPF(gen, prms, pop, outSummarizedGPFfile);
          std::cout << std::endl << "Written generation " << gen + 1 << " to " << outGPFfileName << "." << std::endl;

          // write PRNG engine & distribution state (in case run is to be continued)
          writePRNGState(PRNG_ENG);
          pop.writeDistrStates();

          std::cout << "Written current PRNG & distribution states to prng.state files." << std::endl;

          // update counter for generation output
          if(gen == 0)
          {
            genForOut += genForOutInterval - 1;
          }
          else
          {
            genForOut += genForOutInterval;
          }
        }
        // reproduce, mutate, recombine
        pop.reproduceNoiseGenotypes(PRNG_ENG);
        pop.mutateNoiseGenotypes(PRNG_ENG);
        pop.recombineNoiseGenotypes(PRNG_ENG);
      }

      // close connections to output files
      //outGPFfile.close();
      outSummarizedGPFfile.close();
      // report
      std::cout << std::endl << "Written genotypes/phenotypes/fitness of all generations to " << outGPFfileName << "." << std::endl;
    }

    // resume previous sim run
    else if (continueMode == "-continue")
    {
      std::cout << std::endl << "Continuing previous run. Evolving population of " << prms.POP_SIZE << " networks under selection..." << std::endl;

      genForOut = static_cast<size_t>(lastGen) - 1 + genForOutInterval;
      // reproduce, mutate, recombine
      pop.reproduceNoiseGenotypes(PRNG_ENG);
      pop.mutateNoiseGenotypes(PRNG_ENG);
      pop.recombineNoiseGenotypes(PRNG_ENG);

      for(size_t gen = lastGen; gen < prms.NUM_GENERATIONS; gen++)
      {
        // realize phenotypes & calc fitness for the entire population
        for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
        {
          pop.netlings[cell].realizePhenotype(PRNG_ENG);
          pop.netlings[cell].calculateFitness(sOpt, selPressure);
        }
        // write out if it's an intermediate generations
        if(gen == genForOut)
        {
          //write GPF
          //writeGPF(gen, prms, pop, outGPFfile);
          writeSummarizedGPF(gen, prms, pop, outSummarizedGPFfile);
          std::cout << std::endl << "Written generation " << gen + 1 << " to " << outGPFfileName << "." << std::endl;

          // write PRNG engine & distribution state (in case run is to be continued)
          writePRNGState(PRNG_ENG);
          pop.writeDistrStates();

          std::cout << "Written current PRNG & distribution states to prng.state files." << std::endl;

          // update counter for generation output
          if(gen == 0)
          {
            genForOut += genForOutInterval - 1;
          }
          else
          {
            genForOut += genForOutInterval;
          }
        }
        // reproduce, mutate, recombine
        pop.reproduceNoiseGenotypes(PRNG_ENG);
        pop.mutateNoiseGenotypes(PRNG_ENG);
        pop.recombineNoiseGenotypes(PRNG_ENG);
      }
      // close connections to output files
      //outGPFfile.close();
      outSummarizedGPFfile.close();
      // report
      std::cout << std::endl << "Written genotypes/phenotypes/fitness of all generations to " << outGPFfileName << "." << std::endl;
    }
  }

  ///// NEUTRALITY
  else if (prms.SELPRESS_NOISE_EVOL_ALL == 0)
  {
    std::cout << std::endl << "Evolving population of " << prms.POP_SIZE << " networks under neutrality..." << std::endl;

    // start sim from beginning
    if(continueMode == "no")
    {
      // set fitness of all individuals to 0.5, it will not change throughout evol
      for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
      {
        pop.netlings[cell].fitness = 0.5;
      }
      for(size_t gen = 0; gen < prms.NUM_GENERATIONS; gen++)
      {
        // write out if it's an intermediate generations
        if(gen == genForOut)
        {
          // realize phenotypes for the entire population only in output generations
          for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
          {
            pop.netlings[cell].realizePhenotype(PRNG_ENG);
          }
          //write GPF
          //writeGPF(gen, prms, pop, outGPFfile);
          writeSummarizedGPF(gen, prms, pop, outSummarizedGPFfile);
          std::cout << std::endl << "Written generation " << gen + 1 << " to " << outGPFfileName << "." << std::endl;

          // write PRNG engine & distribution state (in case run is to be continued)
          writePRNGState(PRNG_ENG);
          pop.writeDistrStates();

          std::cout << "Written current PRNG & distribution states to prng.state files." << std::endl;

          // update counter for generation output
          if(gen == 0)
          {
            genForOut += genForOutInterval - 1;
          }
          else
          {
            genForOut += genForOutInterval;
          }
        }
        // reproduce, mutate, recombine
        pop.reproduceNoiseGenotypes(PRNG_ENG);
        pop.mutateNoiseGenotypes(PRNG_ENG);
        pop.recombineNoiseGenotypes(PRNG_ENG);
      }
      // close connections to output files
      //outGPFfile.close();
      outSummarizedGPFfile.close();
      // report
      std::cout << std::endl << "Written genotypes/phenotypes/fitness of all generations to " << outGPFfileName << "." << std::endl;
    }

    // resume previous sim run
    else if (continueMode == "-continue")
    {
      std::cout << std::endl << "Continuing previous run. Evolving population of " << prms.POP_SIZE << " networks under neutrality..." << std::endl;

      // set fitness of all individuals to 0.5, it will not change throughout evol
      for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
      {
        pop.netlings[cell].fitness = 0.5;
      }
      genForOut = static_cast<size_t>(lastGen) - 1 + genForOutInterval;
      // reproduce, mutate, recombine
      pop.reproduceNoiseGenotypes(PRNG_ENG);
      pop.mutateNoiseGenotypes(PRNG_ENG);
      pop.recombineNoiseGenotypes(PRNG_ENG);

      for(size_t gen = lastGen; gen < prms.NUM_GENERATIONS; gen++)
      {
        // write out if it's an intermediate generations
        if(gen == genForOut)
        {
          // realize phenotypes for the entire population only in output generations
          for(size_t cell = 0; cell < prms.POP_SIZE; cell++)
          {
            pop.netlings[cell].realizePhenotype(PRNG_ENG);
          }
          //write GPF
          //writeGPF(gen, prms, pop, outGPFfile);
          writeSummarizedGPF(gen, prms, pop, outSummarizedGPFfile);
          std::cout << std::endl << "Written generation " << gen + 1 << " to " << outGPFfileName << "." << std::endl;

          // write PRNG engine & distribution state (in case run is to be continued)
          writePRNGState(PRNG_ENG);
          pop.writeDistrStates();

          std::cout << "Written current PRNG & distribution states to prng.state files." << std::endl;

          // update counter for generation output
          if(gen == 0)
          {
            genForOut += genForOutInterval - 1;
          }
          else
          {
            genForOut += genForOutInterval;
          }
        }
        // reproduce, mutate, recombine
        pop.reproduceNoiseGenotypes(PRNG_ENG);
        pop.mutateNoiseGenotypes(PRNG_ENG);
        pop.recombineNoiseGenotypes(PRNG_ENG);
      }
      // close connections to output files
      //outGPFfile.close();
      outSummarizedGPFfile.close();
      // report
      std::cout << std::endl << "Written genotypes/phenotypes/fitness of all generations to " << outGPFfileName << "." << std::endl;
    }
  }

  // report elapsed time
  auto simStop = std::chrono::steady_clock::now();
  auto simElapsedTime =  std::chrono::duration<double>(simStop - simStart);
  std::ofstream outDonefile;
  outDonefile.open("enoise.done", std::fstream::out);
  outDonefile << "Simulation completed in " << simElapsedTime.count()/60 << " minutes." << std::endl;
  outDonefile.close();
}
