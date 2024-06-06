# netlings
This program is a simulator of the evolution of expression level mean and variability in gene regulatory networks. 
It is the program used in the following study: 
https://doi.org/10.1101/2023.12.04.569843
For more information about the implementation of the gene regulatory network model and the evolutionary simulation algorithm, please read the corresponding sections in the article.

The full code to reproduce the results of this study can be found at: 
https://gitlab.gwdg.de/molsysevol/supplementarydata_adaptivenoise

## How to install
```bash
# Clone this repository
$ git clone https://github.com/NatashaPuzovic/netlings.git

# Go into the repository, generate makefiles and compile
$ cd netlings
$ cmake .
$ make
``` 

## How to run
```bash
# Print info
$ ./netlings

USAGE
Usage for evolutionary simulations:
$ netlings paramScalarFile paramVectorFile populationFile outFile
 Input (3):
   paramScalarFile    parameters (scalar) for the simulations, including the PRNG seed.
   paramVectorFile    parameters (vector) for the simulations.
   populationFile     all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)
                      of the first generation. The genotypes specified in the firstGenPopFile may be realized into the phenotype or not,
                      it does not matter because they will be realized again in any case.
 Output (5):
   outFile.hist       summarized phenotypes (mean phenotype, variance of phenotypes and mean fitness) per generation.
   outFile.noises     genotypes 1: mean noise genotypes per gene per generation.
   outFile.sbasals    genotypes 2: mean sbasal genotypes per gene per generation.
   outFile.nets       genotypes 3: regulatory networks of all individuals in a generation.
   outFile.pop        all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)
                      of the last generation.
   prng.*             PRNG engine & distribution states, in case the run is interrupted.

Usage for realization of genotypes:
$ netlings paramScalarFile paramVectorFile populationFile -realize
 Input (3):
   paramScalarFile            parameters (scalar) for the simulations, including the PRNG seed.
   paramVectorFile            parameters (vector) for the simulations.
   populationFile             all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)
                              of one generation. The genotypes specified in the populationFile may be realized into the phenotype or not,
                              it does not matter because they will be realized again in any case.
 Output (1):
   populationFile.realized    all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)
                              of one generation, with the phenotypes realized.

OPTIONS
   --help       print help
``` 

## Optional features:
- Optimal expression levels are given by the parameter S_OPT in the parameter file.
- Three selection scenarios available: stabilizing, directional and fluctuating selection acting on gene expression level.
- Stochastic or deterministic realization is specified by the parameter STARTING_NOISE_INT_ALL in the parameter file. 
If STARTING_NOISE_INT_ALL = 0, the realization is deterministic. Conversely, the realization is stochastic.


## Non-optional features:
- Fixed network size (NUM_NODES)
- Output: a summarized table with NUM_OUT_GENS rows and (3 x NUM_NODES + 2) columns,
containing summarized genotypes, phenotypes, phenotypic variances and fitness of the population per generation.

Output:

- History file (.hist)
Dimensions: #rows: NUM_OUT_GENS; #col: 2xNUM_NODES + 2; For 40 genes and 500 generations outputted: 500x122
PVF table design:
| # generation | NUM_NODES x mean population-wide expr. level (s) | NUM_NODES x variance of population-wide expr. level (s) | mean pop fitness |
| ------------ | --- #NUM_NODES ---  | --- #NUM_NODES --- | --- #NUM_NODES --- | ---------- |
| ------------ | --- #NUM_NODES ---  | --- #NUM_NODES --- | --- #NUM_NODES --- | ---------- |

Example PVF table:
| # generation | mean P1 in pop  | mean P2 in pop  | ... | var P1 in pop  | var P2 in pop  | ... | mean fitness |
| ------------ | --------------  | --------------  | ... | -------------- | -------------- | ... | ------------ |
| ---- 1 ----- | ---- 60.5 ----  | ---- 32.6 ----  | ... | ---- 520.2 --- | ---- 32.6 ---- | ... | -- 0.001 --- |
| ---- 50 ---- | ---- 59.5 ----  | ---- 33.4 ----  | ... | ---- 59.5 ---- | ---- 33.4 ---- | ... | -- 0.003 --- |
| ------------ | --------------  | --------------  | ... | -------------- | -------------- | ... | ------------ |
| --- 10000 -- | ---- 60.1 ----  | ---- 34.0 ----  | ... | ---- 60.1 ---- | ---- 34.0 ---- | ... | ---- 0.2 --- |

- Noise genotypes file (.noises)
This file contains the average noise genotype per generation in the output.
Noise genotypes table design:
| # generation | NUM_NODES x mean per-gene population-wide noise genotype (eta) |
| ------------ | ------------------------- #NUM_NODES ------------------------- |
| ------------ | ------------------------- #NUM_NODES ------------------------- |

- SBasal genotypes file (.sbasals)
This file contains the average sbasal genotype per generation in the output.
SBasals genotypes table design:
| # generation | NUM_NODES x mean per-gene population-wide sbasal genotype (s^b) |
| ------------ | ------------------------- #NUM_NODES -------------------------- |
| ------------ | ------------------------- #NUM_NODES -------------------------- |

- Net genotypes file (.nets)
This file contains the unique net genotype per generation in the output. Currently, this does not output
the networks since the network in this version is immutable, and therefore there is only one net genotype throughout the sim.
Net genotypes table design:
| # generation | linearized net genotype |
| ------------ | #NUM_NODES x #NUM_NODES |
| ------------ | #NUM_NODES x #NUM_NODES |

- Population file of pop in the last generation (.pop)
This file contains the information about all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness) 
of the last generation. It has the same format as the first generation .pop file used to start the sim.
It's purpose is to output the full information about the last generation, as opposed to the summarized info that 
is continously output in .noises, .sbasals, .nets files as the sim runs.
Population table design:
| # last generation | # individual in pop | noise genotype of individual | sbasal genotype of individual | linearized net genotype of individual | phenotype of individual | fitness of individual |
| ----------------- | -- # individual --- | -------- #NUM_NODES -------- | --------- #NUM_NODES -------- | ------ #NUM_NODES x #NUM_NODES ------ | ------ #NUM_NODES ----- | ------ realnum ------ |
| ----------------- | -- # individual --- | -------- #NUM_NODES -------- | --------- #NUM_NODES -------- | ------ #NUM_NODES x #NUM_NODES ------ | ------ #NUM_NODES ----- | ------ realnum ------ |

- PRNG and distribution states at last generation outputted (prng.*)

## License
This content is published under a CC-BY-ND Creative Commons License.
