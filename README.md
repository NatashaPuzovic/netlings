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

## Output:

- History file (*.hist)
This file is a record of the phenotypic evolution - the expression level mean and variance per gene through the generations, and the average population-wide fitness in each generation.
Dimensions (rows, columns): (NUM_OUT_GENS, 2xNUM_NODES + 2). For a simulation with 40 genes and 500 generations outputted: 500x82.
History (PVF) table design:

| # generation  | NUM_NODES x mean population-wide gene expr. level (s) | NUM_NODES x variance of population-wide gene expr. level (s)  | mean pop fitness  |
|:---|:---:|:---:|:---:|
| integer num  | real num x #NUM_NODES | real num x #NUM_NODES  |  real num |
| ...  | ... | ... | ... |
| ...  | ... | ... | ... |

Example history file table for 3 genes (threeGenes.hist):

| # generation | mean P1 in pop  | mean P2 in pop  | mean P3 in pop  | var P1 in pop  | var P2 in pop  | var P3 in pop | mean fitness |
|:---|:---:|:---:|:---:|
| 1 | 60.5  | 32.6 | 22.3 | 520.2 | 32.6 | 101.5 | 0.001 |
| 50  | 59.5 | 34.4 | 26.2 | 59.5 | 33.4 | 99.2 | 0.003 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 10000 |  60.1 | 34.0 | ... | 60.1 | 34.0 | 110.5 | 0.2 |

where "#generation" denotes the number of the generation, "mean P1 in pop" denotes mean phenotype (i.e. expression level) of gene 1 in the population in this generation, "var P1 in pop" denotes the variance of the phenotype (i.e. expression level) of gene 1 in the population in this generation, and "mean pop fitness" denotes the mean fitness of the population in this generation.


- Noise genotypes file (*.noises)
This file contains the record of the evolution of noise genotypes - the average noise genotype of each gene per generation in the output.
Dimensions (rows, columns): (NUM_OUT_GENS, NUM_NODES + 1). For a simulation with 40 genes and 500 generations outputted: 500x41.
Noise genotypes table design:

| # generation | NUM_NODES x mean per-gene population-wide noise genotype (eta) |
|:---|:---:|
| integer num | real num x #NUM_NODES |
| ... | ... |
| ... | ... |

Example noise genotypes file table for 3 genes (threeGenes.noises):

| # generation | mean eta1 in pop | mean eta2 in pop | mean eta3 in pop |
|:---|:---:|:---:|:---:|
| 1 | 200 | 200 | 200 |
| 50 | 75.1 | 191.2 | 46.9 |
| ... | ... | ... | ... |
| 10000 | 10.8 | 60.8 | 23.7 |

where "#generation" denotes the number of the generation, "mean eta1 in pop" denotes mean noise genotype of gene 1 (eta1, the inheritable determinant of expression noise) in the population in this generation, calculated as the mean eta1 across all gene1 in the population. It is identically calculated for all genes.


- SBasal genotypes file (*.sbasals)
This file contains the record of the evolution of sbasal (basal/constitutive expression level genotypes - the average sbasal genotype of each gene per generation in the output.
Dimensions (rows, columns): (NUM_OUT_GENS, NUM_NODES + 1). For a simulation with 40 genes and 500 generations outputted: 500x41.
SBasals genotypes table design:

| # generation | NUM_NODES x mean per-gene population-wide sbasal genotype (s^b) |
|:---|:---:|
| integer num | real num x #NUM_NODES |
| ... | ... |
| ... | ... |

Example sbasal genotypes file table for 3 genes (threeGenes.sbasals):

| # generation | mean sb1 in pop | mean sb2 in pop | mean sb3 in pop |
|:---|:---:|:---:|:---:|
| 1 | 13 | 80 | 66 |
| 50 | 27.5 | 78.1 | 68.2 |
| ... | ... | ... | ... |
| 10000 | 47.1 | 83.4 | 70.0 |

where "#generation" denotes the number of the generation, "mean sb1 in pop" denotes mean sbasal genotype of gene 1 (sb1, the inheritable determinant of mean expression level) in the population in this generation, calculated as the mean sb1 across all gene1 in the population. It is identically calculated for all genes.



- Net genotypes file (*.nets)
This file contains the record of the evolution of unique net genotype per generation. Currently, this does not output anything the networks since the network in this version of the program is immutable, and therefore there is only one net genotype throughout the sim.
Net genotypes table design:

| # generation | linearized net genotype |
|:---|:---:|
| integer num | #NUM_NODES x #NUM_NODES |
| ... | ... |
| ... | ... |

- Population file of pop in the last generation (.pop)
This file contains the information about all individuals (#generation|#individual|noise genotype|sbasal genotype|net genotype|phenotype|fitness)  of the last generation. It has the same format as the first generation .pop file used to start the sim.
It's purpose is to output the full information about the last generation, as opposed to the summarized info that is continously output in .noises, .sbasals, .nets files as the sim runs.
Population table design:

| # last generation | # individual in pop | noise genotype of individual | sbasal genotype of individual | linearized net genotype of individual | phenotype of individual | fitness of individual |
|---|---|---|---| 
| integer num | # individual | real num x #NUM_NODES | real num x #NUM_NODES | real num x (#NUM_NODES x #NUM_NODES) | real num x #NUM_NODES | real num |
| ... | ... | ... | ... | ... | ... | ... |
| ... | ... | ... | ... | ... | ... | ... |

Example population file table for 3 genes (threeGenes.pop):

| # last generation | # individual in pop | eta1 | eta2 | eta3 | sb1 | sb2 | sb3 | w11 | w12 | w13 | w21 | w22 | w23 | w31 | w32 | w33 | s1 | s2 | s3 | f |
|---|---|---|---| 
| 10000 | 1 | 55.2 | 0.2 | 10.3 | 17.1 | 93.2 | 56.3 | 0 | 0 | -0.2 | 1.2 | 0 | -0.002 | 0 | 2.7 | 0 | 27.9 | 94.0 | 60.2 | 0.15 |
| 10000 | 2 | 32.3 | 0.2 | 10.3 | 17.1 | 93.2 | 56.3 | 0 | 0 | -0.2 | 1.2 | 0 | -0.002 | 0 | 2.7 | 0 | 19.6 | 93.0 | 49.5 | 0.27 |
| ... | ... | ... | ... | ... | ... | ... |
| 10000 | 1000 | 32.3 | 0.2 | 10.3 | 17.1 | 93.2 | 56.3 | 0 | 0 | -0.2 | 1.2 | 0 | -0.002 | 0 | 2.7 | 0 | 30.3 | 92.8 | 45.2 | 0.27 |


- PRNG and distribution states at last generation outputted (prng.*)

## License
This content is published under a CC-BY-ND Creative Commons License.
