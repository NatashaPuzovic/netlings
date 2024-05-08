# enoise
This program is a simulator of the evolution of expression level variability in gene regulatory networks. 
It is the program used in the following study: 
https://doi.org/10.1371/journal.pcbi.1010982
For more information about the implementation of the gene regulatory network model and the evolutionary simulation algorithm, please read the corresponding sections in the article.

The full code to reproduce the results of this study can be found at: 
https://gitlab.gwdg.de/molsysevol/supplementarydata_expressionnoise

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
$ ./enoise
enoise1.2
Usage: enoise networkFile paramFile outfile [-continue]

 Input:
   networkFile: regulatory matrix that defines the strength of interactions
                of each pair of genes in the network.
   paramFile:   parameters for the simulations, including the PRNG seed.
 Output:
   outfile:     genotype, phenotype and fitness values of individuals in
                specified generations to output.
   prng* files: PRNG engine & distribution states, in case the run is interrupted.

Options:
   -continue    continue an interrupted run
                this option reads the outfile and the PRNG states of a previous run
                and continues the simulation from the last recorded generation.
   --help       print this help
``` 

## Non-optional features:
- All genes are under selection.
- All genes have the same starting level of noise.
- Network size is fixed (NUM_NODES).
- Noise genotypes are mutable.
- Output: a summarized table with NUM_OUT_GENS rows and (2 x NUM_NODES + 2) columns,
containing summarized genotypes, phenotypes and fitness of all individuals.

Dummy output table:
| generation | mean G1 in pop  | mean G2 in pop | ... | var P1 in pop  | var P2 in pop  | ... | mean fitness |
| ---------- | --------------  | -------------- | ... | -------------- | -------------- | ... | ------------ |
| --- 1 ---- | ---- 100 -----  | ---- 100 ----- | ... | ---- 60.5 ---- | ---- 32.6 ---- | ... | -- 0.001 --- |
| --- 50 --- | ---- 90.3 ----  | ---- 30.2 ---- | ... | ---- 59.5 ---- | ---- 33.4 ---- | ... | -- 0.003 --- |
| ---------- | --------------  | -------------- | ... | -------------- | -------------- | ... | ------------ |
| - 10000 -- | ---- 5.02 ----  | ---- 10.5 ---- | ... | ---- 60.1 ---- | ---- 34.0 ---- | ... | ---- 0.2 --- |

## Optional features:
- Optimal expression levels are given by the parameter S_OPT in the parameter file.
- Stochastic or deterministic realization is specified by the parameter STARTING_NOISE_INT_ALL in the parameter file. If STARTING_NOISE_INT_ALL = 0, the realization is deterministic. Conversely, the realization is stochastic.

## License
This content is published under a CC-BY-ND Creative Commons License.
