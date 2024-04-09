# Testing EM algorithm implemented in hmm submodule

## Using HMM API in hmm_test
[hmm_test](https://github.com/mobinasri/flagger/blob/dev-hmm-flagger-v1.0.0/programs/src/hmm_test.c)  is a program that can be used for testing
the HMM API implemented in [hmm.h](https://github.com/mobinasri/flagger/blob/dev-hmm-flagger-v1.0.0/programs/submodules/hmm/hmm.h) and fit HMM parameters to a simulated coverage data.
It can be executed via this experimental docker image, `mobinasri/hmm_flagger_dev:v1.0.0--415ad2b66c6f2114bbda88a8d5f8f61341629f26`.
Here is its help message:
```
Usage: hmm_test
Options:
         --testData, -i		        path to the test file
         --outputDir, -o		directory for saving output files.
         --iterations, -n		maximum number of iterations [Default = 100]
         --convergenceTol, -t		convergence tolerance [Default = 0.001]
         --model, -m			model type can be either 'gaussian', 'negative_binomial', or 'trunc_exp_gaussian' [Default = not defined]
         --coverage, -c		        median coverage for initializing the parameters for the EM algorithm. Note that the inital values will be slightly deviated from the given value by a random factor to make sure the EM algorithm can handle imprecise initial values.
         --collapsedComps, -p		number of components of the collapsed state [Default = 4]
         --writeStatsPerIteration, -w	write emission, transition and posterior statisitics per each iteration.
         --numberOfRegions, -r		number of regions [Default = 1]
         --regionScales, -s		a comma-delimited list of scaling factors for setting initial means of different regions [Default = '1.0'] (As an example with --numberOfRegions 3 --regionScales can be set to '1.0,0.5,1.25')
```

The main input to `hmm_test` is a test tsv file (`--testData, -i`) that contains a simulated coverage data. For generating this input data 
the python script, [simulate_coverage_data.py](https://github.com/mobinasri/flagger/blob/dev-hmm-flagger-v1.0.0/programs/src/simulate_coverage_data.py) can be used.

```
python3 /home/programs/src/simulate_coverage_data.py -h
usage: simulate_coverage_data.py [-h] [--pathToEmission PATHTOEMISSION] [--pathToTransition PATHTOTRANSITION] [--pathOutput PATHOUTPUT]
                                 [--numberOfObservations NUMBEROFOBSERVATIONS] [--regionChangeRate REGIONCHANGERATE]

Simulate coverage data for running hmm_test.

optional arguments:
  -h, --help            show this help message and exit
  --pathToEmission PATHTOEMISSION
                        Path to the tsv file that contains emission parameters for different states and regions.
  --pathToTransition PATHTOTRANSITION
                        Path to the tsv file that contains transition matrices for different regions.
  --pathOutput PATHOUTPUT
                        Path for the output tsv file that contains the simulated observations, states and regions. (Default= "observations.tsv")
  --numberOfObservations NUMBEROFOBSERVATIONS
                        Total number of observations for simulation.(Default = 10000)
  --regionChangeRate REGIONCHANGERATE
                        Rate of changing regions (will be ignored if there is only one region). (Default= 0.001)
```

## How to run simulate_coverage_data.py

### Using Gaussian model

`simulate_coverage_data.py` needs two main tsv file contating the HMM parameters:
- A TSV file for emission parameters: This tsv file should contain the distribution types and parameters for each state of the HMM that we want
  to use for simulating data. If there are multiple regions with different HMM parameters the parameters for each region should be written
  in a separate column.
- A TSV file for transition parameters: This tsv file should contain the distribution types and parameters for each state of the HMM that we want
  to use for simulating data.
  
