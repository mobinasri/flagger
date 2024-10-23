# Testing hmm_flagger with simulated coverage data

## Overview
It is possible to simulate a coverage file with previously defined emission and transition probabilities. For this simulation
the python script, [simulate_coverage_data.py](https://github.com/mobinasri/flagger/blob/dev-hmm-flagger-v1.0.0/programs/src/simulate_coverage_data.py) can be used.

```
python3 /home/programs/src/simulate_coverage_data.py -h
usage: simulate_coverage_data.py [-h] [--pathToEmission PATHTOEMISSION] [--pathToTransition PATHTOTRANSITION] [--pathOutput PATHOUTPUT]
                                 [--numberOfObservations NUMBEROFOBSERVATIONS] [--regionChangeRate REGIONCHANGERATE] [--contigLengths CONTIGLENGTHS]
                                 [--alpha ALPHA]

Simulate coverage data for running hmm_flagger.

optional arguments:
  -h, --help            show this help message and exit
  --pathToEmission PATHTOEMISSION
                        Path to the tsv file that contains emission parameters for different states and regions.
  --pathToTransition PATHTOTRANSITION
                        Path to the tsv file that contains transition matrices for different regions.
  --pathOutput PATHOUTPUT
                        Path for the output tsv file that contains the simulated observations, states and regions. (Default= "observations.cov")
  --numberOfObservations NUMBEROFOBSERVATIONS
                        Total number of observations for simulation.(Default = 10000)
  --regionChangeRate REGIONCHANGERATE
                        Rate of changing regions (will be ignored if there is only one region). (Default= 0.001)
  --contigLengths CONTIGLENGTHS
                        A comma separated list of numbers. The sum of numbers should be equal to the number of observations. For example for
                        --numberOfObservations 100 users can pass --contigLengths 30,40,30 (Default= one contig covering all observations)
  --alpha ALPHA         The dependency factors of the current emission density to the previous emission (Only works for Gaussian). It should be a
                        comma-separated string of 5 numbers for these states respectively err,dup,hap,col,trans. (trans is for transitioning from one
                        state to a different one) [Default = all alpha factors are set to 0]
```

## Input files for simulate_coverage_data.py

`simulate_coverage_data.py` needs two main tsv files contating the HMM parameters:

### Emission TSV
- A TSV file for emission parameters: This tsv file should contain the distribution types and parameters for each state of the HMM that we want
  to use for simulating data. If there are multiple regions with different HMM parameters the parameters for each region should be written
  in a separate column. Below is an example of such a file:
```
#State	Distribution	Components	Parameter	Values_Region_0	Values_Region_1
Err	Gaussian	1	Mean	2.0	3.0
Err	Gaussian	1	Var	2.4	4.0
Err	Gaussian	1	Weight	1.0	1.0
Dup	Gaussian	1	Mean	10.0	15.0
Dup	Gaussian	1	Var	12.0	20.0
Dup	Gaussian	1	Weight	1.0	1.0
Hap	Gaussian	1	Mean	20.0	30.0
Hap	Gaussian	1	Var	24.0	40.0
Hap	Gaussian	1	Weight	1.0	1.0
Col	Gaussian	4	Mean	40.0,60.0,80.0,100.0	60.0,90.0,120.0,150.0
Col	Gaussian	4	Var	48.0,72.0,96.0,120.0	80.0,120.0,160.0,200.0
Col	Gaussian	4	Weight	0.4,0.3,0.2,0.1	0.5,0.4,0.05,0.05
```
The columns of emission tsv file:
1. `#State`: The name of HMM state. It can be either "Err", "Dup", "Hap" or "Col"
2. `Distribution`: The name of the distribution for modeling the emission for each state. It can be either "Gaussian", "Truncated Exponential" or "Negative Binomial"
3. `Components`: The number of components for each state.
4. `Parameter`: The name of the parameter that will be defined in `Values_Region_` columns. It can be either "Mean", "Var", "Weight" or "Trunc_Point" . "Trunc_Point" is valid only for "Truncated Exponential".
5. `Values_Region_X`: (X should be replaced with a number). Each entry in this column is a comma-separated list of numbers (One number don't need a comma). These numbers are the values of the specified parameter for all components per state and region. The region index can be identifed by the ending number in the column name.


### Transition TSV
- A TSV file for transition parameters: This tsv file should contain the probabilities of transitioning between HMM states for all regions. Below is an example of such a file:
```
#Region	State	Err	Dup	Hap	Col	End
0	Err	0.8999	0.01	0.08	0.01	1.0e-4
0	Dup	0.01	0.8999	0.08	0.01	1.0e-4
0	Hap	0.001	0.005	0.9899	0.004	1.0e-4
0	Col	0.02	0.02	0.06	0.8999	1.0e-4
0	Start	2.50e-1	2.50e-01	2.50e-01	2.50e-01	0.0
1	Err	0.8999	0.04	0.05	0.01	1.0e-4
1	Dup	0.01	0.8999	0.08	0.01	1.0e-4  
1	Hap	0.003	0.003	0.9899	0.004	1.0e-4
1	Col	0.02	0.03	0.05	0.8999	1.0e-4
1	Start	2.50e-1	2.50e-01	2.50e-01	2.50e-01	0.0
```

The columns of transition tsv file:
1. `#Region`: The region index. The index should start from 0 and the number of region indices should match the number of `Values_Region_` columns in the emission tsv file. For example in this example we have region indices of 0 and 1 so the number of `Values_Region_` columns should be 2.
2. `State`: The name of HMM state from which the transition begins. It can be either "Err", "Dup", "Hap" or "Col".
3. The remaining columns are "Err", "Dup", "Hap" and "End". Each of them is the HMM state that the related transition ends in. "End" is for termination.

Since the start and end transition values are not important for testing hmm_flagger users can set start probabilities to 0.25 and end probabilites to an arbitrary small number (e.g. 1e-4).

If we keep all rows with the same region index we can extract the transition matrix for that region.

### Run simulate_coverage_data.py

Three examples of emission tsv file are available in flagger git directory `programs/tests/test_files/simulate_coverage/` and also one transition tsv file.
- `emission_exp_gaussian.tsv`
- `emission_gaussian.tsv`
- `emission_negative_binomial.tsv`
- `transition.tsv`

Here we use `emission_gaussian.tsv` and `transition.tsv` to generate a coverage file with 100k observations (bases). These observations are put in two contigs with lengths of 80k and 20k (`--contigLengths 80000,20000`). Setting `--contigLengths` is optional.
```
cd programs/tests/test_files/simulate_coverage/

docker run --rm -v$PWD:/data \
    mobinasri/flagger:v1.1.0-alpha \
    python3 /home/programs/src/simulate_coverage_data.py \
    --pathToEmission emission_gaussian.tsv \
    --pathToTransition transition.tsv \
    --pathOutput test_gaussian_100k.cov \
    --numberOfObservations 100000 \
    --regionChangeRate 0.001 \
    --contigLengths 80000,20000
```

`--regionChangeRate 0.001` means that it will randomly change the region while generating the observations. With a probability of 0.001 the average length of contiguous blocks with the same region will be 1000.

### Run hmm_flagger on the generated coverage file

```
cd programs/tests/test_files/simulate_coverage/
mkdir -p hmm_flagger_runs/gaussian_100k

docker run --rm -it -u$(id -u):$(id -g) -v$PWD:/data \
    mobinasri/flagger:v1.1.0-alpha \
    hmm_flagger_new \
    --input test_gaussian_100k.cov \
    --outputDir hmm_flagger_runs/gaussian_100k \
    --modelType gaussian \
    --trackName gaussian_100k \
    --chunkLen 1000 \
    --windowLen 1 \
    --labelNames Err,Dup,Hap,Col \
    --threads 8 \
    --convergenceTol 1e-4 \
    --collapsedComps 4 \
    --initialRandomDev 0.25
```

`--windowLen 1` is required when running hmm_flagger with simulated coverage file (It means each base is an observation). `--initialRandomDev` is set to 0.25 to show that initial parameter values may not be exact and EM can work properly with inexact initial values to some extent (for example here 25% off).

List of files generated in hmm_flagger_runs/gaussian_100k:

```
ls hmm_flagger_runs/gaussian_100k/
emission_final.tsv	      prediction_summary_final.benchmarking.auN_ratio.tsv    prediction_summary_initial.benchmarking.tsv
emission_initial.tsv	      prediction_summary_final.benchmarking.tsv		     prediction_summary_initial.tsv
final_flagger_prediction.bed  prediction_summary_final.tsv			     transition_final.tsv
loglikelihood.tsv	      prediction_summary_initial.benchmarking.auN_ratio.tsv  transition_initial.tsv
```

`emission_final.tsv` contains the values of emission parameters fit by EM:
```
cat hmm_flagger_runs/gaussian_100k/emission_final.tsv 
#State	Distribution	Components	Parameter	Values_Region_0	Values_Region_1
Err	Gaussian	1	Mean	2.00e+00	3.01e+00
Err	Gaussian	1	Var	2.40e+00	3.99e+00
Err	Gaussian	1	Weight	1.00e+00	1.00e+00
Dup	Gaussian	1	Mean	1.00e+01	1.50e+01
Dup	Gaussian	1	Var	1.20e+01	1.99e+01
Dup	Gaussian	1	Weight	1.00e+00	1.00e+00
Hap	Gaussian	1	Mean	2.00e+01	3.01e+01
Hap	Gaussian	1	Var	2.40e+01	3.99e+01
Hap	Gaussian	1	Weight	1.00e+00	1.00e+00
Col	Gaussian	4	Mean	4.00e+01,6.00e+01,8.00e+01,1.00e+02	6.02e+01,9.03e+01,1.20e+02,1.50e+02
Col	Gaussian	4	Var	4.80e+01,7.19e+01,9.59e+01,1.20e+02	7.97e+01,1.20e+02,1.59e+02,1.99e+02
Col	Gaussian	4	Weight	3.96e-01,2.95e-01,2.17e-01,9.29e-02	5.12e-01,4.00e-01,4.42e-02,4.35e-02
```


`transition_final.tsv` contains the values of transition parameters fit by EM:
```
cat hmm_flagger_runs/gaussian_100k/transition_final.tsv 
#Region	State	Err	Dup	Hap	Col	End
0	Err	8.88e-01	1.41e-02	8.94e-02	7.98e-03	1.00e-04
0	Dup	9.98e-03	8.87e-01	9.41e-02	8.34e-03	1.00e-04
0	Hap	1.12e-03	5.60e-03	9.89e-01	4.59e-03	1.00e-04
0	Col	2.34e-02	1.57e-02	6.61e-02	8.95e-01	1.00e-04
0	Start	2.50e-01	2.50e-01	2.50e-01	2.50e-01	0.00e+00
1	Err	8.99e-01	4.09e-02	5.08e-02	9.07e-03	1.00e-04
1	Dup	1.34e-02	8.99e-01	7.73e-02	1.03e-02	1.00e-04
1	Hap	3.37e-03	2.87e-03	9.90e-01	4.02e-03	1.00e-04
1	Col	1.79e-02	3.46e-02	4.72e-02	9.00e-01	1.00e-04
1	Start	2.50e-01	2.50e-01	2.50e-01	2.50e-01	0.00e+00
```

One can compare these numbers with the truth tsv files passed to `/home/programs/src/simulate_coverage_data.py`.
  
