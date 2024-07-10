# Testing hmm_flagger with simulated coverage data

## Simulating coverage data
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
                        Path for the output tsv file that contains the simulated observations, states and regions. (Default= "observations.tsv")
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

## How to run simulate_coverage_data.py

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
3. The remaining columns are "Err", "Dup", "Hap" and "End". Each of them is a the HMM state that the related transition ends in. "End" is for termination.

Since the start and end transition values are not important for testing hmm_flagger users can set start probabilities to 0.25 and end probabilites to an arbitrary small number (e.g. 1e-4).

If we keep all rows with the same region index we can extract the transition matrix for that region.

  
