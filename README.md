# trajectory-divergence
MATLAB code, demo, and sample data to analyze data for trajectory divergence as in [Russo et al Neuron 2020].

**Author**: Abigail Russo, May 2019

# Datasets
The data included in this demo were recorded in the Churchland lab at Columbia University. Single-unit neural data were  recorded sequentially from a macaque monkey during the performance of a hand-pedaling task. Monkeys grasped a hand crank and cycled through a virtual environment for a number of prescribed cycles as indicated by a visual cue. The data included here correspond to two seven-cycle conditions (forward bottom-start and backward bottom-start for monkey C).  For more information about data processing and the task, see Russo et al  Neuron 2020.

Data should be formatted as for [jPCA](https://www.dropbox.com/sh/2q3m5fqfscwf95j/AAC3WV90hHdBgz0Np4RAKJpYa?dl=0&preview=NOTES.pdf). Each data structure (e.g. D_m1, located in **M1_sampleData.mat**) contains C elements corresponding to the number of conditions (here, 2). The t x n matrix 'A' contains the trial-averaged firing rates as a function of time for each neuron. All other fields are optional. The 'times' field indicates the time course of each sample in A and the 'analyzeTimes' field indicates which of these times should be analyzed for tangling. For more information on formatting, see the **divergAnalysis** code comments and [jPCA documentation](https://www.dropbox.com/sh/2q3m5fqfscwf95j/AAC3WV90hHdBgz0Np4RAKJpYa?dl=0&preview=NOTES.pdf).

* **SMA_sampleData.mat**: trial-avergaged firing rates from 77 neurons (monkey C)
* **M1_sampleData.mat**: trial-avergaged firing rates from 116 neurons (monkey C)
* **kinematic_sampleData.mat**: trial-avergaged kinematic data corresponding to a variety of behavioral parameters such as hand position and velocity (monkey C)

# Code

* **divergeAnalysis.m**: compute trajectory divergence on data formatted as described above
* **divergence_demo.m**: visualize task kinematics, compute divergence on sample datasets, plot results as in Russo et al Neuron 2020, Figure 7.
* **formattedScatter.m**: makes pretty scatter plots

*Supporting functions*: **AxisMMC.m** (makes pretty axes, used by **formattedScatter**)

