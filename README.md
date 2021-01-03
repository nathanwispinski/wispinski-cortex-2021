# Code and data for Wispinski et al. (2020)

Wispinski, N.J., Stone, S.A., Bertrand, J.K., Zuk, A.A.O., Lavoie, E.B., Gallivan, J.P., & Chapman, C.S. (2020). Reaching for known unknowns: Rapid reach decisions accurately reflect the future state of dynamic probabilistic information. *bioRxiv*. https://doi.org/10.1101/2020.07.31.231563.

This project consists of two experiments which aimed to examine the extent to which reach movement trajectories convey internal predictions about the future state of dynamic probabilistic information.

These data and code (along with videos of the task), are also available at https://osf.io/rt5xv/.

For a first look, you should only need to run one analyze_\*.m script in MATLAB. All extra data and code is provided for deeper dives into the project.
If you would simply like to recreate the figures and data used for the paper, just run RUNALL.m.

In brief, cleanData.m takes the .mat files for each subject and outputs a data file called cleanData.mat.
We've included the cleanData.mat data file already if you want to skip this step.
With cleanData.mat, you can run analyze_accuracy.m, analyze_RT.m, or analyze_reachArea.m to pull figures and results for each dependent measure.

## Prerequisites

This code is written to be compatible with MATLAB 2019b and 2020a on Windows 10 OS.

fminsearchbnd
 - Bound-constrained optimization function using MATLAB's fminsearch
 - Written by John D'Errico
 - Accessible at: https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon

CircStat2012a toolbox
 - Circular statistics toolbox for MATLAB
 - Written by Philipp Berens
 - Accessible at: https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

## Code

cleanData.m
cleanData.mat
 - Loads individual .mat files with each participant's data from the Data folder.
 - Preprocesses and cleans data. Outputs cleanData.mat for further use with analyze_\*.m scripts (cleanData.mat is already provided in case you don't want to run cleanData.m).

sine5Cost.m
sine61Cost.m
 - Used by analyze_\*.m scripts.
 - Takes sine wave model parameters and data. Outputs sum of squared deviations between data and model.

analyze_RT.m
analyze_accuracy.m
analyze_reachArea.m
 - Takes cleanData.mat and performs sine wave model analysis. Generates summary tables for use in statisitcs programs, and figures.
 - Script titles denote the dependent variable analyzed in that script (i.e., reaction time, accuracy, or reach area).
 - The 3 analysis scripts are identical except for small differences between the dependent variables (e.g., plot limits, parameter fitting bounds).
 - Note: Because of randomness in the selection of starting parameter estimates for sine wave fitting, values from the analysis scripts may slightly differ from run-to-run or to the values reported in the paper.

RUNALL.m
 - Runs each of the analyze_\*.m files. You could also just run each of the analyze files individually for the same output.
 - For today's impatient scientist.
 
 ## Data

Individual .mat file contents within Data folder:

block [1 x trials]
 - Trial block for each trial (1 through 8).
rxnTime [1 x trials]
 - Reaction time in seconds (i.e., time from beep to button release).
mvmtTime [1 x trials]
 - Movement time in seconds (i.e., time from button release to screen touch).
targetDrawTime [1 x trials]
 - Time from button release to correct circle filling in (often a single frame).
preConditions [1 x trials]
 - Condition on each trial (of 48 unique conditions).
 - Can find what each of 48 unique conditions means in cleanData.m.
LPreTargets [1 x trials]
 - Position L target is in (1 or 2; where 1=Left, 2=Top).
RPreTargets [1 x trials]
 - Position R target is in (3 or 4; where 3=Right, 4=Bottom).
postTargets [1 x trials]
 - Which circle was the correct target.
trigPos [1 x trials]
 - Go-cue happened when the box/arrow hit one of these positions.
StartPos [1 x trials]
 - Spatial location box/arrow begins to move. Random point on a 120-point circle.
 - When boxDir=2 (i.e., stationary trial), startPos doesn't matter.
beepFlag [1 x trials]
 - 1=beep occured.
beepI [1 x trials]
 - Number of frames before go beep.
beepCtr [1 x trials]
 - Used in stimuli generation.
boxEndCtr [1 x trials]
 - Used in stimuli generation.
finalBoxPos [1 x trials]
 - Position the box/arrow was in at the end of the trial in x and y screen dimension.
endPositions [1 x trials]
 - Not used.
Error [4 x trials]
 - TooEarly (1st column).
 - TooSlow (i.e., MT too long; 2nd column).
 - Miss (3rd column).
 - TimeOut (i.e., RT too long; 4th column).
reach_info [struct]
 - Info used during reach trajectory normalization (process described in [Gallivan & Chapman, 2014](https://www.frontiersin.org/articles/10.3389/fnins.2014.00215/full)).
peakV [1 x trials]
 - Peak velocity (not used).
Ttpv [1 x trials]
 - Time to peak velocity (not used).
fdaMat [struct]
 - Reach trajectories space-normalized to 200 points (see [Gallivan & Chapman, 2014](https://www.frontiersin.org/articles/10.3389/fnins.2014.00215/full)).
 - x, y, and z dimensions are most important here.
