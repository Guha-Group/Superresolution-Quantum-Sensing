# Superresolution Quantum Sensing
__Author__: Nico Deshler

## Introduction
This repository contains a simulation pipeline for modelling a multi-stage receiver capable of performing sub-diffraction sensing of vacancy color centers. The receivers employ frequentist and bayesian estimation tools to return an estimate of the immediate brightness of an ensemble of color centers whose separations are below the diffraction limit of an imaging microscope.

## Components
For pedagogical purposes, we break up the repository into two estimation contexts. The first context employs the multi-stage sensing protocol for a system consisting of 2 color centers separated by a sub-diffraction length along one dimension. Scripts associated with this system can be found in the 
```2-Source/``` directory. The second context employs our multi-stage sensing protocol for a system consisting of $K$ color centers, all of which reside in a sub-diffraction spot in two dimensions. Scripts associated with this system can be found in the ```K-Source``` directory.

## Required Packages and Installation Instructions
The K-Source sensing algorithm uses [manopt](https://www.manopt.org/) - a manifold optimization package. This package is used to numerically find the YKL measurements for optimal brightness sensing. It is also used for solving the maximum likelihood estimate of the source brightnesses. Please follow the steps below to install manopt so that it can be used for our codebase.
1. Download manopt
2. Move the zipped package to whichever sub-directory you wish to store matlab packages in. If working on a Windows machine I recommend ```C:\Users\YourUserName\Documents\MATLAB\Toolboxes```. If such a directory does not already exist, make one.
3. Unzip the package
4. Now we're going to add manopt to matlab's ```startup.m``` routine for your system. That way you won't have to worry about shuffling manopt around to particular directories where it is called - it's methods just get used as if they were built-in matlab functions.
5. In the matlab command line enter ```userpath``` and navigate to the resulting path. This is your user directory where MATLAB looks for ```startup.m``` by default.
6. Check if ```startup.m``` already exists. If not, make it.
7. Open ```startup.m``` in the matlab editor and add the following line: ```addpath(genpath('C:\Users\YourUserName\Documents\MATLAB\Toolboxes\manopt'));``` which adds the search path to wherever you unzipped manopt.
8. Save the file and close matlab
9. Reopen matlab and test the path by calling ```manopt_version``` in the matlab command line

## Running the Code

Now that manopt is installed, you can download our github repo either by clicking the green ```<CODE>``` button on the top right of the page or cloning the repository through a git terminal with the command:
```git clone https://github.com/Guha-Group/Superresolution-Quantum-Sensing.git```

To run our receiver protocol on a randomly generated sub-diffraction scene, navigate to, 
```Superresolution-Quantum-Sensting/K-Source/```
and run the script ```TestReceivers.m```. This will run a direct imaging receiver and our proposed receiver on a localization and brightness estimation task. The parameters of scene include:

```num_sources```  : the number of sources

```d_min```        : the minimum pair-wise source separation

```contrast```     : the ratio of the minimum brightness to the maximum brightness

 One can also adjust the number of photons used for the calibration stage $M$ and the sensing stage $N$ as well as the splitting fraction which determines the proportion of photons allocated to direct imaging versus SPADE measurements. Note that throughout the code, we define $\sigma=1$. This value is not to be adjusted - it is merely for mathematical housekeeping. In this regard, all distance-valued parameters used throughout the codebase implicitly are in units of $\sigma$.
