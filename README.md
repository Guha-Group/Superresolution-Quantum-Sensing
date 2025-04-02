# Superresolution Quantum Sensing
__Author__: Nico Deshler

## Introduction
This repository contains a simulation pipeline for modelling a multi-stage receiver capable of performing sub-diffraction sensing of vacancy color centers. The receivers employ frequentist and bayesian estimation tools to return an estimate of the immediate brightness of an ensemble of color centers whose separations are below the diffraction limit of an imaging microscope.

## Components
For pedagogical purposes, we break up the repository into two estimation contexts. The first context employs the multi-stage sensing protocol for a system consisting of 2 color centers separated by a sub-diffraction length along one dimension. Scripts associated with this system can be found in the 
```2-Source/``` directory. The second context employs our multi-stage sensing protocol for a system consisting of $K$ color centers, all of which reside in a sub-diffraction spot in two dimensions. Scripts associated with this system can be found in the ```K-Source``` directory.
