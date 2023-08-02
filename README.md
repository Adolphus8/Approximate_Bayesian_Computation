# Distance Metrics and Approximate Bayesian Computation Tutorial
This is a repository on containing MATLAB codes designed to demonstrate the implementation of Approximate Bayesian Computation (ABC) and the different distance metrics which can be implemented for Bayesian model updating and model validation.

## 1) Numerical study

The numerical study serves to demonstrate how the following distance metrics, employed in ABC, quantify the stiatistical distance (or difference) between any two sample sets and evaluate the strengths and weaknesses of each distance metric from there:
1) Euclidean distance metric (EDMe.m)
2) Bhattacharyya distance metric (BDMe.m)
3) Bray-Curtis distance metric (BCMe.m)
4) 1-Wasserstein distance metric (areaMe.m)

## 2) The 2014 NASA-LaRC Uncertainty Quantification Challenge

The [2014 NASA-LaRC Uncertainty Quantification Challenge](https://uqtools.larc.nasa.gov/nda-uq-challenge-problem-2014/) serves as the first case study to evluate the compare the Bayesian model updating results on a black-box model using ABC and the subsequent model validation performance between the different distance metrics discussed. 
The focus of the case study is on problems A1 and A2 of the challenge.

Note that the user is required to run the following codes in sequence:
1) NASA_LaRC_Challenge_Part1.m
2) NASA_LaRC_Challenge_Part2.m
3) NASA_LaRC_Challenge_Part3.m

## 3) The 2008 Thermal Problem

The [2008 Thermal Problem](https://www.sciencedirect.com/science/article/pii/S004578250700504X) serves as the second case study to evluate the compare the Bayesian model updating results on a white-box, physics-based model using ABC and the subsequent model validation performance between the different distance metrics discussed. 
The focus of the case study is on the ensemble validation procedure.

Note that the user is required to download the [RiskCalc MATLAB package](https://github.com/Institute-for-Risk-and-Uncertainty/pba-for-matlab.git) as a folder and add it to path before running the following codes in sequence:
1) Thermal_Problem_Part1.m
2) Thermal_Problem_Part2.m

## Author:
* Name: Adolphus Lye
* Contact: snrltsa@nus.edu.sg
* Affiliation: Singapore Nuclear Research and Safety Initiatives, National University of Singapore
