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

## 4) Jenson-Shannon

The Jenson-Shannon divergence is a distance function that is studied for the purpose of stochastic model updating via ABC. There within lies two folders:
1) Numerical-Example: Run the MATLAB file titled: "Illustrative_example.m" which provides a numerical set-up to demonstrate the mathematical properties of the Jenson-Shannon divergence distance function.
2) NASA-LaRC-Challenge: Run the MATLAB file titled: "NASA_LaRC_Challenge_Part1.m", followed by "NASA_LaRC_Challenge_Part2.m" to evluate the Bayesian model updating results on a black-box model using ABC and the subsequent model validation performance. The focus of the case study is on problems A1 and A2 of the challenge.

## Reference:
* A. Lye, S. Ferson, and S. Xiao (2024). Comparison between distance functions for Approximate Bayesian Computation towards Stochastic model updating and Model validation under limited data. *ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems Part A: Civil Engineering, 10*, 03124001. doi: [10.1061/AJRUA6.RUENG-1223](https://ascelibrary.org/doi/10.1061/AJRUA6.RUENG-1223)
* A. Lye, S. Ferson, and S. Xiao (2024). Stochastic Model Updating Using The Jenson-Shannon Divergence For Calibration and Validation Under Limited Data. *In the Proceedings of the 34th European Safety And Reliability Conference, 1*. Link to paper: [Click here](https://esrel2024.com/wp-content/uploads/articles/part1/stochastic-model-updating-using-jenson-shannon-divergence-for-calibration-and-validation-under-limited-data.pdf)

## Author:
* Name: Adolphus Lye
* Contact: snrltsa@nus.edu.sg
* Affiliation: Singapore Nuclear Research and Safety Institute, National University of Singapore
