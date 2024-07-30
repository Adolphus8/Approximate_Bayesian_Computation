function g=x_and_d_to_g(x,d) %#codegen
% INPUTS: 
% x: an nx5 matrix of intermediate parameters where n is the 
% number of samples.
% d: a 1x14 vector of design variables. The baseline design is
% dbaseline=[0.0189,-0.2410,-0.1003,-2.1800,0.8191,-0.1642,-0.0981,...
%            -0.4362,-0.5958,0.3230,0.0053,-0.1858,1.5875,-0.0378];
% Note that x_and_d_to_g(x,dbaseline)=x_to_g(x)
%
% OUTPUT:
% g: an nx8 matrix where the jth column corresponds to the jth requirement
% function and the ith row corresponds to the ith sample of x
%

% CODED BY:  Luis G. Crespo NIA/NASA-DSCB
% CODED ON:  31 August, 2012

