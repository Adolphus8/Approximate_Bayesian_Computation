function [f] = makePDF(post_samps, xin)
%% Function-handle to construct PDF via Kernel Density Estimates:
%
% Inputs:
% post_samp: the N x 1 x dm vector of input samples;
% xin:       the 1 x n input vector of x values to compute the PDF;
%
% Output:
% f:         the 1 x n output vector of the PDF estimates;
%% Define the parameters:
dm = size(post_samps,3); % The no. of distance metrics considered
Nbins = 10;

for j = 1:dm
width = (max(post_samps(:,1,j)) - min(post_samps(:,1,j)))./Nbins;
pd = fitdist(post_samps(:,1,j), 'Kernel', 'Width', width);
f(j,:) = pdf(pd, xin)./max(pdf(pd, xin));
end

end