function [logL] = loglikelihood(theta, data, model, dm, width, Nsim)
%% The Log-likelihood function defined for the Bayesian model updating procedure:
%
% Inputs:
% theta: the N x dim matrix of input samples;
%
% theta_1 = E[p1]; - Mean of p1
% theta_2 = V[p1]; - Variance of p1
% theta_3 = delta; - Epistemic interval of p2
% theta_4 = E[p4]; - Mean of p4
% theta_5 = V[p4]; - Variance of p4
% theta_6 = E[p5]; - Mean of p5
% theta_7 = V[p5]; - Variance of p5
% theta_8 = rho;   - Correlation coefficient between p4 and p5
%
% data:  the 25 x 1 vector of training data of x1;
% model: the black-box function-handle of h1;
% dm:    the scalar index to decide the choice of distance metric;
%        - 1 for Euclidean distance metric
%        - 2 for Bhattacharyya distance metric
%        - 3 for Bray-Curtis distance metric
%        - 4 for 1-Wasserstein distance metric
% width: the scalar value of the width parameter of the loglikelihood function;
% Nsim:  the scalar number of simulated realisations of the stochastic model output;
%
% Output:
% logL:  the N x 1 vector of loglikelihood values;
%% Define the parameters:
N = size(theta,1); % The input sample size

%% Define the models:

% Beta shape parameter models:
alpha = @(x) ((x(:,1)./x(:,2)).*(1 - x(:,1)) - 1) .* x(:,1);
beta  = @(x) ((x(:,1)./x(:,2)).*(1 - x(:,1)) - 1) .* (1 - x(:,1));

% Multivariate Normal shape parameter models:
mu =  @(x) [x(:,1), x(:,2)]; 
cov = @(x) [x(:,1), (x(:,3).*sqrt(x(:,1)).*sqrt(x(:,2))); (x(:,3).*sqrt(x(:,1)).*sqrt(x(:,2))), x(:,2)];

% Define the black-box model inputs:
p1 = @(x) betarnd(alpha(x), beta(x), Nsim, 1);
p2 = @(x) x .* ones(Nsim,1);

% Define the distance metric:
if dm == 1
d = @(x,y) EDMe(x,y);
elseif dm == 2
d = @(x,y) BDMe(x,y,[]);    
elseif dm == 3
d = @(x,y) BCMe(x,y);  
elseif dm == 4
d = @(x,y) areaMe(x,y);  
elseif dm == 5
d = @(x,y) JSdiv(x,y);  
end

%% Computation procedure:

logL = zeros(N,1); 
for i = 1:N
input = [p1(theta(i,[1:2])), p2(theta(i,3)), unifrnd(0, 1, Nsim, 1), mvnrnd(mu(theta(i, [4,6])), cov(theta(i, [5,7,8])), Nsim)];
x1_sim = model(input); dist = d(x1_sim,data);

logL(i) = - (dist./width).^2;

if isinf(logL(i)) || isnan(logL(i))
logL(i) = -1e+100;    
end
end

end

function ed = EDMe(sample_1,sample_2)
% Return the Euclidean distance between y_sim and y_exp
%%
% INPUT: 
% sample_1: N x dim matrix of simulated model output y_sim samples;
% sample_2: N x dim matrix of data y_exp samples;
%
% OUTPUT: 
% ed:       The scalar value of the Euclidean distance;
%%
% Define the variables:
[Nsamp1, dim1] = size(sample_1); [Nsamp2, dim2] = size(sample_2);

if dim1 ~= dim2
error('No. of column(s) of the two samples must be equal to each other.')    
end

mean_diff_squared = zeros(dim1,1);
for i = 1:dim1
mean_diff_squared(i) = (mean(sample_1(:,i)) - mean(sample_2(:,i))).^2;    
end

ed = sqrt(sum(mean_diff_squared));
end

function bd = BDMe(sample_1, sample_2, Nbin)
% Return the Bhattacharrya distance between y_sim and y_exp based on the binning algorithm
%%
% INPUT: 
% sample_1: N x dim matrix of simulated model output y_sim samples;
% sample_2: N x dim matrix of data y_exp samples;
% Nbin:     Scalar number of bins
%
% OUTPUT: 
% bd:       The scalar value of the Bhattacharrya distance;
%%
% Define the variables:
[Nsamp1, dim1] = size(sample_1); [Nsamp2, dim2] = size(sample_2);

if dim1 ~= dim2
error('No. of column(s) of the two samples must be equal to each other.')    
end

%% 
% Compute number of bins using the Empirical formula:

if isempty(Nbin) == 1
    
nbin = ceil(max(Nsamp1,Nsamp2)./10);
if nbin <= 10
nbin = 11;
end

Nbin = nbin;

end

%%
% Compute the PMF function for each sample:

max_1 = max(sample_1); min_1 = min(sample_1);
max_2 = max(sample_2); min_2 = min(sample_2);

ub = max([max_1; max_2]); lb = min([min_1; min_2]);

% This treatment is necessary to avoid unexpected error in the following frequenty counting process:
ub = ub + abs(ub*0.0001); lb = lb - abs(lb*0.0001); 

intervals = zeros(Nbin + 1, dim1);
for icolumn = 1:dim1
intervals(:, icolumn) = linspace(lb(icolumn), ub(icolumn), Nbin + 1);
end

% Compute the histogram counts:
count_1 = histcounts(sample_1, intervals); count_2 = histcounts(sample_2, intervals);
count_ratio_1 = count_1/Nsamp1; count_ratio_2 = count_2/Nsamp2;

% Compute the distance metric:
dis = zeros(Nbin, 1);
m = 1;
for i = 1:Nbin
dis(m) = sqrt(count_ratio_1(i)*count_ratio_2(i));
m = m + 1;
end

bd = -log(sum(dis));
end

function bc = BCMe(sample_1, sample_2)
% Return the Bray-Curtis distance between y_sim and y_exp based on the adaptive binning algorithm
%%
% INPUT: 
% sample_1: N x dim matrix of simulated model output y_sim samples;
% sample_2: N x dim matrix of data y_exp samples;
%
% OUTPUT: 
% bc:       The scalar value of the Bray-Curtis distance;
%%
% Define the variables:
[Nsamp1, dim1] = size(sample_1); [Nsamp2, dim2] = size(sample_2);

if dim1 ~= dim2
error('No. of column(s) of the two samples must be equal to each other.')    
end

%% 
% Initiate the Adaptive-binning algorithm:

EDme = EDMe(sample_1,sample_2); % Compute the Euclidean distance

delta_vec = zeros(dim1,1);
for i = 1:dim1
delta_mat = zeros(Nsamp1, Nsamp1);    
for j = 1:Nsamp1
for k = 1:Nsamp1
delta_mat(j,k) = abs(sample_1(j,i) - sample_1(k,i));
end
end
delta_vec(i,1) = max(delta_mat, [], 'all');
end

delta_sim = max(delta_vec);
bin_width = (log(delta_sim + 1)./max([Nsamp1^(1/3), Nsamp2^(1/3)])) .* exp(EDme); % Compute the bin width
Nbin = ceil(delta_sim./bin_width); % Compute number of bins

if Nbin < 2
Nbin = 2;
elseif Nbin > ceil(max(Nsamp1, Nsamp2)./10)
Nbin = ceil(max(Nsamp1, Nsamp2)./10);
else
Nbin = Nbin;
end

%%
% Compute the PMF function for each sample:

max_1 = max(sample_1); min_1 = min(sample_1);
max_2 = max(sample_2); min_2 = min(sample_2);

ub = max([max_1; max_2]); lb = min([min_1; min_2]);

% This treatment is necessary to avoid unexpected error in the following frequenty counting process:
ub = ub + abs(ub*0.0001); lb = lb - abs(lb*0.0001); 

intervals = zeros(Nbin + 1, dim1);
for icolumn = 1:dim1
intervals(:, icolumn) = linspace(lb(icolumn), ub(icolumn), Nbin + 1);
end

% Compute the histogram counts:
count_1 = histcounts(sample_1, intervals); count_2 = histcounts(sample_2, intervals);
count_ratio_1 = count_1/Nsamp1; count_ratio_2 = count_2/Nsamp2;

% Compute the distance metric:
dis_num = zeros(Nbin, 1); dis_den = zeros(Nbin, 1);
m = 1;
for i = 1:Nbin
dis_num(m) = abs(count_ratio_1(i) - count_ratio_2(i)); % Compute numerator term
dis_den(m) = count_ratio_1(i) + count_ratio_2(i);      % Compute denominator term
m = m + 1;
end

bc = sum(dis_num)./sum(dis_den);
end

function AM = areaMe(d1,d2,varargin)
%AREAME Computes the stochastic area metric between two data sets.
%   This version works with data sets of different sizes.
% 
%    %-------------------------------%
%    Author: Marco De Angelis  
%    _______ Created Oct 2020
%    _______ github.com/marcodeangelis
%    _______ University of Liverpool
%    %-------------------------------%
%
%   This version implements the presence of bounds. This can be useful for
%   example, when the area metric is needed only on a portion of the 
%   domain spanned by the data.
%
%   For example if your data span between -5 and 5, and you wanna know what
%   is the area metric for positive values only, you would do the following: 
%                                           
%   >>> AM = areaMe(d1,d2,[0,Inf]) % bounded [0,Inf]
% 
%   The following statement is equivalent to the unbounded area metric:
% 
%   >>> AM = areaMe(d1,d2,[-inf,Inf]) % unbounded [-Inf,Inf]
% 
%   For the unbounded case never use the above command, for efficiency use
%   this instead:
% 
%   >>> AM = areaMe(d1,d2).
% 
%%

n1 = length(d1);
n2 = length(d2);

bounded = false;
if nargin==3
bounded = true; % bounds have been provided
bounds = varargin{1}; % bounds on the area metric. This sends all of the cdf values to the left of the left-bound to zero; and all of the cdf values to the right of the right-bound to one.
end

[x1,y1]=ecdf_Lpool(d1); % Subroutine
[x1,y1]=stairs(x1,y1); % In-built matlab function
x1 = [x1(1);x1]; % this will make sure the data steps from 0 -> 1 on first datum
y1 = [0;y1]; % this will make sure the data steps from 0 -> 1 on first datum

[x2,y2]=ecdf_Lpool(d2); % Subroutine
[x2,y2]=stairs(x2,y2); % In-built matlab function
x2 = [x2(1);x2]; % this will make sure the data steps from 0 -> 1 on first datum
y2 = [0;y2]; % this will make sure the data steps from 0 -> 1 on first datum

assert(2*n1==length(y1),'The stairs function isnt being used correctly.')

d12 = sort([d1;d2]);

if bounded
    inside_bounds = d12>bounds(1) & d12<bounds(2);
    d12_bounded = d12(inside_bounds);
    if bounds(1)>-Inf
        d12_bounded=[bounds(1);d12_bounded];
    end
    if bounds(2)<Inf
        d12_bounded=[d12_bounded;bounds(2)];
    end
    d12 = d12_bounded;
end

if n1==n2 && ~bounded
    AM = sum(abs(x2-x1))/(2*n1);
elseif n1>n2 || bounded
    y1q = cdfinterpolator(x1,y1,d12); % Subroutine
    y2q = cdfinterpolator(x2,y2,d12); % Subroutine
    xdif = diff(d12);
    ydif = abs(y2q-y1q);
    ydif = ydif(2:end); 
    AM = sum(xdif.*ydif);
end

end

function yq = cdfinterpolator(x,y,xq) % This can be pretty slow on large datasets
n=length(xq);
yq = zeros(n,1);
for i=1:n
    pos = xq(i)>x;
    if any(pos)
        yqi=y(pos);
        yq(i)=yqi(end);
    elseif all(pos)
        yq(i)=1;
    end
end
end

function [xs,ps] = ecdf_Lpool(x)
%ECDF_LPOOL Summary of this function goes here
%%%
%   A quick script for finding the ecdf of data. 
%
%                   
%%%
n = length(x);
p = 1/n;
xs = sort(x(:));
ps = linspace(p,1,n)';
end

function [output] = JSdiv(samp1, samp2)
% Return the Jenson-Shannon divergence between y_sim and y_exp based on the adaptive binning algorithm
%%
% INPUT: 
% samp1 (model_samp): Nsim x dim matrix of simulated model output y_sim samples;
% samp2 (data_samp):  Nobs x dim matrix of y_obs samples;
%
% OUTPUT: 
% output:     The scalar value of the Jenson-Shannon divergence;
%%
% Define the variables:
[Nsamp1, dim1] = size(samp1); [Nsamp2, dim2] = size(samp2);

if dim1 ~= dim2
error('No. of column(s) of the two samples must be equal to each other.')    
end

%% 
% Initiate the Adaptive-binning algorithm:

EDme = EDMe(samp1,samp2); % Compute the Euclidean distance

delta_vec = zeros(dim1,1);
for i = 1:dim1
delta_mat = zeros(Nsamp1, Nsamp1);    
for j = 1:Nsamp1
for k = 1:Nsamp1
delta_mat(j,k) = abs(samp1(j,i) - samp1(k,i));
end
end
delta_vec(i,1) = max(delta_mat, [], 'all');
end

delta_sim = max(delta_vec);
bin_width = (log(delta_sim + 1)./max([Nsamp1^(1/3), Nsamp2^(1/3)])) .* exp(EDme); % Compute the bin width
Nbin = ceil(delta_sim./bin_width); % Compute number of bins

if Nbin < 2
Nbin = 2;
elseif Nbin > ceil(max(Nsamp1, Nsamp2)./10)
Nbin = ceil(max(Nsamp1, Nsamp2)./10);
else
Nbin = Nbin;
end

%%
% Compute the PMF function for each sample:

max_1 = max(samp1); min_1 = min(samp1);
max_2 = max(samp2); min_2 = min(samp2);

ub = max([max_1; max_2]); lb = min([min_1; min_2]);

% This treatment is necessary to avoid unexpected error in the following frequenty counting process:
ub = ub + abs(ub*0.0001); lb = lb - abs(lb*0.0001); 

intervals = zeros(Nbin + 1, dim1);
for icolumn = 1:dim1
intervals(:, icolumn) = linspace(lb(icolumn), ub(icolumn), Nbin + 1);
end

% Compute the histogram counts:
if dim1 == 1     % 1-D case
count_1 = histcounts(samp1, intervals); count_2 = histcounts(samp2, intervals); 
count_ratio_1 = count_1/Nsamp1; count_ratio_2 = count_2/Nsamp2;

% Compute the KL-divergence metric:
dis_unit = zeros(Nbin, 1); 

for i = 1:Nbin
M = 0.5 .* (count_ratio_1(i) + count_ratio_2(i));
dum1 = 0.5.*((count_ratio_1(i) .* log((count_ratio_1(i))./M)));
dum2 = 0.5.*((count_ratio_2(i) .* log((count_ratio_2(i))./M)));

if isnan(dum1)
dum1 = 0;
end

if isnan(dum2)
dum2 = 0;
end

dis_unit(i) = dum1 + dum2;

end

output = sum(dis_unit);

elseif dim1 == 2 % 2-D case
count_1 = histcounts2(samp1(:,1), samp1(:,2), intervals(:,1), intervals(:,2)); 
count_2 = histcounts2(samp2(:,1), samp2(:,2), intervals(:,1), intervals(:,2)); 
count_ratio_1 = count_1/Nsamp1; count_ratio_2 = count_2/Nsamp2;

% Compute the KL-divergence metric:
dis_unit = zeros(Nbin, Nbin); 
for i = 1:Nbin
for j = 1:Nbin
M = 0.5 .* (count_ratio_1(i,j) + count_ratio_2(i,j));
dum1 = 0.5.*((count_ratio_1(i,j) .* log((count_ratio_1(i,j))./M)));  
dum2 = 0.5.*((count_ratio_2(i,j) .* log((count_ratio_2(i,j))./M)));

if isnan(dum1)
dum1 = 0;
end

if isnan(dum2)
dum2 = 0;
end

dis_unit(i,j) = dum1 + dum2;
end
end

output = sum(dis_unit(:));

end
end


