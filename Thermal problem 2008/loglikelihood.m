function [logL] = loglikelihood(theta, data, bounds, dm, width, Nsim)
%% The Log-likelihood function defined for the Bayesian model updating procedure:
%
% Inputs:
% theta:  the N x dim matrix of input samples;
% data:   the 20 x 1 vector of training data of x1;
% bounds: the 1 x 2 vector of bounds on the Staircase Density Function;
% dm:     the scalar index to decide the choice of distance metric;
%         - 1 for Euclidean distance metric
%         - 2 for Bhattacharyya distance metric
%         - 3 for Bray-Curtis distance metric
%         - 4 for 1-Wasserstein distance metric
% width:  the scalar value of the width parameter of the loglikelihood function;
% Nsim:   the scalar number of simulated realisations of the stochastic model output;
%
% Output:
% logL:  the N x 1 vector of loglikelihood values;
%% Define the parameters:
N = size(theta,1); % The input sample size

%% Define the distance metric:
if dm == 1
d = @(x,y) EDMe(x,y);
elseif dm == 2
d = @(x,y) BDMe(x,y,[]);    
elseif dm == 3
d = @(x,y) BCMe(x,y);  
elseif dm == 4
d = @(x,y) areaMe(x,y);  
end

%% Computation procedure:

logL = zeros(N,1); 
for i = 1:N

x_sim = SDF_rnd(theta, bounds, 2, Nsim); dist = d(x_sim, data);
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

function [samples] = SDF_rnd(theta, bounds, objective, Nsamp)
% This is the function handle of the Staircase Density Function Random Number Generator:

% Inputs:
% theta:     N x dim vector of epistemic Staricase Density Function parameters;
% bounds:    Bounds of the aleatory parameters;
% objective: Numerical objective function flag for the optimization problem;
% Nsamp:     Numerical value of the number of samples to obtain from the joint distribution defined by the Staircase Density Function;

% Output:
% samples:   The Nsamp x dim vector of sample output;

%% Error check:
assert(size(theta,1)==1);

%% Define the function:
objective_func = objective;
N = Nsamp; 
dim = size(theta,2)./4;

theta_a = cell(1);
samples = zeros(N, dim);
for ia = 1:dim
theta_a{ia} = theta(1 + 4*(ia - 1):4*ia);
    
% Fit a staircase density:
[z_i, l, ~] = staircasefit(bounds, theta_a{ia}, objective_func);
l(l < 0) = 0;
    
% Generate samples from the staircase density:
samples(:, ia) = staircasernd(N, z_i, l, bounds);
end

end

function [z_i, l, c_i] = staircasefit(bounds, theta, objective, varargin)
% Calculation of the staircase random variables
%
%     INPUT : bounds    -- prior distribution of aleatory parameters
%             theta     -- samples of the epistemic parameters
%             objective -- objective function flag for the optimization problem
%
%     OUTPUT : z_i -- partitioning points
%              l   -- staircase heights
%              c_i -- centers of the bins
%

n_b = 50;   % n. of bins of staircase RVs
if nargin >= 4
    n_b = varargin{1};
end

z_i = linspace(bounds(1), bounds(2), n_b + 1);   % partitioning points
kappa = diff(bounds)/n_b;                        % subintervals
c_i = z_i(1:end - 1) + kappa/2;                  % centers of the bins

[feasible, ~] = isfeasible(bounds, theta);

theta(3) = theta(3)*theta(2)^(3/2);
theta(4) = theta(4)*theta(2)^2;

if feasible
    Aeq = [kappa*ones(size(c_i));
        kappa*c_i;
        kappa*c_i.^2 + kappa^3/12;
        kappa*c_i.^3 + kappa^3*c_i/4;
        kappa*c_i.^4 + kappa^3*c_i.^2/2 + kappa^5/80];
    beq = [1;
        theta(1);
        theta(1)^2 + theta(2);
        theta(3)+3*theta(1)*theta(2) + theta(1)^3;
        theta(4) + 4*theta(3)*theta(1) + 6*theta(2)*theta(1)^2 + theta(1)^4];
    
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 10000*n_b;
    options.MaxIterations = 1000;
    options.ConstraintTolerance = 1e-6;
    options.StepTolerance = 1e-12;
    switch objective
        case 1
            J = @(l) kappa*log(l)*l';
        case 2
            J = @(l) l*l';
        case 3
            J = @(l) -omega*log(l)';
    end
    
    % Do fmincon starting from a uniform distribution over bounds
    l = fmincon(J, 1/(diff(bounds))*ones(size(c_i)), [], [], Aeq, beq, zeros(size(c_i)), [], [], options);
else
    error('unfeasible set of parameters')
end

l(l < 0) = 0;   % ignore negative values

end

function [Lfeasible, constraints] = isfeasible(bounds, theta)
% Return a column of the prior pdf.
%
%     INPUT : bouds -- prior distribution of aleatory parameters
%             theta -- theta(:, 1): mean
%                      theta(:, 2): variance
%                      theta(:, 3): the third-order central moment
%                      theta(:, 4): the fourth-order central moment
%
%     OUTPUT : Lfeasible   -- the prior pdf
%              constraints -- the moment constraints
%

Nsample = size(theta, 1);
Lfeasible = zeros(Nsample, 1);

theta(:, 3) = theta(:, 3).*theta(:, 2).^(3/2);
theta(:, 4) = theta(:, 4).*theta(:, 2).^2;

for isample = 1:Nsample
    u = bounds(1) + bounds(2) - 2*(theta(isample, 1));
    v = (theta(isample, 1) - bounds(1))*(bounds(2) - theta(isample, 1));
    
    constraints = [bounds(1) - theta(isample, 1);   % g2
        theta(isample, 1) - bounds(2);   % g3
        -theta(isample, 2);   % g4
        theta(isample, 2) - v;   % g5
        theta(isample, 2)^2 - theta(isample, 2)*(theta(isample, 1) - bounds(1))^2 - theta(isample, 3)*...
        (theta(isample, 1) - bounds(1));   % g6
        theta(isample, 3)*(bounds(2) - theta(isample, 1)) - theta(isample, 2)*...
        (bounds(2) - theta(isample, 1))^2 + theta(isample, 2)^2;   % g7
        4*theta(isample, 2)^3 + theta(isample, 3)^2 - theta(isample, 2)^2*diff(bounds)^2;   % g8
        6*sqrt(3)*theta(isample, 3) - diff(bounds)^3;   % g9
        -6*sqrt(3)*theta(isample, 3) - diff(bounds)^3;   % g10
        -theta(isample, 4);   % g11
        12*theta(isample, 4) - diff(bounds)^4;   % g12
        (theta(isample, 4) - v*theta(isample, 2) - u*theta(isample, 3))*...
        (v - theta(isample, 2)) + (theta(isample, 3)- u*theta(isample, 2))^2;   % g13
        theta(isample, 3)^2 + theta(isample, 2)^3 - theta(isample, 4)*theta(isample, 2)];   % g14
    
    Lfeasible(isample) = all(constraints <= 0);
end
end

function x = staircasernd(N, z_i, l, bounds)
% Return a column of parameters sampled from the prior pdf
%
%     INPUT : N       -- n. of samples
%             z_i     -- partitioning points
%             l       -- staircase heights
%             bounds  -- prior distribution of aleatory parameters
%
%     OUTPUT : x -- matrix of samples from x_pdf
%

n_b = length(l);   % n. of bins of staircase RVs
idx = (randsample(length(l), N, true, diff(bounds)/n_b*l));   % select random stair
x = unifrnd(z_i(idx), z_i(idx + 1))';   % generate uniform in each stair
end


