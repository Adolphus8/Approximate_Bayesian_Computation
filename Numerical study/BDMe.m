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