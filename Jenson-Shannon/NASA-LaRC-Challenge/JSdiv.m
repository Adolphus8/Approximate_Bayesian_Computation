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

% Compute the JS-divergence metric:
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

% Compute the JS-divergence metric:
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
dim1 = size(sample_1,2); dim2 = size(sample_2,2);

if dim1 ~= dim2
error('No. of column(s) of the two samples must be equal to each other.')    
end

mean_diff_squared = zeros(dim1,1);
for i = 1:dim1
mean_diff_squared(i) = (mean(sample_1(:,i)) - mean(sample_2(:,i))).^2;    
end

ed = sqrt(sum(mean_diff_squared));
end
