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
