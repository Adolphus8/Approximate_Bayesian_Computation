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