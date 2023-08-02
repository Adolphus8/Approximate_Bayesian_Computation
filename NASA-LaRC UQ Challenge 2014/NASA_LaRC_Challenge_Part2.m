%% NASA LaRC Challenge 2014: Subproblem A
%
% The codes provided here describe the approach to analyse the results from
% the TEMCMC sampling approach via ABC.

clc; clear;
%% Load the training and validation data sets:
load('x1samples1.mat'); load('x1samples2.mat'); load('NASA_LaRC_Challenge_Part1.mat', 'TEMCMC', 'timeTEMCMC')

%% Analysis of samples:
post_samps = zeros(1000,8,4);
for i = 1:4
TEMCMC_struct = TEMCMC{i};
post_samps(:,:,i) = TEMCMC_struct.samples;
end

% Plot the Histogram representation of the posteriors given each distance metric:
figure;
label = {'$E[p_1]$', '$V[p_1]$', '$p_{2}$', '$E[p_4]$', '$V[p_4]$', '$E[p_5]$', '$V[p_5]$', '$\rho$'};
for i = 1:8
subplot(2,4,i)
hold on; box on; grid on;
nbins = 10;
for j = 1:4
histogram(post_samps(:,i,j), nbins, 'Normalization','pdf')
end
set(gca, 'Fontsize', 18); xlabel(label{i}, 'Interpreter', 'latex'); ylabel('Count')
end
legend('Euclidean', 'Bhattacharrya', 'Bray-Curtis', '1-Wasserstein', 'linewidth', 2)


% Plot the approximated Fuzzy-set function via KDE:
Nplot = 500;
xin = [linspace(0.6, 0.8, Nplot); linspace(0.02, 0.04, Nplot); linspace(0, 1, Nplot); linspace(-5, 5, Nplot);... 
       linspace(1/400, 4, Nplot); linspace(-5, 5, Nplot); linspace(1/400, 4, Nplot); linspace(-1, 1, Nplot)];
f = zeros(size(post_samps,3), Nplot, size(post_samps,2));

for i = 1:size(post_samps,2)
f(:,:,i) = makePDF(post_samps(:,i,:), xin(i,:));
end

figure;
alpha_lvl = 0.9; % Define the alpha level
c = {'r', 'g', 'b', [0.9290 0.6940 0.1250]}; fz = 18;
label = {'$E[p_1]$', '$V[p_1]$', '$p_{2}$', '$E[p_4]$', '$V[p_4]$', '$E[p_5]$', '$V[p_5]$', '$\rho$'};
for i = 1:size(f,3)
subplot(2,4,i)
hold on; box on; grid on;
for j = 1:4
plot(xin(i,:), f(j,:,i), 'color', c{j}, 'linewidth', 2); set(gca, 'Fontsize', fz);
xlabel(label{i}, 'Interpreter', 'latex'); ylabel('Normalised PDF value')
plot(xin(i,:), alpha_lvl.*ones(1,Nplot), 'k--', 'linewidth', 2, 'HandleVisibility', 'off')
end
end
legend('Euclidean', 'Bhattacharrya', 'Bray-Curtis', '1-Wasserstein', 'linewidth', 2)

%% Compute the result of the epistemic interval at the determined alpha-level:

interval_res = zeros(size(post_samps,2), 2, size(post_samps,3));
interval_samps = zeros(size(post_samps,2), 2, size(post_samps,3));
for i = 1:size(post_samps,3)
for j = 1:size(post_samps,2)
idx = knnsearch(f(i,:,j)', alpha_lvl, 'K', 2); interval_res(j,:,i) = [min(xin(j, idx)), max(xin(j, idx))];
interval_samps(j,:,i) = prctile(post_samps(:,j,i), [45, 55]);
end
end