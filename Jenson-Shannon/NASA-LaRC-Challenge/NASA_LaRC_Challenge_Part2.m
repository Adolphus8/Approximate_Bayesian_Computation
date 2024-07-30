%% NASA LaRC Challenge 2014: Subproblem A
%
% The codes provided here describe the approach to analyse the results from
% the TEMCMC sampling approach via ABC.

clc; clear;
%% Load the training and validation data sets:
load('x1samples1.mat'); load('x1samples2.mat'); load('NASA_LaRC_Challenge_Part1.mat', 'TEMCMC', 'timeTEMCMC')

%% Analysis of samples:

post_samps = TEMCMC.samples;

% Plot the Histogram representation of the posteriors given each distance metric:
figure;
label = {'$E[p_1]$', '$V[p_1]$', '$p_{2}$', '$E[p_4]$', '$V[p_4]$', '$E[p_5]$', '$V[p_5]$', '$\rho$'};
for i = 1:8
subplot(2,4,i)
hold on; box on; grid on;
nbins = 10;
histogram(post_samps(:,i), nbins, 'Normalization','pdf')
set(gca, 'Fontsize', 18); xlabel(label{i}, 'Interpreter', 'latex'); ylabel('Count')
end

% Plot the approximated Fuzzy-set function via KDE:
Nplot = 10000;
xin = [linspace(0.6, 0.8, Nplot); linspace(0.02, 0.04, Nplot); linspace(0, 1, Nplot); linspace(-5, 5, Nplot);... 
       linspace(1/400, 4, Nplot); linspace(-5, 5, Nplot); linspace(1/400, 4, Nplot); linspace(-1, 1, Nplot)];
f = zeros(Nplot, size(post_samps,2));

for i = 1:size(post_samps,2)
f(:,i) = makePDF(post_samps(:,i,:), xin(i,:));
end

figure;
alpha_lvl = 0.9; % Define the alpha level
fz = 18;
label = {'$E[p_1]$', '$V[p_1]$', '$p_{2}$', '$E[p_4]$', '$V[p_4]$', '$E[p_5]$', '$V[p_5]$', '$\rho$'};
for i = 1:size(f,2)
subplot(2,4,i)
hold on; box on; grid on;
plot(xin(i,:)', f(:,i), 'b', 'linewidth', 2); set(gca, 'Fontsize', fz);
xlabel(label{i}, 'Interpreter', 'latex'); ylabel('Normalised PDF value'); ylim([0, 1])
plot(xin(i,:), alpha_lvl.*ones(1,Nplot), 'k--', 'linewidth', 2, 'HandleVisibility', 'off')
end

%% Compute the result of the epistemic interval at the determined alpha-level:

interval_res = zeros(size(post_samps,2), 2); 
for i = 1:size(post_samps,2)
idx = knnsearch(f(:,i), alpha_lvl, 'K', 2); interval_res(i,:) = [min(xin(i, idx)), max(xin(i, idx))];
end

%% Double-Loop Monte Carlo for JS-divergence:

Ne = 1000;                         % No. of epistemic realizations to generate from the epistemic space
Na = 1000;                         % No. of aleatory realizations from the stochastic black-box model
blackbox_model = @(p) p_to_x1(p); % The Black-box model describing x1 = h1(p1,p2,p3,p4,p5)

bounds_JS = interval_res;
out_JS = DLMC(bounds_JS, Ne, Na, blackbox_model);

pbox_JS = out_JS.pbox;            % Compute the area of the P-box
pbox_area = areaMe(pbox_JS(:,1), pbox_JS(:,2));
sprintf('The area of the P-box obtained via the JS-divergence distance metric is = %4f', pbox_area)

%% Verify and Validate the model output:

figure;
hold on; box on; grid on;
[y1,x1] = ecdf(x1sams1); stairs(x1, y1, 'b', 'LineWidth', 2); [y1,x1] = ecdf(x1sams2); stairs(x1, y1, 'r', 'LineWidth', 2);
[y1,x1] = ecdf(combined_data); stairs(x1, y1, 'g', 'LineWidth', 2);
[y1,x1] = ecdf(pbox_JS(:,1)); [y2,x2] = ecdf(pbox_JS(:,2)); 
stairs(x1, y1, 'k', 'LineWidth', 2); stairs(x2, y2, 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
plot([min(x1), min(x2)], [0,0], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); plot([max(x1), max(x2)], [1,1], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
legend('Training data', 'Validation data', 'Combined data', 'P-box', 'linewidth', 2, 'location', 'southeast'); 
xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20); xlim([0,1])

%% Analyse the Verification, Validation, and Combined scores:
samps = out_JS.samples; val_stats = zeros(2,3);

for k = 1:3

if k == 1
ref_data = x1sams1;
elseif k == 2
ref_data = x1sams2;
elseif k == 3
ref_data = combined_data;
end

val_score = zeros(Ne,1);
for j = 1:Ne
val_score(j) = areaMe(samps(:,j), ref_data);
end
val_stats(:,k) = [mean(val_score), std(val_score)]';

end

%% Save the data:
save('NASA_LaRC_Challenge_Part2')
