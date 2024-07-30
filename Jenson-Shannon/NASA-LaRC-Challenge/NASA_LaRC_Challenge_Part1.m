%% NASA LaRC Challenge 2014: Subproblem A
%
% The codes provided here describe the approach to address the following
% subproblems: A1 and A2

clc; clear;
%% Load the training and validation data sets:
load('x1samples1.mat'); load('x1samples2.mat'); combined_data = [x1sams1; x1sams2];
nbins = 10; % No. of bins for the histograms

figure;
subplot(2,1,1)
hold on; box on; grid on;
[y1,x1] = ecdf(x1sams1); stairs(x1, y1, 'b', 'LineWidth', 2);
[y1,x1] = ecdf(x1sams2); stairs(x1, y1, 'r', 'LineWidth', 2);
[y1,x1] = ecdf(combined_data); stairs(x1, y1, 'g', 'LineWidth', 2);
legend('Training data', 'Validation data', 'Combined data', 'linewidth', 2, 'location', 'southeast'); xlim([0, 0.5])
xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20)

subplot(2,3,4); 
hold on; box on; grid on;
histogram(x1sams1, nbins)
legend('Training data', 'linewidth', 2, 'location', 'northeast'); xlim([0, 0.5]); ylim([0, 11])
xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20)

subplot(2,3,5); 
hold on; box on; grid on;
histogram(x1sams2, nbins, 'FaceColor', 'r')
legend('Validation data', 'linewidth', 2, 'location', 'northeast'); xlim([0, 0.5]); ylim([0, 11])
xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20)

subplot(2,3,6); 
hold on; box on; grid on;
histogram(combined_data, nbins, 'FaceColor', 'g')
legend('Combined data', 'linewidth', 2, 'location', 'northeast'); xlim([0, 0.5]); ylim([0, 11])
xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20)

%% Bayesian Model Updating set-up:
% There are in total 8 parameters to update:
% theta_1 = E[p1]; - Mean of p1
% theta_2 = V[p1]; - Variance of p1
% theta_3 = delta; - Epistemic interval of p2
% theta_4 = E[p4]; - Mean of p4
% theta_5 = V[p4]; - Variance of p4
% theta_6 = E[p5]; - Mean of p5
% theta_7 = V[p5]; - Variance of p5
% theta_8 = rho;   - Correlation coefficient between p4 and p5

% Define the Prior: 
lb = [3/5, 1/50, 0, -5, 1/400, -5, 1/400, -1]; % Vector of lower bounds
ub = [4/5, 1/25, 1, 5, 4, 5, 4, 1];            % Vector of upper bounds

priorpdf1 = @(theta_1) unifpdf(theta_1, lb(1), ub(1)); % Prior PDF for theta_1
priorpdf2 = @(theta_2) unifpdf(theta_2, lb(2), ub(2)); % Prior PDF for theta_2
priorpdf3 = @(theta_3) unifpdf(theta_3, lb(3), ub(3)); % Prior PDF for theta_3
priorpdf4 = @(theta_4) unifpdf(theta_4, lb(4), ub(4)); % Prior PDF for theta_4
priorpdf5 = @(theta_5) unifpdf(theta_5, lb(5), ub(5)); % Prior PDF for theta_5
priorpdf6 = @(theta_6) unifpdf(theta_6, lb(6), ub(6)); % Prior PDF for theta_6
priorpdf7 = @(theta_7) unifpdf(theta_7, lb(7), ub(7)); % Prior PDF for theta_7
priorpdf8 = @(theta_8) unifpdf(theta_8, lb(8), ub(8)); % Prior PDF for theta_8

priorpdf = @(theta) priorpdf1(theta(:,1)) .* priorpdf2(theta(:,2)) .* priorpdf3(theta(:,3)) .* priorpdf4(theta(:,4)) .* ... 
                    priorpdf5(theta(:,5)) .* priorpdf6(theta(:,6)) .* priorpdf7(theta(:,7)) .* priorpdf8(theta(:,8));
 
priorrnd = @(N) [unifrnd(lb(1), ub(1), N, 1), unifrnd(lb(2), ub(2), N, 1), unifrnd(lb(3), ub(3), N, 1), ...
                 unifrnd(lb(4), ub(4), N, 1), unifrnd(lb(5), ub(5), N, 1), unifrnd(lb(6), ub(6), N, 1), ...
                 unifrnd(lb(7), ub(7), N, 1), unifrnd(lb(8), ub(8), N, 1)];

% Define the Log-likelihood function:
blackbox_model = @(p) p_to_x1(p);      % The Black-box model describing x1 = h1(p1,p2,p3,p4,p5)
Nsim = 30;                             % The number of stochastic model output realisations
loglike = @(theta, width) loglikelihood(theta, x1sams1, blackbox_model, 5, width, Nsim); % dm is the index for the type of distance metric

% Define the sample size:
N = 1000;

%% Run the TEMCMC sampler

width = 4e-3; % Width parameter vector
logL = @(theta) loglike(theta, width);
tic;
TEMCMC = TEMCMCsampler('nsamples', N, 'loglikelihood', logL, 'priorpdf', priorpdf, 'priorrnd', priorrnd);
timeTEMCMC = toc;

%% Save the data:
save('NASA_LaRC_Challenge_Part1', 'TEMCMC', 'timeTEMCMC')
