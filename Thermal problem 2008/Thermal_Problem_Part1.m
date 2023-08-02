%% Material characterisation exercise:
% 
% For the first part of the challenge, the objective is to characterise the
% variability of the thermal properties of the given specimen based on a
% given limited data set.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; riskCalc4Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the data set:

temperature_data = [20, 250, 500, 750, 1000];

k_data = [0.0496, 0.0628, 0.0602, 0.0657 0.0631; 0.0530, 0.0620, 0.0546 0.0713, 0.0796; ...
          0.0493, 0.0537, 0.0638, 0.0694, 0.0692; 0.0455, 0.0561, 0.0614, 0.0732, 0.0739];
pd_k = fitdist(k_data(:),'Normal')

% Normal distribution
% mu =    0.06187   [0.0575502, 0.0661898]
% sigma = 0.00923011   [0.00701941, 0.0134812]

v_data = [3.76E+05, 3.87E+05, 4.52E+05, 4.68E+05, 4.19E+05; 3.38E+05, 4.69E+05, 4.10E+05 4.24E+05, 4.38E+05; ...
          3.50E+05, 4.19E+05, 4.02E+05, 3.72E+05, 3.45E+05; 4.13E+05 4.28E+05, 3.94E+05, 3.46E+05, 3.95E+05];
pd_v = fitdist(v_data(:),'Normal')

% Normal distribution
% mu =  402250   [383758, 420742]
% sigma = 39511.3   [30048, 57709.1]

%% Illustrating the data for Thermal conductivity & Volumetric heat capacity:

f = 20; nbins = 4;

figure; 
subplot(2,2,1)
hold on; box on; grid on;
[y1,x1] = ecdf(k_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$k$ $[W/{m} \cdot {^oC}]$', 'Interpreter', 'latex'); xlim([0.04, 0.09]); ylabel('ECDF value'); set(gca, 'Fontsize', f)

subplot(2,2,2)
hold on; box on; grid on;
[y1,x1] = ecdf(v_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$\rho$ $C_p$ $[J/{m^3} \cdot {^oC}]$', 'Interpreter', 'latex'); xlim([3E+05, 5E+05]); ylabel('ECDF value'); set(gca, 'Fontsize', f)

subplot(2,2,3)
hold on; box on; grid on;
histogram(k_data(:), nbins)
xlabel('$k$ $[W/{m} \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('Count'); ylim([0, 10]); set(gca, 'Fontsize', f)

subplot(2,2,4)
hold on; box on; grid on;
histogram(v_data(:), nbins)
xlabel('$\rho$ $C_p$ $[J/{m^3} \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('Count'); ylim([0, 10]); set(gca, 'Fontsize', f)

%% Perform Bayesian model updating via ABC on the distribution of the Thermal conductivity & Volumetric heat capacity:

% The distribution model used is the Normal distribution

% Priors for the Thermal conductivity & Volumetric heat capacity:
Nsim = 30;                                    % The number of stochastic model output realisations
N = 1000;                                     % Define the sample size from posterior
Nbatch = 1;                                % Define the number of sampling batches
TEMCMC_k = cell(Nbatch,4); timeTEMCMC_k = zeros(Nbatch,4); 
TEMCMC_rho = cell(Nbatch,4); timeTEMCMC_rho = zeros(Nbatch,4); 

bounds_k = [0.02, 0.08; 0.001, 0.1];
priorpdf_k = @(x) unifpdf(x(:,1), bounds_k(1,1), bounds_k(1,2)) .*...
                  unifpdf(x(:,2), bounds_k(2,1), bounds_k(2,2));  % Prior for Thermal conductivity hyper parameters;
priorrnd_k = @(N) [unifrnd(bounds_k(1,1), bounds_k(1,2), N, 1), unifrnd(bounds_k(2,1), bounds_k(2,2), N, 1)]; 
model_beta = @(theta,Nsim) normrnd(theta(:,1), theta(:,2), Nsim, 1);
logL_k2 = @(theta, dm, width) loglikelihood_func(theta, k_data(:), dm, width, Nsim, model_beta);  % dm is the index for the type of distance metric

for i = 1:4
width_k = [0.001, 0.06, 0.06, 0.0015];          % The width factor of the approximate likelihood function
for j = 1:Nbatch
sprintf('Now running simulation batch = %d', i)
logL = @(theta) logL_k2(theta, i, width_k(i));
tic;
TEMCMC_k{j,i} = TEMCMCsampler('nsamples', N, 'loglikelihood', logL, 'priorpdf', priorpdf_k, 'priorrnd', priorrnd_k);
timeTEMCMC_k(j,i) = toc;
end
end

bounds_rho = [2e5, 8e5; 1e4, 8e4];
priorpdf_rho = @(x) unifpdf(x(:,1), bounds_rho(1,1), bounds_rho(1,2)) .*...
                    unifpdf(x(:,2), bounds_rho(2,1), bounds_rho(2,2));  % Prior for Volumetric heat capacity hyper parameters;
priorrnd_rho = @(N) [unifrnd(bounds_rho(1,1), bounds_rho(1,2), N, 1), unifrnd(bounds_rho(2,1), bounds_rho(2,2), N, 1)];
model_beta = @(theta,Nsim) normrnd(theta(:,1), theta(:,2), Nsim, 1);
logL_rho2 = @(theta, dm, width) loglikelihood_func(theta, v_data(:), dm, width, Nsim, model_beta); % dm is the index for the type of distance metric

for i = 1:4
width_rho = [0.025, 0.07e-5, 0.036e-5, 0.04]*1e5;   % The width factor of the approximate likelihood function
for j = 1:Nbatch
sprintf('Now running simulation batch = %d', i);
logL = @(theta) logL_rho2(theta, i, width_rho(i));
tic;
TEMCMC_rho{j,i} = TEMCMCsampler('nsamples', N, 'loglikelihood', logL, 'priorpdf', priorpdf_rho, 'priorrnd', priorrnd_rho);
timeTEMCMC_rho(j,i) = toc;
end
end

%% Posterior analysis:

k_post = {TEMCMC_k{1,1}, TEMCMC_k{1,2}, TEMCMC_k{1,3}, TEMCMC_k{1,4}};
rho_post = {TEMCMC_rho{1,1}, TEMCMC_rho{1,2}, TEMCMC_rho{1,3}, TEMCMC_rho{1,4}};

k_post_samps = zeros(N,2,4); rho_post_samps = zeros(N,2,4); 
for j = 1:4
k_TEMCMC = k_post{j}; k_post_samps(:,:,j) = k_TEMCMC.samples;
rho_TEMCMC = rho_post{j}; rho_post_samps(:,:,j) = rho_TEMCMC.samples;
end

figure; nbins = 30; f = 18;
label = {'$\mu_{k}$', '$\sigma_{k}$', '$\mu_{\rho C_p}$', '$\sigma_{\rho C_p}$'};

subplot(2,2,1)
hold on; box on; grid on;
for j =1:4
histogram(k_post_samps(:,1,j), nbins); xline(0.06187, 'k--', 'linewidth', 2)
end
xlabel(label{1}, 'Interpreter', 'latex'); ylabel('Count'); set(gca, 'Fontsize', f);

subplot(2,2,2)
hold on; box on; grid on;
for j =1:4
histogram(k_post_samps(:,2,j), nbins); xline(0.00923011, 'k--', 'linewidth', 2)
end
xlabel(label{2}, 'Interpreter', 'latex'); ylabel('Count'); set(gca, 'Fontsize', f);

subplot(2,2,3)
hold on; box on; grid on;
for j =1:4
histogram(rho_post_samps(:,1,j), nbins); xline(402250, 'k--', 'linewidth', 2)
end
xlabel(label{3}, 'Interpreter', 'latex'); ylabel('Count'); set(gca, 'Fontsize', f);

subplot(2,2,4)
hold on; box on; grid on;
for j =1:4
histogram(rho_post_samps(:,2,j), nbins); xline(39511.3, 'k--', 'linewidth', 2) 
end
xlabel(label{4}, 'Interpreter', 'latex'); ylabel('Count'); set(gca, 'Fontsize', f);

%% Plot the approximated Fuzzy-set function via KDE:
Nplot = 500;
xin = [linspace(0.02, 0.08, Nplot); linspace(0.02, 0.08, Nplot); linspace(0.02, 0.08, Nplot); linspace(0.02, 0.08, Nplot);...
       linspace(0, 0.06, Nplot); linspace(0, 0.06, Nplot); linspace(0, 0.06, Nplot); linspace(0, 0.06, Nplot);... 
       linspace(3.6e5, 4.6e5, Nplot); linspace(3.6e5, 4.6e5, Nplot); linspace(3.6e5, 4.6e5, Nplot); linspace(3.6e5, 4.6e5, Nplot);...
       linspace(0, 8e4, Nplot); linspace(0, 8e4, Nplot); linspace(0, 8e4, Nplot); linspace(0, 8e4, Nplot)];
post_samps = [k_post_samps(:,1,1), k_post_samps(:,1,2), k_post_samps(:,1,3), k_post_samps(:,1,4), ...
              k_post_samps(:,2,1), k_post_samps(:,2,2), k_post_samps(:,2,3), k_post_samps(:,2,4), ...
              rho_post_samps(:,1,1), rho_post_samps(:,1,2), rho_post_samps(:,1,3), rho_post_samps(:,1,4), ...
              rho_post_samps(:,2,1), rho_post_samps(:,2,2), rho_post_samps(:,2,3), rho_post_samps(:,2,4)];

f = zeros(size(post_samps,3), Nplot, 16);
for i = 1:size(post_samps,2)
f(:,:,i) = makePDF(post_samps(:,i,:), xin(i,:));
end

figure;
alpha_lvl = 0.9; % Define the alpha level
c = {'r', 'g', 'b', [0.9290 0.6940 0.1250]}; fz = 18;
label = {'$\mu_{k}$ $[W/{m} \cdot {^oC}]$', '$\sigma_{k}$ $[W/{m} \cdot {^oC}]$', ...
         '$\mu_{\rho C_p}$ $[J/{m^3} \cdot {^oC}]$', '$\sigma_{\rho C_p}$ $[J/{m^3} \cdot {^oC}]$'};

subplot(2,2,1)
hold on; box on; grid on;
for j = 1:4
plot(xin(j,:), f(1,:,j), 'color', c{j}, 'linewidth', 2); set(gca, 'Fontsize', fz);
end
xlabel(label{1}, 'Interpreter', 'latex'); ylabel('Normalised PDF value'); xlim([0.04, 0.08])
plot(xin(4,:), alpha_lvl.*ones(1,Nplot), 'k--', 'linewidth', 2, 'HandleVisibility', 'off')

subplot(2,2,2)
hold on; box on; grid on;
for j = 1:4
plot(xin(j+4,:), f(1,:,j+4), 'color', c{j}, 'linewidth', 2); set(gca, 'Fontsize', fz);
end
xlabel(label{2}, 'Interpreter', 'latex'); ylabel('Normalised PDF value')
plot(xin(8,:), alpha_lvl.*ones(1,Nplot), 'k--', 'linewidth', 2, 'HandleVisibility', 'off')

subplot(2,2,3)
hold on; box on; grid on;
for j = 1:4
plot(xin(j+8,:), f(1,:,j+8), 'color', c{j}, 'linewidth', 2); set(gca, 'Fontsize', fz);
end
xlabel(label{3}, 'Interpreter', 'latex'); ylabel('Normalised PDF value')
plot(xin(12,:), alpha_lvl.*ones(1,Nplot), 'k--', 'linewidth', 2, 'HandleVisibility', 'off')

subplot(2,2,4)
hold on; box on; grid on;
for j = 1:4
plot(xin(j+12,:), f(1,:,j+12), 'color', c{j}, 'linewidth', 2); set(gca, 'Fontsize', fz);
end
xlabel(label{4}, 'Interpreter', 'latex'); ylabel('Normalised PDF value'); 
plot(xin(16,:), alpha_lvl.*ones(1,Nplot), 'k--', 'linewidth', 2, 'HandleVisibility', 'off')
legend('Euclidean', 'Bhattacharyya', 'Bray-Curtis', '1-Wasserstein', 'linewidth', 2)

%% Initiate the RiskCalc package

riskCalc4Matlab

%% Plot the P-boxes for the aleatory variables:
k_pbox_e = normal([6.028e-2, 6.353e-2],[1.323e-3, 1.154e-2]);  rho_pbox_e = normal([399479, 403487],[15390.8, 24368.7]);
k_pbox_b = normal([5.824e-2, 6.112e-2],[8.297e-3, 1.118e-2]);  rho_pbox_b = normal([389459, 405892],[47775.6, 55631.3]);
k_pbox_bc = normal([6.016e-2, 6.497e-2],[8.778e-3, 1.347e-2]);  rho_pbox_bc = normal([399479, 407695],[34468.9, 44088.2]); 
k_pbox_1w = normal([6.148e-2, 6.257e-2],[8.417e-3, 9.379e-3]);  rho_pbox_1w = normal([401483, 407495],[32384.8, 38316.6]);


figure; f = 18; c = {'r', 'g', 'b', [0.9290 0.6940 0.1250]};
subplot(2,4,1)
hold on; box on; grid on;
lb_pboxR = k_pbox_e.u; ub_pboxR = k_pbox_e.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(k_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$k$ $[W/m \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([0, 0.1])
title('P-box (Euclidean)')

subplot(2,4,2)
hold on; box on; grid on;
lb_pboxR = k_pbox_b.u; ub_pboxR = k_pbox_b.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(k_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$k$ $[W/m \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([0, 0.1])
title('P-box (Bhattacharyya)')

subplot(2,4,3)
hold on; box on; grid on;
lb_pboxR = k_pbox_bc.u; ub_pboxR = k_pbox_bc.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(k_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$k$ $[W/m \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([0, 0.1])
title('P-box (Bray-Curtis)')

subplot(2,4,4)
hold on; box on; grid on;
lb_pboxR = k_pbox_1w.u; ub_pboxR = k_pbox_1w.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(k_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$k$ $[W/m \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([0, 0.1])
title('P-box (1-Wasserstein)')

subplot(2,4,5)
hold on; box on; grid on;
lb_pboxR = rho_pbox_e.u; ub_pboxR = rho_pbox_e.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(v_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$\rho$ $C_p$ $[J/{m^3} \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([2e5, 6e5])
title('P-box (Euclidean)')

subplot(2,4,6)
hold on; box on; grid on;
lb_pboxR = rho_pbox_b.u; ub_pboxR = rho_pbox_b.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(v_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$\rho$ $C_p$ $[J/{m^3} \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([2e5, 6e5])
title('P-box (Bhattacharyya)')

subplot(2,4,7)
hold on; box on; grid on;
lb_pboxR = rho_pbox_bc.u; ub_pboxR = rho_pbox_bc.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(v_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$\rho$ $C_p$ $[J/{m^3} \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([2e5, 6e5])
title('P-box (Bray-Curtis)')

subplot(2,4,8)
hold on; box on; grid on;clc
lb_pboxR = rho_pbox_1w.u; ub_pboxR = rho_pbox_1w.d;
[yp1, xp1] = ecdf(lb_pboxR); stairs(xp1, yp1, 'color', c{3}, 'linewidth', 2);
[yp2, xp2] = ecdf(ub_pboxR); stairs(xp2, yp2, 'color', c{3}, 'linewidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(v_data(:)); stairs(x1, y1, 'r', 'LineWidth', 2);
xlabel('$\rho$ $C_p$ $[J/{m^3} \cdot {^oC}]$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', f); xlim([2e5, 6e5])
title('P-box (1-Wasserstein)')

%% Save data:

save('Thermal_Problem_Part3')