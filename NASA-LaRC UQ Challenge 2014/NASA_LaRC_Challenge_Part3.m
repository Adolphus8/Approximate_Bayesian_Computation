%% NASA LaRC Challenge 2014: Subproblem A
%
% The codes provided here describe the approach to validate the results
% from the Bayesian model updating procedure via ABC

clc; clear;
%% Load the training and validation data sets:
load('x1samples1.mat'); load('x1samples2.mat');
combined_data = [x1sams1; x1sams2];

%% Define the parameters and model:

Ne = 500;                        % No. of epistemic realizations to generate from the epistemic space
Na = 500;                        % No. of aleatory realizations from the stochastic black-box model
blackbox_model = @(p) p_to_x1(p); % The Black-box model describing x1 = h1(p1,p2,p3,p4,p5)

%% Double-Loop Monte Carlo for Euclidean distance metric:

bounds_ED = [0.658, 0.677;...
             0.024, 0.029;...
             0.088, 0.234;...
            -2.174, -0.812;...
             3.239, 3.872;...
             1.413, 3.336;...
             0.635, 1.312;...
             0.663, 0.980];

out_ED = DLMC(bounds_ED, Ne, Na, blackbox_model);

%% Double-Loop Monte Carlo for Bhattacharyya distance metric:

bounds_BD = [0.649, 0.668;...
             0.031, 0.034;...
             0.078, 0.202;...
             1.092, 2.355;...
             0.091, 0.531;...
            -4.078, -3.116;...
             1.429, 1.965;...
             0.190, 0.511];

out_BD = DLMC(bounds_BD, Ne, Na, blackbox_model);

%% Double-Loop Monte Carlo for Bray-Curtis distance metric:

bounds_BC = [0.618, 0.631;...
             0.035, 0.037;...
             0.046, 0.144;...
             0.271, 1.373;...
             3.095, 3.495;...
            -2.856, -2.455;...
             0.155, 0.555;...
             0.760, 0.975];

out_BC = DLMC(bounds_BC, Ne, Na, blackbox_model);

%% Double-Loop Monte Carlo for 1-Wasserstein distance metric:

bounds_1W = [0.626, 0.636;...
             0.034, 0.036;...
             0.082, 0.194;...
            -1.513, -0.431;...
             1.501, 1.909;...
            -4.299, -3.657;...
             0.315, 0.836;...
            -0.940, -0.728];

out_1W = DLMC(bounds_1W, Ne, Na, blackbox_model);

%% Save the data:
save('NASA_LaRC_Challenge_Part3', 'out_ED', 'out_BD', 'out_BC', 'out_1W')

%% Compute the Area of P-box:
clc;

pbox_ED = out_ED.pbox; pbox_BD = out_BD.pbox; pbox_BC = out_BC.pbox; pbox_1W = out_1W.pbox;

pbox_area = [areaMe(pbox_ED(:,1), pbox_ED(:,2)), areaMe(pbox_BD(:,1), pbox_BD(:,2)), areaMe(pbox_BC(:,1), pbox_BC(:,2)), areaMe(pbox_1W(:,1), pbox_1W(:,2))];
sprintf('The area of the P-box obtained via the Eucidean distance metric is = %4f', pbox_area(1))
sprintf('The area of the P-box obtained via the Bhattacharyya distance metric is = %4f', pbox_area(2))
sprintf('The area of the P-box obtained via the Bray-Curtis distance metric is = %4f', pbox_area(3))
sprintf('The area of the P-box obtained via the 1-Wasserstein distance metric is = %4f', pbox_area(4))

%% Verify and Validate the model output:

figure;
subplot(2,2,1)
hold on; box on; grid on;
[y1,x1] = ecdf(x1sams1); stairs(x1, y1, 'b', 'LineWidth', 2); [y1,x1] = ecdf(x1sams2); stairs(x1, y1, 'r', 'LineWidth', 2);
[y1,x1] = ecdf(combined_data); stairs(x1, y1, 'g', 'LineWidth', 2);
[y1,x1] = ecdf(pbox_ED(:,1)); [y2,x2] = ecdf(pbox_ED(:,2)); 
stairs(x1, y1, 'k', 'LineWidth', 2); stairs(x2, y2, 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
plot([min(x1), min(x2)], [0,0], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); plot([max(x1), max(x2)], [1,1], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
legend('Training data', 'Validation data', 'Combined data', 'P-box (Euclidean)', 'linewidth', 2, 'location', 'southeast'); 
xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20); xlim([0,1])

subplot(2,2,2)
hold on; box on; grid on;
[y1,x1] = ecdf(x1sams1); stairs(x1, y1, 'b', 'LineWidth', 2, 'handlevisibility', 'off'); [y1,x1] = ecdf(x1sams2); stairs(x1, y1, 'r', 'LineWidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(combined_data); stairs(x1, y1, 'g', 'LineWidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(pbox_BD(:,1)); [y2,x2] = ecdf(pbox_BD(:,2)); 
stairs(x1, y1, 'k', 'LineWidth', 2); stairs(x2, y2, 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
plot([min(x1), min(x2)], [0,0], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); plot([max(x1), max(x2)], [1,1], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
legend('P-box (Bhattacharyya)', 'linewidth', 2, 'location', 'southeast'); xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20); xlim([0,1])

subplot(2,2,3)
hold on; box on; grid on;
[y1,x1] = ecdf(x1sams1); stairs(x1, y1, 'b', 'LineWidth', 2, 'handlevisibility', 'off'); [y1,x1] = ecdf(x1sams2); stairs(x1, y1, 'r', 'LineWidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(combined_data); stairs(x1, y1, 'g', 'LineWidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(pbox_BC(:,1)); [y2,x2] = ecdf(pbox_BC(:,2)); 
stairs(x1, y1, 'k', 'LineWidth', 2); stairs(x2, y2, 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
plot([min(x1), min(x2)], [0,0], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); plot([max(x1), max(x2)], [1,1], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
legend('P-box (Bray-Curtis)', 'linewidth', 2, 'location', 'southeast'); xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20); xlim([0,1])

subplot(2,2,4)
hold on; box on; grid on;
[y1,x1] = ecdf(x1sams1); stairs(x1, y1, 'b', 'LineWidth', 2, 'handlevisibility', 'off'); [y1,x1] = ecdf(x1sams2); stairs(x1, y1, 'r', 'LineWidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(combined_data); stairs(x1, y1, 'g', 'LineWidth', 2, 'handlevisibility', 'off');
[y1,x1] = ecdf(pbox_1W(:,1)); [y2,x2] = ecdf(pbox_1W(:,2)); 
stairs(x1, y1, 'k', 'LineWidth', 2); stairs(x2, y2, 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
plot([min(x1), min(x2)], [0,0], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); plot([max(x1), max(x2)], [1,1], 'k', 'LineWidth', 2, 'handlevisibility', 'off'); 
legend('P-box (1-Wasserstein)', 'linewidth', 2, 'location', 'southeast'); xlabel('$x_{1}$', 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 20); xlim([0,1])

%% Analyse the Verification, Validation, and Combined scores:
samps(:,:,1) = out_ED.samples; samps(:,:,2) = out_BD.samples; samps(:,:,3) = out_BC.samples; samps(:,:,4) = out_1W.samples; 

val_stats = zeros(2,4,3);

for k = 1:3

if k == 1
ref_data = x1sams1;
elseif k == 2
ref_data = x1sams2;
elseif k == 3
ref_data = combined_data;
end

for i = 1:4
val_score = zeros(Ne,1);
for j = 1:Ne
val_score(j) = areaMe(samps(:,j,i), ref_data);
end
val_stats(:,i,k) = [mean(val_score), std(val_score)]';
end
end