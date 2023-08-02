%% Model Ensemble Validation exercise:
% 
% For the first part of the challenge, the objective is to characterise the
% variability of the thermal properties of the given specimen based on a
% given limited data set.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data:

load('Thermal_Problem_Part3')

%% Define the ensemble validation data set:

t = [0:100:1000]';
T_data = [25.0, 99.5, 130.7, 154.4, 174.3, 191.7, 207.3, 221.7, 235.0, 247.6, 259.3; ...
          25.0, 106.6, 140.4, 165.9, 187.2, 205.8, 222.4, 237.6, 251.7, 264.9, 277.4; ...
          25.0, 96.2, 126.1, 148.7, 167.7, 184.3, 199.3, 213.0, 225.7, 237.6, 248.9; ...
          25.0, 101.3, 133.1, 157.2, 177.2, 194.8, 210.6, 225.0, 238.4, 251.0, 262.9]';

figure; col = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
subplot(3,1,1)
hold on; box on; grid on;
for i = 1:4
plot(t, T_data(:,i), '--', 'LineWidth', 2, 'Color',col{i}, 'handlevisibility', 'off') 
plot(t, T_data(:,i), 's', 'MarkerSize', 10, 'MarkerEdgeColor',col{i}, 'MarkerFaceColor', col{i}) 
end
legend('Experiment #1', 'Experiment #2', 'Experiment #3', 'Experiment #4', 'linewidth', 2, 'location', 'southeast')
xlabel('$t$ $[s]$', 'Interpreter', 'latex'); ylabel('$T(x = 0 m)$ $[^o C]$', 'Interpreter', 'latex'); set(gca, 'Fontsize', 18)

for j = 1:10
title_cell = {'$T$ $[^o C]$ $|$ $t=100s$', '$T$ $[^o C]$ $|$ $t=200s$', '$T$ $[^o C]$ $|$ $t=300s$', '$T$ $[^o C]$ $|$ $t=400s$',...
              '$T$ $[^o C]$ $|$ $t=500s$', '$T$ $[^o C]$ $|$ $t=600s$', '$T$ $[^o C]$ $|$ $t=700s$', '$T$ $[^o C]$ $|$ $t=800s$',...
              '$T$ $[^o C]$ $|$ $t=900s$', '$T$ $[^o C]$ $|$ $t=1000s$'};
subplot(3,5,j+5)
hold on; box on; grid on;
[y1,x1] = ecdf(T_data(j+1,:)); stairs(x1, y1, 'k', 'LineWidth', 2); 
xlabel(title_cell{j}, 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 18); xlim([50, 300]); 
end

%% Define the parameters and model:

Ne = 500;                        % No. of epistemic realizations to generate from the epistemic space
Na = 500;                        % No. of aleatory realizations from the Temperature model

%% Double-Loop Monte Carlo for Euclidean distance metric:

bounds_ED = [6.028e-2, 6.353e-2;...
             1.323e-3, 1.154e-2;...
             399479, 403487;...
             15390.8, 24368.7];

out_ED = DLMC(bounds_ED, Ne, Na);

%% Double-Loop Monte Carlo for Bhattacharyya distance metric:

bounds_BD = [5.824e-2, 6.112e-2;...
             8.297e-3, 1.118e-2;...
             389459, 405892;...
             47775.6, 55631.3];

out_BD = DLMC(bounds_BD, Ne, Na);

%% Double-Loop Monte Carlo for Bray-Curtis distance metric:

bounds_BC = [6.016e-2, 6.497e-2;...
             8.778e-3, 1.347e-2;...
             399479, 407695;...
             34468.9, 44088.2];

out_BC = DLMC(bounds_BC, Ne, Na);

%% Double-Loop Monte Carlo for 1-Wasserstein distance metric:

bounds_1W = [6.148e-2, 6.257e-2;...
             8.417e-3, 9.379e-3;...
             401483, 407495;...
             32384.8, 38316.6];

out_1W = DLMC(bounds_1W, Ne, Na);

%% Save the data:

save('Thermal_Problem_Part4', 'out_ED', 'out_BD', 'out_BC', 'out_1W')

%% Model validation:
samp_ED = out_ED.samples; samp_BD = out_BD.samples; samp_BC = out_BC.samples; samp_1W = out_1W.samples; 
pbox_ED = out_ED.pbox;    pbox_BD = out_BD.pbox;    pbox_BC = out_BC.pbox;    pbox_1W = out_1W.pbox; 
EDsamp = samp_ED(:,:); BDsamp = samp_BD(:,:); BCsamp = samp_BC(:,:); Wsamp = samp_1W(:,:); 

model_bounds = zeros(size(EDsamp,1),2,4);
for i = 1:size(EDsamp,1)
model_bounds(i,:,1) = [min(EDsamp(i,:)), max(EDsamp(i,:))]; model_bounds(i,:,2) = [min(BDsamp(i,:)), max(BDsamp(i,:))];
model_bounds(i,:,3) = [min(BCsamp(i,:)), max(BCsamp(i,:))]; model_bounds(i,:,4) = [min(Wsamp(i,:)), max(Wsamp(i,:))];
model_bounds(1,:,1) = [25, 25]; model_bounds(1,:,2) = [25, 25]; model_bounds(1,:,3) = [25, 25]; model_bounds(1,:,4) = [25, 25];
end

figure; col = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
c = {'r', 'g', 'b', [0.9290 0.6940 0.1250]}; title_cell = {'Euclidean', 'Bhattacharrya', 'Bray-Curtis', '1-Wasserstein'};
for j = 1:4
subplot(2,2,j)
hold on; box on; grid on;
for i = 1:4
plot(t, T_data(:,i), '--', 'LineWidth', 2, 'Color',col{i}, 'handlevisibility', 'off') 
plot(t, T_data(:,i), 's', 'MarkerSize', 10, 'MarkerEdgeColor',col{i}, 'MarkerFaceColor', col{i}) 
end
tin = [0:10:1000]';
plot(tin, model_bounds(:,1,j), 'color', 'k', 'LineWidth', 2); plot(tin, model_bounds(:,2,j), 'color', 'k', 'LineWidth', 2, 'handlevisibility', 'off') 
xlabel('$t$ $[s]$', 'Interpreter', 'latex'); ylabel('$T(x = 0 m)$ $[^o C]$', 'Interpreter', 'latex'); set(gca, 'Fontsize', 18); ylim([0, 500])
title(title_cell{j})
end
legend('Experiment #1', 'Experiment #2', 'Experiment #3', 'Experiment #4', 'linewidth', 2, 'location', 'southeast')

permute_ED = permute(samp_ED, [2,3,1]); permute_BD = permute(samp_BD, [2,3,1]); 
permute_BC = permute(samp_BC, [2,3,1]); permute_1W = permute(samp_1W, [2,3,1]); 

figure; 
c = {'r', 'g', 'b', [0.9290 0.6940 0.1250]};
for j = 1:10

idx = find(tin == j*100);
val_bounds = zeros(size(permute_ED,1),2,4);
for i = 1:size(permute_ED,1)
val_bounds(i,:,1) = [min(permute_ED(i,:,idx)), max(permute_ED(i,:,idx))]; val_bounds(i,:,2) = [min(permute_BD(i,:,idx)), max(permute_BD(i,:,idx))];
val_bounds(i,:,3) = [min(permute_BC(i,:,idx)), max(permute_BC(i,:,idx))]; val_bounds(i,:,4) = [min(permute_1W(i,:,idx)), max(permute_1W(i,:,idx))];
end
title_cell = {'$T$ $[^o C]$ $|$ $t=100s$', '$T$ $[^o C]$ $|$ $t=200s$', '$T$ $[^o C]$ $|$ $t=300s$', '$T$ $[^o C]$ $|$ $t=400s$',...
              '$T$ $[^o C]$ $|$ $t=500s$', '$T$ $[^o C]$ $|$ $t=600s$', '$T$ $[^o C]$ $|$ $t=700s$', '$T$ $[^o C]$ $|$ $t=800s$',...
              '$T$ $[^o C]$ $|$ $t=900s$', '$T$ $[^o C]$ $|$ $t=1000s$'};
subplot(3,4,j)
hold on; box on; grid on;
[y1,x1] = ecdf(T_data(j+1,:)); stairs(x1, y1, 'k', 'LineWidth', 2); 
for k = 1:4
[y1,x1] = ecdf(val_bounds(:,1,k)); stairs(x1, y1, 'color', c{k}, 'LineWidth', 2);
[y1,x1] = ecdf(val_bounds(:,2,k)); stairs(x1, y1, 'color', c{k}, 'LineWidth', 2, 'handlevisibility', 'off');
end
xlabel(title_cell{j}, 'Interpreter', 'latex'); ylabel('ECDF value'); set(gca, 'Fontsize', 18); xlim([0, 500]);
end
legend('Data', 'P-box (Euclidean)', 'P-box (Bhattacharyya)', 'P-box (Bray-Curtis)', 'P-box (1-Wasserstein)', 'linewidth', 2)

validation_stats = zeros(10,3,4);
for j = 1:10
idx = find(tin == j*100);
validation_stats(j,1,1) = areaMe(pbox_ED(:,1,idx), pbox_ED(:,2,idx)); validation_stats(j,1,2) = areaMe(pbox_BD(:,1,idx), pbox_BD(:,2,idx));
validation_stats(j,1,3) = areaMe(pbox_BC(:,1,idx), pbox_BC(:,2,idx)); validation_stats(j,1,4) = areaMe(pbox_1W(:,1,idx), pbox_1W(:,2,idx));

stats_ED = zeros(Ne,1); stats_BD = zeros(Ne,1); stats_BC = zeros(Ne,1); stats_1W = zeros(Ne,1);
for i = 1:Ne
stats_ED(i) = areaMe(permute_ED(:,i,idx), [T_data(j+1,:)]'); stats_BD(i) = areaMe(permute_BD(:,i,idx), [T_data(j+1,:)]');
stats_BC(i) = areaMe(permute_BC(:,i,idx), [T_data(j+1,:)]'); stats_1W(i) = areaMe(permute_1W(:,i,idx), [T_data(j+1,:)]');
end
validation_stats(j,2:3,1) = [mean(stats_ED), std(stats_ED)]; validation_stats(j,2:3,2) = [mean(stats_BD), std(stats_BD)];
validation_stats(j,2:3,3) = [mean(stats_BC), std(stats_BC)]; validation_stats(j,2:3,4) = [mean(stats_1W), std(stats_1W)];
end

figure; c = {'r', 'g', 'b', [0.9290 0.6940 0.1250]}; ylab = {'P-box area $[^o C]$', 'Mean $[^o C]$', 'Stdev. $[^o C]$'};
for jdx = 1:3
subplot(3,1,jdx)
hold on; box on; grid on;
for i = 1:4
plot((100:100:1000)', validation_stats(:,jdx,i), 'LineWidth', 2, 'color',c{i}, 'handlevisibility', 'off') 
plot((100:100:1000)', validation_stats(:,jdx,i), 's', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor',c{i}, 'MarkerFaceColor', c{i}) 
set(gca, 'Fontsize', 18); xlabel('$t$ $[s]$', 'Interpreter', 'latex'); ylabel(ylab{jdx}, 'Interpreter', 'latex');
end
end
legend('Euclidean', 'Bhattacharyya', 'Bray-Curtis', '1-Wasserstein', 'linewidth', 2)

