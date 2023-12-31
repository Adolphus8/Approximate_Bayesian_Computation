%% Illustrative Example:
%
% Examples to demonstrate the implementation of the respective Stochastic
% distance metrics and compare the results against a defined pair of random
% sample sets:
%
clc; clear; 
%% Define the parameters:
Nsamps = 10000; % No. of samples

x1 = normrnd(2, 1, Nsamps, 1);  % Random variable x1 ~ N(2,1)
x2 = normrnd(8, 2, Nsamps, 1);  % Random variable x2 ~ N(8,2)
x3 = normrnd(15, 3, Nsamps, 1); % Random variable x3 ~ N(15,3)
x4 = normrnd(15, 1, Nsamps, 1); % Random variable x4 ~ N(15,1)
x5 = normrnd(25, 1, Nsamps, 1); % Random variable x4 ~ N(15,1)

% Plot their combined histograms:
figure;
subplot(2,1,1)
hold on; box on; grid on;
nbins = 50; % No. of bins for the histograms
histogram(x1, nbins, 'FaceColor','b','Normalization','pdf')
histogram(x2, nbins, 'FaceColor','g','Normalization','pdf')
histogram(x3, nbins, 'FaceColor','r','Normalization','pdf')
histogram(x4, nbins, 'FaceColor',[0.9290 0.6940 0.1250],'Normalization','pdf')
histogram(x5, nbins, 'FaceColor','c','Normalization','pdf')
legend('$\mathbf{x}_1$', '$\mathbf{x}_2$', '$\mathbf{x}_3$', '$\mathbf{x}_4$', '$\mathbf{x}_5$', 'linewidth', 2, 'Interpreter', 'latex')
xlabel('Sample values'); ylabel('PDF value'); set(gca, 'Fontsize', 20); xlim([-5, 30])

subplot(2,1,2)
hold on; box on; grid on;
[yo, xo] = ecdf(x1); stairs(xo, yo, 'b', 'linewidth', 2);
[yo, xo] = ecdf(x2); stairs(xo, yo, 'g', 'linewidth', 2);
[yo, xo] = ecdf(x3); stairs(xo, yo, 'r', 'linewidth', 2);
[yo, xo] = ecdf(x4); stairs(xo, yo, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2);
[yo, xo] = ecdf(x5); stairs(xo, yo, 'c', 'linewidth', 2);
legend('$\mathbf{x}_1$', '$\mathbf{x}_2$', '$\mathbf{x}_3$', '$\mathbf{x}_4$', '$\mathbf{x}_5$', 'linewidth', 2, 'Interpreter', 'latex')
xlabel('Sample values'); ylabel('CDF value'); set(gca, 'Fontsize', 20); xlim([-5, 30])

%% Comapre the results by each distance metric:

samps = [x1, x2, x3, x4, x5];
EDMe_Res = zeros(size(samps,2), size(samps, 2)); % Results for Euclidean distance metric
BDMe_Res = zeros(size(samps,2), size(samps, 2)); % Results for Bhattachryya distance metric
BCMe_Res = zeros(size(samps,2), size(samps, 2)); % Results for Bray-Curtis distance metric
WMe_Res = zeros(size(samps,2), size(samps, 2));  % Results for 1-Wasserstein distance metric

timeED = zeros(size(samps,2), size(samps, 2));
timeBD = zeros(size(samps,2), size(samps, 2));
timeBC = zeros(size(samps,2), size(samps, 2));
timeW = zeros(size(samps,2), size(samps, 2));

for i = 1:size(samps,2)
for j = 1:size(samps,2)
tic; EDMe_Res(i,j) = EDMe(samps(:,i), samps(:,j)); timeED(i,j) = toc;
tic; BDMe_Res(i,j) = BDMe(samps(:,i), samps(:,j),[]); timeBD(i,j) = toc;
tic; BCMe_Res(i,j) = BCMe(samps(:,i), samps(:,j)); timeBC(i,j) = toc;
tic; WMe_Res(i,j) = areaMe(samps(:,i), samps(:,j)); timeW(i,j) = toc;
end
end

%% Obtain the statistics of the computation times by each distance metric:

stats_ED = [mean(timeED(:)), std(timeED(:))]; stats_BD = [mean(timeBD(:)), std(timeBD(:))];
stats_BC = [mean(timeBC(:)), std(timeBC(:))]; stats_W = [mean(timeW(:)), std(timeW(:))];

%% Save the data:
save('Illustrative_example')