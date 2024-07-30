function [output] = DLMC(bounds, Ne, Na, model)
%% The Double-Loop Monte Carlo function:
%
% Inputs:
% bounds: the dim x 2 matrix of input bounds;
% Ne:     the scalar value of the number of epistemic realizations;
% Na:     the scalar value of the number of aleatory realizations;
% model: the black-box function-handle of h1;
%
% output:
% output.samples:    the Na x Ne matrix of sample outputs;
% output.pbox:       the Na x 2 matrix of sample outputs for P-box;
% output.time:       the total time elapsed by the DLMC procedure;
%
% Note:
% theta_1 = E[p1]; - Mean of p1
% theta_2 = V[p1]; - Variance of p1
% theta_3 = delta; - Epistemic interval of p2
% theta_4 = E[p4]; - Mean of p4
% theta_5 = V[p4]; - Variance of p4
% theta_6 = E[p5]; - Mean of p5
% theta_7 = V[p5]; - Variance of p5
% theta_8 = rho;   - Correlation coefficient between p4 and p5
%
%% Define the epistemic hyper-rectangle:
tic;
epistemic_samps = [unifrnd(bounds(1,1), bounds(1,2), Ne, 1), unifrnd(bounds(2,1), bounds(2,2), Ne, 1),...
                   unifrnd(bounds(3,1), bounds(3,2), Ne, 1), unifrnd(bounds(4,1), bounds(4,2), Ne, 1),...
                   unifrnd(bounds(5,1), bounds(5,2), Ne, 1), unifrnd(bounds(6,1), bounds(6,2), Ne, 1),...
                   unifrnd(bounds(7,1), bounds(7,2), Ne, 1), unifrnd(bounds(8,1), bounds(8,2), Ne, 1)]; % Ne x 8 matrix


%% Define the models:

% Beta shape parameter models:
alpha = @(x) ((x(:,1)./x(:,2)).*(1 - x(:,1)) - 1) .* x(:,1);
beta  = @(x) ((x(:,1)./x(:,2)).*(1 - x(:,1)) - 1) .* (1 - x(:,1));

% Multivariate Normal shape parameter models:
mu =  @(x) [x(:,1), x(:,2)]; 
cov = @(x) [x(:,1), (x(:,3).*sqrt(x(:,1)).*sqrt(x(:,2))); (x(:,3).*sqrt(x(:,1)).*sqrt(x(:,2))), x(:,2)];

% Define the black-box model inputs:
p1 = @(x) betarnd(alpha(x), beta(x), Na, 1);
p2 = @(x) x .* ones(Na,1);

%% Computation procedure:

samples = zeros(Na,Ne); pbox = zeros(Na,2);
for i = 1:Ne
input = [p1(epistemic_samps(i,[1:2])), p2(epistemic_samps(i,3)), unifrnd(0, 1, Na, 1), mvnrnd(mu(epistemic_samps(i, [4,6])), cov(epistemic_samps(i, [5,7,8])), Na)];
samples(:,i) = model(input);
end
samples = sort(samples);

for i = 1:Na
pbox(i,:) = [min(samples(i,:)), max(samples(i,:))]; 
end

timeDLMC = toc;
sprintf('Total time elapsed for the DLMC procedure is = %3f', timeDLMC)

%% Generate the outputs:
output.samples = samples;
output.pbox = pbox;
output.time = timeDLMC;

end

