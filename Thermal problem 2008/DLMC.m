function [output] = DLMC(bounds, Ne, Na)
%% The Double-Loop Monte Carlo function:
%
% Inputs:
% bounds: the dim x 2 matrix of input bounds;
% Ne:     the scalar value of the number of epistemic realizations;
% Na:     the scalar value of the number of aleatory realizations;
% model: the black-box function-handle of h1;
%
% output:
% output.samples: the Na x Ne matrix of sample outputs;
% output.pbox:    the Na x 2 matrix of sample outputs for P-box;
% output.time:    the total time elapsed by the DLMC procedure;
%
%% Define the epistemic hyper-rectangle:
tic;
epistemic_samps = [unifrnd(bounds(1,1), bounds(1,2), Ne, 1),...
                   unifrnd(bounds(2,1), bounds(2,2), Ne, 1),...
                   unifrnd(bounds(3,1), bounds(3,2), Ne, 1),...
                   unifrnd(bounds(4,1), bounds(4,2), Ne, 1)]; % Ne x 4 matrix

%% Computation procedure:
tin = [0:10:1000]';
samples = zeros(length(tin), Na, Ne); pbox = zeros(Na, 2, length(tin));

for i = 1:Ne
model_input = [normrnd(epistemic_samps(i,1), epistemic_samps(i,2), Na, 1), normrnd(epistemic_samps(i,3), epistemic_samps(i,4), Na, 1)];

for j = 1:Na
if model_input(j,1) < 0
model_input(j,1) = 0.00001;
end

samps = Temperature_model(tin, model_input(j,1), model_input(j,2));
samples(:,j,i) = sort(samps);
end
end

for i = 1:Na
mat_sim = permute(samples, [2,3,1]);
for k = 1:length(tin)
pbox(i,:,k) = [min(mat_sim(i,:,k)), max(mat_sim(i,:,k))]; 
end
end

timeDLMC = toc;
sprintf('Total time elapsed for the DLMC procedure is = %3f', timeDLMC)

%% Generate the outputs:
output.samples = samples;
output.pbox = pbox;
output.time = timeDLMC;

end

%% Define the model function:
function [output] = Temperature_model(t,k,v)
%% Function-handle for the Temperature model of the Sandia Thermal problem:
% See Eq. (2) of the Sandia Thermal problem question paper.
%
% Reference:
% K. J. Dowding, M. Pilch, and R. G. Hills (2008). Formulation of the thermal problem. 
% Computer Methods in Applied Mechanics and Engineering, 197(29-32), 2385–2389. 
% doi: 10.1016/j.cma.2007.09.029 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:
% x: Scalar input of x-coordinates along the thickness L of the slab [m];
% t: Time parameter [s];
% k: Thermal conductivity of the slab [W/m.deg]
% v: Volumetric heat capactiy [J/m^3.deg]
% q: Heat flux [W/m^2]
% L: Thickness of slab [m] 
%
% Output:
% output: Output temperature [deg];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the fixed parameter:

Ti = 25;     % Initial temperature of slab [deg]
x = 0;       % Position on slab [m]
q = 1000;    % Heat flux [W/m^2]
L = 2.54e-2; % Thickness of slab [m]     

%% Define the model:

% Define the summation term:
f = @(n,x,t) ((1/n).^2).*(exp(-(n.^2)*(pi.^2)*(k/v).*t./L.^2))*cos(n.*pi.*x./L);

% Define the output:
if t == 0
output = Ti;
else
output = Ti + (q*L/k)*(((k/v).*t./L^2) + (1/3) - (x./L) + 0.5*(x./L).^2 - ...
        (2/pi^2)*(f(1,x,t) + f(2,x,t) + f(3,x,t) + f(4,x,t) + f(5,x,t) + f(6,x,t)));
end
end
