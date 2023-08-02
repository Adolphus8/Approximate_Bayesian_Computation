function [output] = Temperature_model(x,t,k,v,q,L)
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