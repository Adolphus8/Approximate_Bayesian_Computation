function AM = areaMe(d1,d2,varargin)
%AREAME Computes the stochastic area metric between two data sets.
%   This version works with data sets of different sizes.
% 
%    %-------------------------------%
%    Author: Marco De Angelis  
%    _______ Created Oct 2020
%    _______ github.com/marcodeangelis
%    _______ University of Liverpool
%    %-------------------------------%
%
%   This version implements the presence of bounds. This can be useful for
%   example, when the area metric is needed only on a portion of the 
%   domain spanned by the data.
%
%   For example if your data span between -5 and 5, and you wanna know what
%   is the area metric for positive values only, you would do the following: 
%                                           
%   >>> AM = areaMe(d1,d2,[0,Inf]) % bounded [0,Inf]
% 
%   The following statement is equivalent to the unbounded area metric:
% 
%   >>> AM = areaMe(d1,d2,[-inf,Inf]) % unbounded [-Inf,Inf]
% 
%   For the unbounded case never use the above command, for efficiency use
%   this instead:
% 
%   >>> AM = areaMe(d1,d2).
% 
%%

n1 = length(d1);
n2 = length(d2);

bounded = false;
if nargin==3
bounded = true; % bounds have been provided
bounds = varargin{1}; % bounds on the area metric. This sends all of the cdf values to the left of the left-bound to zero; and all of the cdf values to the right of the right-bound to one.
end

[x1,y1]=ecdf_Lpool(d1); % Subroutine
[x1,y1]=stairs(x1,y1); % In-built matlab function
x1 = [x1(1);x1]; % this will make sure the data steps from 0 -> 1 on first datum
y1 = [0;y1]; % this will make sure the data steps from 0 -> 1 on first datum

[x2,y2]=ecdf_Lpool(d2); % Subroutine
[x2,y2]=stairs(x2,y2); % In-built matlab function
x2 = [x2(1);x2]; % this will make sure the data steps from 0 -> 1 on first datum
y2 = [0;y2]; % this will make sure the data steps from 0 -> 1 on first datum

assert(2*n1==length(y1),'The stairs function isnt being used correctly.')

d12 = sort([d1;d2]);

if bounded
    inside_bounds = d12>bounds(1) & d12<bounds(2);
    d12_bounded = d12(inside_bounds);
    if bounds(1)>-Inf
        d12_bounded=[bounds(1);d12_bounded];
    end
    if bounds(2)<Inf
        d12_bounded=[d12_bounded;bounds(2)];
    end
    d12 = d12_bounded;
end

if n1==n2 && ~bounded
    AM = sum(abs(x2-x1))/(2*n1);
elseif n1>n2 || bounded
    y1q = cdfinterpolator(x1,y1,d12); % Subroutine
    y2q = cdfinterpolator(x2,y2,d12); % Subroutine
    xdif = diff(d12);
    ydif = abs(y2q-y1q);
    ydif = ydif(2:end); 
    AM = sum(xdif.*ydif);
end

end


function yq = cdfinterpolator(x,y,xq) % This can be pretty slow on large datasets
n=length(xq);
yq = zeros(n,1);
for i=1:n
    pos = xq(i)>x;
    if any(pos)
        yqi=y(pos);
        yq(i)=yqi(end);
    elseif all(pos)
        yq(i)=1;
    end
end
end

function [xs,ps] = ecdf_Lpool(x)
%ECDF_LPOOL Summary of this function goes here
%%%
%   A quick script for finding the ecdf of data. 
%
%                   
%%%
n = length(x);
p = 1/n;
xs = sort(x(:));
ps = linspace(p,1,n)';
end