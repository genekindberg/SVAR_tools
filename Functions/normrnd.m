function r = normrnd(mu,sigma,varargin)
%NORMRND Random arrays from the normal distribution.
%   R = NORMRND(MU,SIGMA) returns an array of random numbers chosen from a
%   normal distribution with mean MU and standard deviation SIGMA.  The size
%   of R is the common size of MU and SIGMA if both are arrays.  If either
%   parameter is a scalar, the size of R is the size of the other
%   parameter.
%
%   R = NORMRND(MU,SIGMA,M,N,...) or R = NORMRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMSTAT,
%   RANDOM, RANDN.

%   NORMRND uses Marsaglia's ziggurat method.

%   References:
%      [1] Marsaglia, G. and Tsang, W.W. (1984) "A fast, easily implemented
%          method for sampling from decreasing or symmetric unimodal density
%          functions", SIAM J. Sci. Statist. Computing, 5:349-359.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 2
    error(message('stats:normrnd:TooFewInputs'));
end

[err, sizeOut] = statsizechk(2,mu,sigma,varargin{:});
if err > 0
    error(message('stats:normrnd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

ty = internal.stats.dominantType(mu, sigma);
r = randn(sizeOut,'like',ty) .* sigma + mu;
