function [out, R_, S] = mymvlap(O,mu,R,N,b)
% OUT = MYMVLAP(MU,R,N,B)
% Generates observations from a Multivariate Laplace distribution.
%
% MU = (1 x d)  :  Mean vector
% R  = (d x d)  :  Correlation matrix (if needed, R is tweaked internally
%                  into R_ in order to make it positive definite)
% N             :  Number of (1 x d) observations to be generated
% B  = 1/LAMBDA :  LAMBDA is the parameter defined in Taesu & Eltoft's paper
%
% d             :  Length (dimensionality) of the distribution. This is
%                  infered from MU.
%
% OUT = (N x d) :  The output sample, from Multivariate Laplace distn.
% R_            :  The ground-truth correlation matrix (may differ from R).
% S             :  The ground-truth dispersion matrix.
%
% The ouput covariance is: Cov(out) = ( B*((d+1)/2)/(det(R_)^(1/d)) )*R_
% Alternatively, Cov(out) = B*((d+1)/2)*S
%
% Test Code:
% This should display 2 matrices with all entries close to 1.
% R = [1 .2 .5 -.8; .2 1 .35 -.6; .5 .35 1 -.05; -.8 -.6 -.05 1];
% d = size(R,1);
% mu = zeros(1,d);
% N = 32968;
% [out, R_, S] = mymvlap(mu,R,N,2);
% corr(out)./R_, cov(out)./(2*(d+1)*S./2)
%

% ut = utils;

d = length(mu);
a = (d+1)/2; % This reduces the Multivariate K distn. to Multivariate Laplace

if ~exist('b','var') || isempty(b)
    b = 2/(d+1);    % Makes Cov(out) = S when a = (d+1)/2
end

if      nargout == 1,   out = O.mymvk(mu,R,N,a,b);
elseif  nargout == 2,   [out, R_] = O.mymvk(mu,R,N,a,b);
else                    [out, R_, S] = O.mymvk(mu,R,N,a,b);
end

end