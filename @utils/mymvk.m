function [out, R_, S] = mymvk(mu,R,N,a,b)
% OUT = MYMVK(MU,R,N,A,B)
% Generates observations from a Multivariate K distribution.
%
% MU = (1 x d)  :  Mean vector
% R  = (d x d)  :  Correlation matrix (if needed, R is tweaked internally
%                  into R_ in order to make it positive definite)
% N             :  Number of (1 x d) observations to be generated
% A  = ALPHA+1  :  ALPHA is the parameter defined in Taesu & Eltoft's paper
% B  = 1/LAMBDA :  LAMBDA is the parameter defined in Taesu & Eltoft's paper
%
% OUT = (N x d) :  The output sample, from Multivariate K distn.
% R_            :  The ground-truth correlation matrix (may differ from R).
% S             :  The ground-truth dispersion matrix.
%
% d             :  Length (dimensionality) of the distribution. This is
%                  infered from MU.
%
% The ouput covariance is: Cov(out) = ( (A*B)/(det(R_)^(1/d)) )*R_
% Alternatively, Cov(out) = (A*B)*S
%
% Test Code:
% This should display 2 matrices with all entries close to 1.
% R = [1 .2 .5 -.8; .2 1 .35 -.6; .5 .35 1 -.05; -.8 -.6 -.05 1];
% d = size(R,1);
% mu = zeros(1,d);
% N = 32968;
% [out, R_, S] = mymvk(mu,R,N,6,.1); corr(out)./R_, cov(out)./(.6*S)
%

d = length(mu);

if ~exist('a','var') || isempty(a)
    a = (d+1)/2;        % This reduces to Multivariate Laplace
end
if ~exist('b','var') || isempty(b)
    b = 2/(d+1);        % Makes Cov(out) = S when a = (d+1)/2
end

% Tweak R to make it positive definite:
R_ = R;%genvalidR(R); % If R is too bad, it will replace it with a random one...

% Make abs(det(S)) = 1
S = R_.*10^(-sum(log10(abs(eig(R_))))/d);% equivalent to mult. by 1/(abs(det(R_))^(1/d));

% N observations from a zero-mean, unit-std., uncorr., d-variate Gaussian distn.
xG = mvnrnd(zeros(1,d),eye(d),N);
% N observations from a gamma distribution with parameters a and b.
% To use Taesu & Eltoft notation, consider a = ALPHA+1 and b = 1/LAMBDA
zGm = gamrnd(a, b,[N 1]); % If a*b = 1, Cov(out) = S

% det(S) must be 1!!!!!!!!!!!!!!!!!!
out = repmat(mu',[N 1]) + repmat(sqrt(zGm),[1 d]).*(xG*(sqrtm(S)'));

end