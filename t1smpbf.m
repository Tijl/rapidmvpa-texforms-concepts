function bf10 = t1smpbf(t,n,r)
%
% bf10 = t1smpbf(t,n,[r=0.707])
%
% Calculates JZS Bayes Factor for a one-sample t-test given t and sample size n.
% The optional input r is the scale factor which defaults to 0.707.
% This quantifies the evidence in favour of the alternative hypothesis. 
% See Rouder et al, 2009, Psychon Bull Rev for details.
%
% original code obtained from https://doi.org/10.6084/m9.figshare.1357917.v1
% adapted to support a vector of t values
%

% Default scale factor
if nargin < 3
    r = 0.707;
end

% Function to be integrated
F = @(g,t,n,r) (1+n.*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+n.*g.*r.^2).*(n-1))).^(-n./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));

% Bayes factor calculation
% bf01 = (1 + t^2/(n-1))^(-n/2) / integral(@(g) F(g,t,n,r),0,Inf);
% Bayes factor calculation that handles arrays of t-values
bf01 = (1 + t.^2/(n-1)).^(-n/2) ./ integral(@(g) F(g,t,n,r),0,Inf,'ArrayValued',true);

% Invert Bayes Factor
% bf10 = 1 / bf01;
bf10 = 1 ./ bf01;
