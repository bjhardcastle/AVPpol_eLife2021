function [rho pval phi_0 amax] = circlin_corr(phi, x, bounds)
%
% [rho pval ] = circlin_corr(phi, x)
%   Correlation coefficient between one circular and one linear random
%   variable.
%
%   Input:
%     phi     sample of angles in radians
%     x       sample of linear random variable
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%     phi_0   phase offset
%     amax    slope
%
% References:
%     Kempter et al, 2012 - Quantifying circular–linear associations: Hippocampal phase precession
%
% BJH 26/6/2019
% modified from circ_corrcl in circstats toolbox for Matlab, by P Berens


if size(phi,2) > size(phi,1)
    phi = phi';
end

if size(x,2) > size(x,1)
    x = x';
end

if length(phi)~=length(x)
    error('Input dimensions do not match.')
end

n = length(phi);

if nargin < 3 || isempty(bounds)
    bounds = [-1 1];
end


R = @(a)-sqrt(mean(cos(phi - 2*pi*a*x))^2 + mean(sin(phi - 2*pi*a*x))^2);
amax =  fminbnd(R,bounds(1),bounds(2));

S1 = sum(sin(phi - 2*pi*amax*x));
C1 = sum(cos(phi - 2*pi*amax*x));
phi_0 = arctan4quad(S1,C1);


% Convert linear variable, x, into circular variable, theta:
theta_fit = mod(2*pi*abs(amax)*x,2*pi);

%     % Eqn 4, Kempter et al 2012
%     S2 = sum(sin(phi));
%     C2 = sum(cos(phi));
%     phi_mean = arctan4quad(S2,C2);
%     S3 = sum(sin(theta_fit));
%     C3 = sum(cos(theta_fit));
%     theta_mean = arctan4quad(S3,C3);
% 
%     N = sum(sin(phi - phi_mean).*sin(theta_fit-theta_mean));
%     D = sqrt( sum( (sin(phi - phi_mean)).^2) .* sum( (sin(theta_fit - theta_mean)).^2) );
%     rho = N/D;

[rho, pval] = circ_corrcc(phi', theta_fit');
% 
% lincorr = corr(phi,x);
% 
% if sign(lincorr)==-1 && sign(rho)==1
% rho = abs(rho)*sign(lincorr);
% amax = abs(amax)*sign(lincorr);
% end

end

function value = arctan4quad(O,A)
% calculate quadrant-specific inverse tangeant from two lengths, opposite and adjacent
% eqn A.6 Kempter et al 2012
if A>0 && O>=0
    value = atan(O/A);
elseif A==0 && O>0
    value = pi/2;
elseif A<0
    value = atan(O/A) + pi;
elseif A>=0 && O<0
    value = atan(O/A) + 2*pi;
elseif A==0 && O==0
    value = 0;
end
end


% eq3


%
% % compute correlation coefficent for sin and cos independently
% rxs = corr(x,sin(alpha));
% rxc = corr(x,cos(alpha));
% rcs = corr(sin(alpha),cos(alpha));
%
% % compute angular-linear correlation (equ. 27.47)
% rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));
% % rho = sign(rxs)*sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));
%
% % compute pvalue
% pval = 1 - chi2cdf(n*rho^2,2);
%
