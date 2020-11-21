function [rho pval phi_0 amax] = circlin_wcorr(phi, x, bounds, weights)
%
% [rho pval ] = circlin_wcorr(phi, x)
%   Correlation coefficient between one circular and one linear random
%   variable, with weighting vector.
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
    bounds = [-4 4];
end
if nargin < 4 || isempty(weights)
   weights = ones(size(x));
else
    assert(isequal(size(weights),size(x)),'Vector of weights must be same dimensions as data')
end
    weights = weights./sum(weights);

R = @(a)-sqrt(mean(weights.*cos(phi - 2*pi*a*x))^2 + mean(weights.*sin(phi - 2*pi*a*x))^2);
%{
amax =  fminbnd(R,bounds(1),bounds(2));

S1 = nansum(weights.*sin(phi - 2*pi*amax*x));
C1 = nansum(weights.*cos(phi - 2*pi*amax*x));
phi_0 = arctan4quad(S1,C1);


% Convert linear variable, x, into circular variable, theta:
theta_fit = mod(2*pi*abs(amax)*x,2*pi);

%     % Eqn 4, Kempter et al 2012
    S2 = nansum(weights.*sin(phi));
    C2 = nansum(weights.*cos(phi));
    phi_mean = arctan4quad(S2,C2);
    S3 = nansum(weights.*sin(theta_fit));
    C3 = nansum(weights.*cos(theta_fit));
    theta_mean = arctan4quad(S3,C3);

    N = sum(weights.*sin(phi - phi_mean).*sin(theta_fit-theta_mean));
    D = sqrt( sum( weights.*(sin(phi - phi_mean)).^2) .* sum(weights.* (sin(theta_fit - theta_mean)).^2) );
    rho = N/D;

[rho, pval] = circ_wcorrcc(phi', theta_fit',weights');

if sign(amax)==-sign(rho)
    rho = abs(rho)*sign(amax);
end
%}

amax1 =  fminbnd(R,0,max(abs(bounds)));

S1 = nansum(weights.*sin(phi - 2*pi*amax1*x));
C1 = nansum(weights.*cos(phi - 2*pi*amax1*x));
phi_01 = arctan4quad(S1,C1);


% Convert linear variable, x, into circular variable, theta:
theta_fit = mod(2*pi*abs(amax1)*x,2*pi);

%     % Eqn 4, Kempter et al 2012
    S2 = nansum(weights.*sin(phi));
    C2 = nansum(weights.*cos(phi));
    phi_mean = arctan4quad(S2,C2);
    S3 = nansum(weights.*sin(theta_fit));
    C3 = nansum(weights.*cos(theta_fit));
    theta_mean = arctan4quad(S3,C3);

    % Original numerator from A.17
%     N = sum(weights.*sin(phi - phi_mean).*sin(theta_fit-theta_mean));
       
% Suggested modification for uniform data ( x ) : resultants of one
% parameter minus the other
N = sqrt (  mean(cos((phi -theta_fit))).^2 +  mean(sin((phi -theta_fit))).^2 )  -  sqrt (  mean(cos((phi +theta_fit))).^2 +  mean(sin((phi +theta_fit))).^2 ) ;

% % As above except : resultant of one, minus resultant of the other
%     N = (sqrt (  mean(cos((phi))).^2 +  mean(sin((phi))).^2 ) - sqrt (  mean(cos((theta_fit))).^2 +  mean(sin((theta_fit))).^2 )) - (sqrt (  mean(cos((phi))).^2 +  mean(sin((phi))).^2 ) + sqrt (  mean(cos((theta_fit))).^2 +  mean(sin((theta_fit))).^2 ));
% 
    D = sqrt( sum( weights.*(sin(phi - phi_mean)).^2) .* sum(weights.* (sin(theta_fit - theta_mean)).^2) );
    rhoI1 = N/D; % This is more robust when sampling of theta (x) is uniform (which it is). However it leads to correlation values greater than +/-1, so we only use it to compare pos vs neg fit

[rho1, pval1] = circ_wcorrcc(phi', theta_fit',weights');
 rss1 = sum( weights.*(sin(phi - phi_mean)).^2) * (1-rhoI1^2);

% r1 = abs(rhoI1)-abs(rho1)

amax2 =  fminbnd(R,-max(abs(bounds)),0);

S1 = nansum(weights.*sin(phi - 2*pi*amax2*x));
C1 = nansum(weights.*cos(phi - 2*pi*amax2*x));
phi_02 = arctan4quad(S1,C1);


% Convert linear variable, x, into circular variable, theta:
theta_fit = mod(2*pi*abs(amax2)*x,2*pi);

%     % Eqn 4, Kempter et al 2012
    S2 = nansum(weights.*sin(phi));
    C2 = nansum(weights.*cos(phi));
    phi_mean = arctan4quad(S2,C2);
    S3 = nansum(weights.*sin(theta_fit));
    C3 = nansum(weights.*cos(theta_fit));
    theta_mean = arctan4quad(S3,C3);
    
    % Original numerator from A.17

%     N = sum(weights.*sin(phi - phi_mean).*sin(theta_fit-theta_mean));

    % Suggested modification for uniform data ( x ) : resultants of one
% parameter minus the other
N = sqrt (  mean(cos((phi -theta_fit))).^2 +  mean(sin((phi -theta_fit))).^2 )  -  sqrt (  mean(cos((phi +theta_fit))).^2 +  mean(sin((phi +theta_fit))).^2 ) ;

% As above except : resultant of one, minus resultant of the other
%     N = (sqrt (  mean(cos((phi))).^2 +  mean(sin((phi))).^2 ) - sqrt (  mean(cos((theta_fit))).^2 +  mean(sin((theta_fit))).^2 )) - (sqrt (  mean(cos((phi))).^2 +  mean(sin((phi))).^2 ) + sqrt (  mean(cos((theta_fit))).^2 +  mean(sin((theta_fit))).^2 ));

    D = sqrt( sum( weights.*(sin(phi - phi_mean)).^2) .* sum(weights.* (sin(theta_fit - theta_mean)).^2) );
    rhoI2 = N/D; % This is more robust when sampling of theta (x) is uniform (which it is). However it leads to correlation values greater than +/-1, so we only use it to compare pos vs neg fit

[rho2, pval2] = circ_wcorrcc(phi', theta_fit',weights');
 rss2 = sum( weights.*(sin(phi - phi_mean)).^2) * (1-rhoI2^2);

% r2 = abs(rhoI2)-abs(rho2)
% 
% rI12 =  abs(rhoI2) -  abs(rhoI1)
% r12 = abs(rho2) -  abs(rho1)


if pval2<0.05 && pval1>=0.05
    choice = 2;
elseif pval2>=0.05 && pval1<0.05
    choice = 1;
elseif sign(rho1) == sign(rho2)  && sign(rhoI1) == sign(rhoI2) && sign(rhoI1) == sign(rho1)
    if sign(rho1) == -1
        choice = 2;
    else
         choice = 1;
    end
elseif pval2<0.05 && pval1<0.05 
    if abs(rhoI2)> abs(rhoI1) &&  abs(rho2)>abs(rho1)
        choice = 2;
    elseif abs(rhoI1)>abs(rhoI2)  &&  abs(rho1)>abs(rho2)
        choice = 1;
    elseif pval2< pval1
        choice = 2;
    elseif pval1<pval2
        choice = 1;
    else
        choice = 0;    
    end
else
    choice = 0;
end


if choice == 2
    rho = rho2;
    if abs(rho) >1
        rho = sign(rho);
    end
    pval = pval2;
    phi_0 = phi_02;
    amax = amax2;
elseif choice == 1
   rho = rho1;
    if abs(rho) >1
        rho = sign(rho);
    end
    pval = pval1;
    phi_0 = phi_01;
    amax = amax1;
elseif choice == 0
    
    amax =  fminbnd(R,-max(abs(bounds)),max(abs(bounds)));

S1 = nansum(weights.*sin(phi - 2*pi*amax*x));
C1 = nansum(weights.*cos(phi - 2*pi*amax*x));
phi_0 = arctan4quad(S1,C1);


% Convert linear variable, x, into circular variable, theta:
theta_fit = mod(2*pi*abs(amax)*x,2*pi);

%     % Eqn 4, Kempter et al 2012
    S2 = nansum(weights.*sin(phi));
    C2 = nansum(weights.*cos(phi));
    phi_mean = arctan4quad(S2,C2);
    S3 = nansum(weights.*sin(theta_fit));
    C3 = nansum(weights.*cos(theta_fit));
    theta_mean = arctan4quad(S3,C3);

    N = sum(weights.*sin(phi - phi_mean).*sin(theta_fit-theta_mean));
    
    D = sqrt( sum( weights.*(sin(phi - phi_mean)).^2) .* sum(weights.* (sin(theta_fit - theta_mean)).^2) );
    
    rho = N/D;

[rho, pval] = circ_wcorrcc(phi', theta_fit',weights');

if sign(amax)==-sign(rho)
    rho = abs(rho)*sign(amax);
end

end

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

function [rho pval] = circ_wcorrcc(alpha1, alpha2, W)
%
% [rho pval ts] = circ_corrcc(alpha1, alpha2)
%   Circular correlation coefficient for two circular random variables.
%
%   Input:
%     alpha1	sample of angles in radians
%     alpha2	sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%   Topics in circular statistics, S.R. Jammalamadaka et al., p. 176
%
% PHB 6/7/2008
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
W = W./sum(W);

if size(alpha1,2) > size(alpha1,1)
	alpha1 = alpha1';
end

if size(alpha2,2) > size(alpha2,1)
	alpha2 = alpha2';
end

if size(W,2) > size(W,1)
	W = W';
end

if length(alpha1)~=length(alpha2)
  error('Input dimensions do not match.')
end

% compute mean directions
n = length(alpha1);
alpha1_bar = circ_mean(alpha1,W);
alpha2_bar = circ_mean(alpha2,W);

% compute correlation coeffcient from p. 176
num = sum(W.*sin(alpha1 - alpha1_bar) .* sin(alpha2 - alpha2_bar))./sum(W);
den = sqrt( (sum(W.*sin(alpha1 - alpha1_bar).^2)./sum(W)) .* (sum(W.*sin(alpha2 - alpha2_bar).^2)./sum(W)) );
rho = num / den;	

% compute pvalue
l20 = mean(W.*sin(alpha1 - alpha1_bar).^2);
l02 = mean(W.*sin(alpha2 - alpha2_bar).^2);
l22 = mean(W.*(sin(alpha1 - alpha1_bar).^2) .* W.*(sin(alpha2 - alpha2_bar).^2));

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));

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
