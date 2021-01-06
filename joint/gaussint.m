function SS = gaussint(funct,alim,blim,mu,variance,b);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Function used to calculate integrals via
% Hermite Integration. (from -\infty to \infty)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin<6; b = 1; end;
if nargin<5; variance = 1; end;
if nargin<4; mu = 0; end;
if nargin<2; alim = -1; end;
if nargin<3; blim = -alim; end;

npts = 48;
load gauss48pts.dat;
X = gauss48pts(:,1);
W = gauss48pts(:,2);
clear gauss48pts;

SS=0;
D1 = (blim-alim)/2;
D2 = (blim+alim)/2;

% CALCULATE EXPECTED VALUES
for j=1:(npts/2)
 x1= D1*X(j)+D2;
 x2= -D1*X(j)+D2;

 switch funct
  case 1		% Normal function
    FF1 = inv(sqrt(2*pi*variance))*exp(-((x1-mu)^2)/(2*variance));
    FF2 = inv(sqrt(2*pi*variance))*exp(-((x2-mu)^2)/(2*variance));

  case 2		% E(x) from normal function
    gx1 = x1;
    FF1 = gx1*inv(sqrt(2*pi*variance))*exp(-((x1-mu)^2)/(2*variance));
    gx2 = x2;
    FF2 = gx2*inv(sqrt(2*pi*variance))*exp(-((x2-mu)^2)/(2*variance));

  case 3  		% Calculate E(exp(b*x)), where x is normally dist.
    gx1 = exp(b*x1);
    FF1 = gx1*inv(sqrt(2*pi*variance))*exp(-((x1-mu)^2)/(2*variance));
    gx2 = exp(b*x2);
    FF2 = gx2*inv(sqrt(2*pi*variance))*exp(-((x2-mu)^2)/(2*variance));

  case 4		% 
    FF1 = exp(b*x1);
    FF2 = exp(b*x2);

  case 5		%
    FF1 = exp(b*x1)/(1+exp(b*x1));
    FF2 = exp(b*x2)/(1+exp(b*x2));

  case 6		% Calculate E(logit(b*x)) where x is normally dist
    gx1 = exp(b*x1)/(1+exp(b*x1));
    FF1 = gx1*inv(sqrt(2*pi*variance))*exp(-((x1-mu)^2)/(2*variance));
    gx2 = exp(b*x2)/(1+exp(b*x2));
    FF2 = gx2*inv(sqrt(2*pi*variance))*exp(-((x2-mu)^2)/(2*variance));

  case 7		% E(x^2) from normal function
    gx1 = x1^2;
    FF1 = gx1*inv(sqrt(2*pi*variance))*exp(-((x1-mu)^2)/(2*variance));
    gx2 = x2^2;
    FF2 = gx2*inv(sqrt(2*pi*variance))*exp(-((x2-mu)^2)/(2*variance));

  otherwise	% Not a valid option
    disp('ERROR: invalid integral function to evaluate');
    return;
 end;

 SS=SS+W(j)*(FF1+FF2);
end;

SS = SS*D1;

switch funct
  case 1
    trueval = 1;
  case 2
    trueval = mu;
  case 3
    trueval = exp(b*mu + variance*b^2/2);
  case 4
    trueval = exp(b*mu + variance*b^2/2);
  case 5
    trueval = NaN;
  case 6
    trueval = NaN;
  case 7
    trueval = mu^2 + variance;
  otherwise
    trueval = NaN;
end;
    
disp(['True Value = ' num2str(trueval)]);
disp(['Est. Value = ' num2str(SS)]);
