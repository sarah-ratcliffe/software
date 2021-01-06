function [out,stderr,resid] = rpat0(data,X,Z,W,patt,drop,g,maxit,tol,ind);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe  (copyright 2002)
% Last updated: 2 April 2002 
%
% COMMAND: [out,stderr] = rpat0(data,X,Z,W,patt,drop,g,maxit,tol,ind);
% ACTION:  Calculates longitudinal drop-out model with Gaussian data and
%          continuous drop-out with shared (random) pattern effect vector. 
%          Additionally, fixed effects and random effects (nested within
%          pattern) are available for the response data; and fixed effects
%          for the drop-out process. 
%          (Initial values must be in the Matlab database rpatinit.mat).
%
%          Assumes variance of u_2i = 0;
%
% INPUT:  data = [i Y];
%                i: subject number (must be unique for each subject, across all groups)
%                Y: response variables for each subject
%            X = fixed effects design matrix for each y(i). (stacked)
%            Z = subject random effects design matrix for each y(i). (stacked)
%	     W = pattern random effects design matrix for each y(i). (stacked)
%         patt = pattern associated with each subject. (n x 1)
%         drop = [delta Ti x2]; dropout information for all subjects
%                delta: vector of indicators for informative dropout
%                       0: noninform/complete/censored data
%                       1: informative dropout/ event occurred
%                Ti: dropout time or censoring time, depending on the 
%			value of delta  (n x 1)
%                x2:  fixed effects for t_i model  (n X q)
%            g = link function for R_i
%                1: linear link function
%                2: log link function [default]
%        maxit = maximum number of iterations to perform [default=100]
%          tol = convergence tolerance [default=1e-3]
%          ind = output indicator:
%                0: estimates at each iteration [default]
%                1: only iteration number is printed
%                2: same as 1 but also converged values at different tol values
%                   saved (for use with simulation studies).
%                3: no printed output
%
% OUTPUT:  out = parameter estimates.
%                out.alpha1 = parameter estimates for y(i) fixed effects.
%                out.SigBeta = covariance matrix for subject random effects 
%                out.Sigu1 = covariance matrix for y random pattern effects.
%                out.sigu2 = variance of R random pattern effect.
%                out.sig2e = variance of the errors for the y(i) model.
%                out.alpha2 = parameter estimates for R(i) fixed effects.
%                out.b = parameter estimate for R(i) pattern effect.
%                out.sig2R = variance of the errors for the R(i) model.
%                out.iterations = number of iterations performed.
%          stderr = standard errors of fixed effects (alpha's) estimates.
%          resid  = residuals from the fitted model.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Define initial values
if exist('rpatinit.mat')==2
  load rpatinit;
  disp('Initial values loaded');
else
  alpha1 = ones(size(X,2),1);
  [r,c] = size(drop);
  alpha2 = ones(c-2,1);
  SigBeta = eye(size(Z,2));
  Sigu1 = eye(size(W,2));
  sigu2 = 0;
  sig2e = 1;
  sig2R = 1;
  b = ones(size(W,2)+1,1);
end;

sigu2 = 0;

%% Error Checks
if nargin<10;  ind = 0;    end;		% Estimates printed at each iteration
if nargin<9;   tol = 1e-3; end;		% Convergence Tolerance
if nargin<8; maxit = 100;  end;		% Maximum number of iterations
if nargin<7;  g = 2; 
elseif (nargin==0) | (nargin>10)
  disp('ERROR: incorrect number of parameters');
  disp('       rpat(data=[i Y], X, Z, W, patt, drop=[delta Ti x2], g, maxit, tol)');
  return;
else
end;

[sum_ni,cd] = size(data);
if cd~=2  % Not 2 columns
  disp('ERROR: data must contain 2 columns: [i Y]');
  return;
end;

[ni,n] = counti(data(:,1));

[rx,p] = size(X);
if (rx~=sum_ni)
  disp('ERROR: X must have same number of rows at data');
  return;
end;
[rz,nz] = size(Z);
if (rz~=sum_ni)
  disp('ERROR: Z must have same number of rows at data');
  return;
end;
[rw,nw] = size(W);
if (rw~=sum_ni)
  disp('ERROR: W must have same number of rows at data');
  return;
end;
clear cd rx rz rw;

% Sort data so all subjects from same pattern are together
[patts,temp] = sortr(patt,[data X Z W]);
Y = temp(:,2);            xend = 2+p;
X = temp(:,3:xend);       zend = xend+nz;
Z = temp(:,xend+1:zend);  wend = zend+nw;
W = temp(:,zend+1:wend);
clear xend zend wend patts temp;
[patt,indpatt] = sort(patt);
drop = drop(indpatt,:);

% Get pattern counts and make sure patterns are 1,2,3,...
[mi,m] = counti(patt);
sum_mi = sum(mi);
patt = [];
for i=1:m
  patt = [patt; i*ones(mi(i),1)];
end;

% Get start indexes for each pattern / subject within the repeated data.
nstart = [1 cumsum(ni)+1];	% subject start indices
pstart = [1 cumsum(mi)+1];	% Pattern start indices
istart = nstart(pstart);	% Pattern groups at subject level indices
  
[rd,cd] = size(drop);
if rd~=n
  disp('ERROR: dropout informative not consistent with number of subjects');
  return;
end
if cd<3
  disp('ERROR: drop variable must contain at least 3 columns [delta Ti w]');
  return;
end;

delta = drop(:,1);
Ti = drop(:,2);
x2 = drop(:,3:cd);
q = cd - 2;
clear rd cd;

if g==1    % Linear link function - do nothing
elseif g==2  % Log link function
  Ti = log(Ti);
else
  disp('ERROR: invalid link function g');
  return;
end;


%% Set up various subject matrices
for i=1:m
 Xi = X(istart(i):istart(i+1)-1,:);
 Wi = W(istart(i):istart(i+1)-1,:);
 x2i = x2(pstart(i):pstart(i+1)-1,:);
 sum_X1tW(:,:,i) = Xi'*Wi;
 sum_x21(:,i) = x2i'*ones(mi(i),1);
end
sum_X1ty = X'*Y;
sum_x2t1 = sum(sum_x21,1);
inv_X1tX1 = inv(X'*X);
inv_x2x2t = inv(x2'*x2);
inv_ni    = inv(sum_ni);
inv_n     = inv(n);
  

%%%%%%%%%%%%%%%%%%%%% Parameter Estimation via EM algorithm %%%%%%%%%%%%%%%%%%
conver = 0;
if ind==2;
  conver3 = 0;
  conver4 = 0;
end;
iter = 0;
if ind==0;
  disp(['Iter  alpha1  SigBeta ' blanks(10*(nz-1)+4) 'Sigu1' blanks(10*(nw-1)+4) 'sig2e    alpha2    b    sigu2   sig2R']);
elseif ind==1 | ind==2 | ind==4;
  fprintf('\nIteration: ');
else
end

while conver==0		% Do until convergence or max iterations

%% Display stuff
if ind==0;
  disp('  ');
  fprintf(' %1.0f  %1.4f ',iter,alpha1(1,1));
  fprintf('  %0.4f ',SigBeta(1,:));
  fprintf('  %0.4f ',Sigu1(1,:));
  fprintf('  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f \n',sig2e,alpha2(1,1),b(1,1),sigu2,sig2R);
  for i=2:max([p nz nw size(alpha2,1)]);
    fprintf(' ');
    if i<=p;
      fprintf('   %0.4f ',alpha1(i,1));
    else
      fprintf('          ');
    end
    if i<=nz;
      fprintf('  %0.4f ',SigBeta(i,:));
    else
      fprintf(blanks(nz*length(num2str(SigBeta(1,1)))));
    end;
    if i<=nw;
      fprintf('  %0.4f ',Sigu1(i,:));
    else
      fprintf(blanks(nw*length(num2str(Sigu1(1,1)))+3));
    end;
    fprintf(blanks(length(num2str(sig2e))+7));
    if i<=size(alpha2,1);
      fprintf('  %0.4f',alpha2(i,1));
    else
      fprintf('        ');
    end
    if i<=nw;
      fprintf('  %0.4f \n',b(i,1));
    else
      fprintf('\n');
    end;
  end
  fprintf('\n');
elseif ind==1 | ind==2 | ind==4
  fprintf(' %1.0f',iter);
else
end;
iter = iter+1;

%% Reset sums to zero for next iteration
sum_al1  = sum_X1ty;
sum_al2  = 0;
sum_b1   = 0;
sum_b2   = 0;
sum_beta = 0;
sum_ete  = 0;
sum_sR2  = 0;
sum_su2  = 0;
u = zeros(nw+1,m);
ucu = 0;
Z_Beta = [];
W_u = [];
eresidR = [];
Rexp = [];

%%%%%%%%%%%%%%%%%%%%%%%%% E-step %%%%%%%%%%%%%%%%%%%%%%%%%

bSb = b'*Sigu1*b;
Sigma_u = [Sigu1 (Sigu1*b); (b'*Sigu1) (bSb+sigu2)];

for i=1:m			%% Group/pattern level stuff %%
 Yi = Y(istart(i):istart(i+1)-1);
 Xi = X(istart(i):istart(i+1)-1,:);
 Zi = Z(istart(i):istart(i+1)-1,:);
 Wi = W(istart(i):istart(i+1)-1,:);
 t  = Ti(pstart(i):pstart(i+1)-1);
 x2i = x2(pstart(i):pstart(i+1)-1,:)';
 faili = delta(pstart(i):pstart(i+1)-1);

 mu_yi = Xi * alpha1;
 mu_Ri = x2i'*alpha2;

 Sigma_Ri = (bSb + sigu2)*ones(mi(i)) + sig2R*eye(mi(i));
 Sigma_yR = Wi*Sigu1*b*ones(1,mi(i));

 invQi = [];
 DZi   = [];
 DXi   = [];
 sumij = 0;

 %% Calculations for invSigma_yi
 for j=1:mi(i)			% Subject within group
   sub = pstart(i) -1 + j;
   nij = ni(sub);
   Zij = Z(nstart(sub):nstart(sub+1)-1,:);
   invQij = inv(Zij*SigBeta*Zij' + sig2e*eye(nij));
   invQi  = [invQi zeros(sumij,nij); zeros(nij,sumij) invQij];
   if j==1
      DZi = [Zij];
   else
      DZi = [DZi zeros(sumij,nz); zeros(nij,nz*(j-1)) Zij];
   end;
   DXi    = [DXi; (1-faili(j))*ones(nij,p)];	% censored subjects
   sumij  = sumij + nij;
 end;
 temp = invQi*Wi;
 invSigma_yi = invQi - temp*inv(Wi'*temp + inv(Sigu1))*temp';

 F = invSigma_yi*Sigma_yR;
 invE = inv(Sigma_Ri - Sigma_yR'*F);
 temp = F*invE;
 invC22 = [(invSigma_yi+temp*F') (-temp); (-temp') (invE)];
 clear temp invE;


 %% Calculations for E(Ri | ...) and E(Ri Ri^T | ...)
 mu_Rgy   = mu_Ri + F'*mu_yi;
 sig2_Rgy = Sigma_Ri - Sigma_yR'*F;
 sig_Rgy  = ((sig2_Rgy)^(1/2));
 isig_Rgy = inv(sig_Rgy);
 temp     = isig_Rgy*(t-mu_Rgy);			% temp ~ N(0,I);
%temp = (1-faili).*temp;
 H        = normpdf(temp,0,1)./(1-normcdf(temp,0,1));	% H ~ N(0,I);
H = (1-faili).*H;

 Ri    = t.*faili + (1-faili).*(mu_Rgy + sig_Rgy*H);
 Faili = diag(1-faili);
% RiRi  = Ri*Ri' + ((1-faili)*(1-faili)').*(sig_Rgy*diag((1 + temp.*H - H.^2))*sig_Rgy');
 RiRi  = Ri*Ri' + ((1-faili)*(1-faili)').*(sig_Rgy*diag((1-faili).*(1 + temp.*H - H.^2))*sig_Rgy');
 Rexp = [Rexp; Ri];

 clear mu_Rgy sig2_Rgy sig_Rgy isig_Rgy temp H Faili F;


 %% Calculations for E(u_i | ..) & E(u_i u_i' | ...)
 tempy = Yi - mu_yi;
 tempR = Ri - mu_Ri;
 temp1 = Sigu1*Wi';
 C12 = [temp1 (Sigu1*b*ones(1,mi(i))); (b'*temp1) ((bSb+sigu2)*ones(1,mi(i)))];
 temp3 = C12*invC22;
 ui = temp3*[tempy; tempR];		% E(u_i | ...)
 u(:,i) = ui;
 cov_ui  = Sigma_u - temp3*C12';
 Ai2 = RiRi - Ri*mu_Ri' - mu_Ri*Ri' + mu_Ri*mu_Ri';
 Ai = [(tempy*tempy') (tempy*tempR'); (tempR*tempy') Ai2];
 uiui = temp3*Ai*temp3' + cov_ui;	% E(u_i u_i' | ...)

 sum_al1 = sum_al1 - sum_X1tW(:,:,i)*ui(1:nw,1);   %% sums used in M-step %%
 sum_al2 = sum_al2 + x2i*Ri - sum_x21(:,i)*ui(nw+1,1);
% sum_b1  = sum_b1 + mi(i)*uiui(nw+1,1:nw);
 sum_b1  = sum_b1 + sum(tempR)*ui(1:nw,1)';
 sum_b2  = sum_b2 + mi(i)*uiui(1:nw,1:nw);
 W_u = [W_u; Wi*ui(1:nw)];

 clear temp3 cov_ui Ai2 C12 ui uiui;


 %% Calculations for E(u_1' b b' u_1 + u_2^2 | ...) & E(u_2 b' u_1 | ...)
 C12 = [(b'*temp1) (bSb*ones(1,mi(i))); (b'*temp1) ((bSb+sigu2)*ones(1,mi(i)))];
 temp3 = C12*invC22;
 Ebui = temp3*[tempy; tempR];		% E([b'u_1; u_2] | ...)
 Sigma_bu = [bSb bSb; bSb (bSb+sigu2)];
 cov_bui = Sigma_bu - temp3*C12';
 tempi  = temp3'*temp3;
 D11    = tempi(1:sumij,1:sumij);
 D12    = tempi(1:sumij,sumij+1:sumij+mi(i));
 D22    = tempi(sumij+1:sumij+mi(i),sumij+1:sumij+mi(i));
 Esu2_1 = tempy'*D11*tempy + tempy'*D12*tempR + tempR'*D12'*tempy + trace(D22*RiRi) - Ri'*D22*mu_Ri - mu_Ri'*D22*(Ri - mu_Ri) + trace(cov_bui);
 Esu2_2 = temp3*Ai*temp3' + cov_bui;

 sum_su2 = sum_su2 + Esu2_1 + Esu2_2(1,2);	%% sum used in M-step %%

 clear C12 temp3 Ebui Sigma_bu cov_bui tempi D11 D12 D22 Esu2_1 Esu2_2;


 %% Calculations for E(beta_i|y_i,R_i,...) & E(beta_i beta_i^T |y_i,R_i,...)
 D11 = kron(SigBeta,eye(mi(i)));
 D12 = [(D11*DZi') zeros(nz*mi(i),mi(i))];
 temp3 = D12*invC22;
 Ebeta = temp3*[tempy; tempR];
 cov_b = D11 - temp3*D12';
 Ebbt  = temp3*Ai*temp3' + cov_b;

 for j=1:mi(i);					%% sums used in M-step
   indj = nz*(j-1);
   sum_beta = sum_beta + Ebbt(indj+1:indj+nz,indj+1:indj+nz);
   sub = pstart(i) -1 + j;
   nij = ni(sub);
   Zij = Z(nstart(sub):nstart(sub+1)-1,:);
   Z_Beta = [Z_Beta; Zij*Ebeta(j,:)'];
 end;
 sum_al1 = sum_al1 - Xi'*DZi*Ebeta;

 clear D11 D12 temp3 Ebeta cov_b Ebbt indj;


 %% Calculations for E(e_i' e_i | ...)	- y error term
 D11 = sig2e*eye(sumij);
 D12 = [D11 zeros(sumij,mi(i))];
 Eerr = D12*invC22*[tempy; (Ri-mu_Ri)];
 cov_e = D11 - D12*invC22*D12';
 F11   = invC22(1:sumij,1:sumij);
 F12   = invC22(1:sumij,sumij+1:sumij+mi(i));
 D22   = F12'*F12;
 Eete  = sig2e^2*(tempy'*F11*F11*tempy + tempy'*F11*F12*tempR + tempR'*F12'*F11*tempy + trace(D22*RiRi) - (Ri'*D22*mu_Ri) + (mu_Ri'*D22*(mu_Ri - Ri))) + trace(cov_e);

 sum_ete = sum_ete + Eete;			%% sum used in M-step
 
 clear D11 D12 D22 cov_e F11 Eete;


 %% Calculations for E(err_i' err_i | ...)  - R error term
 D11 = sig2R*eye(mi(i));
 D12 = [zeros(mi(i),sumij) D11];
 Eerr = D12*invC22*[tempy; (Ri-mu_Ri)];
 eresidR = [eresidR; Eerr];
 cov_e = D11 - D12*invC22*D12';
 F22 = invC22(sumij+1:sumij+mi(i),sumij+1:sumij+mi(i));
 D22 = F22*F22;
 Eete  = sig2R^2*(tempy'*F12*F12'*tempy + tempy'*F12*F22*tempR + tempR'*F22*F12'*tempy + trace(D22*RiRi) - (Ri'*D22*mu_Ri) + (mu_Ri'*D22*(mu_Ri - Ri))) + trace(cov_e);

 sum_sR2 = sum_sR2 + Eete;			%% sum used in M-step

 clear D11 D12 D22 cov_e F12 F22 Eete;
 
 
 %% Standard Errors for fixed effects 
 % all subjects
 Ui = [Xi zeros(sumij,q); zeros(mi(i),p) x2i'];
 ucu = ucu + Ui'*invC22*Ui;
 % censored subjects only
 Ui = [(Xi.*DXi) zeros(sumij,q); zeros(mi(i),p) (x2i'.*(faili*ones(1,q)))]; 
 ucu = ucu + Ui'*invC22*Ai*invC22*Ui;

end;				% End of group calculations


%%%%%%%%%%%%%%%%%%%%%%%%% M-step %%%%%%%%%%%%%%%%%%%%%%%%%
al1_old = alpha1;
al2_old = alpha2;
b_old = b;
SB_old = SigBeta;
Su1_old = Sigu1;
sigu2_old = sigu2;
su2_old = sigu2;
se_old = sig2e;
sR_old = sig2R;

alpha1 = inv_X1tX1 * sum_al1;
alpha2 = inv_x2x2t * sum_al2;
b = sum_b1 * inv(sum_b2);
b = b';
SigBeta = inv_n * sum_beta;
SigBeta = (SigBeta + SigBeta')/2;
Sigu1 = inv_n * sum_b2;
Sigu1 = (Sigu1 + Sigu1')/2;
%sigu2 = inv_n * sum_su2;
sig2e = inv_ni * sum_ete;
sig2R = inv_n * sum_sR2;


%%% Check for Convergence %%%
c_sum = 0;
if relerr(alpha1,al1_old,tol);  c_sum = c_sum+1; end;
if relerr(alpha2,al2_old,tol);  c_sum = c_sum+1; end;
if relerr(b,b_old,tol);         c_sum = c_sum+1; end;
if relerr(SigBeta,SB_old,tol);  c_sum = c_sum+1; end;
if relerr(Sigu1,Su1_old,tol);	c_sum = c_sum+1; end;
if relerr(sigu2,su2_old,tol);   c_sum = c_sum+1; end;
if relerr(sig2e,se_old,tol);    c_sum = c_sum+1; end;
if relerr(sig2R,sR_old,tol);    c_sum = c_sum+1; end;

if c_sum==8;
  disp(['Convergence achieved after ' num2str(iter) ' iterations']);
  conver = 1;
end;
if iter==maxit;
  disp('Maximum number of iterations reached');
  conver = 1;
end;
if sum(isnan([alpha1; alpha2; b; sigu2; sig2e; sig2R; sum(SigBeta,2); sum(Sigu1,2)]))>0 ;
  disp(['One or more estimates = NaN']);
  conver = 1;
end; 

if rem(iter,10)==0
  save rpat0working alpha1 alpha2 b SigBeta Sigu1 sigu2 u sig2e sig2R iter;
end
end		% End of while loop


fprintf('\n');
out.alpha1 = alpha1;
out.alpha2 = alpha2;
out.b = b;
out.SigBeta = SigBeta;
out.Sigu1 = Sigu1;
out.sigu2 = sigu2;
out.Sigma_u  = Sigma_u;
out.u = u;
out.sig2e = sig2e;
out.sig2R = sig2R;
out.iterations = iter;

%%% Residuals %%
resid.Y_obs = Y;
resid.Y_hat = X*alpha1 + W_u + Z_Beta;
resid.Y_err = Y - resid.Y_hat;
resid.R_obs = Ti;
resid.R_exp = Rexp;
resid.drop = delta;
resid.R_hat = x2*alpha2 + (b'*u(nw+1,patt))';
resid.R_obs_err = Ti - resid.R_hat;
resid.R_exp_err = eresidR;

%%% Standard Errors %%%
inv_Hess = inv(ucu);
stdall = sqrt(diag(inv(ucu)));
stderr.alpha1 = stdall(1:length(alpha1));
stderr.alpha2 = stdall(length(alpha1)+1:length(stdall));

if (ind==4) | (ind==0);
  disp('Final Estimates');
  disp(['Iter  alpha1  SigBeta ' blanks(10*(nz-1)+4) 'Sigu1' blanks(10*(nw-1)+4) 'sig2e    alpha2     b    sigu2   sig2R']);
  disp('  ');
  fprintf(' %1.0f  %0.4f ',iter,alpha1(1,1));
  fprintf('  %0.4f ',SigBeta(1,:));
  fprintf('  %0.4f ',Sigu1(1,:));
  fprintf('  %0.4f  %0.4f  %0.4f  %0.4f  %0.4f \n',sig2e,alpha2(1,1),b(1,1),sigu2,sig2R);
  for i=2:max([p nz nw size(alpha2,1)]);
    if i<=p;
      fprintf('   %0.4f ',alpha1(i,1));
    else
      fprintf('          ');
    end
    if i<=nz;
      fprintf('  %0.4f ',SigBeta(i,:));
    else
      fprintf(blanks(nz*length(num2str(SigBeta(1,1)))+3));
    end;
    if i<=nw;
      fprintf('  %0.4f ',Sigu1(i,:));
    else
      fprintf(blanks(nw*length(num2str(Sigu1(1,1)))+3));
    end;
    fprintf(blanks(length(num2str(sig2e))+2));
    if i<=size(alpha2,1);
      fprintf('  %0.4f',alpha2(i,1));
    else
      fprintf('        ');
    end
    if i<=nw;
      fprintf('  %0.4f \n',b(i,1));
    else
      fprintf('\n');
    end;
  end
  fprintf('\n');
end;
