function [out,stderr,std2,resid] = surv(data,X,Z,patt,drop,maxit,tol,ind);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe  (copyright 2001)
%
% COMMAND: [out,stderr] = surv(data,X,Z,patt,drop,maxit,tol,ind);
%
% ACTION: ECM algorithm using expected newton-raphson algorithm.
%
%	   (Initial values must be in the Matlab database survinit.mat).
% INPUT:  data = [i Y];
%                i: subject number (must be unique for each subject, across all groups)
%                Y: response variables for each subject
%            X = fixed effects corresponding to each y(i). (stacked)
%            Z = random effects corresponding to each y(i). (stacked)
%         patt = pattern associated with each subject. (n x 1)
%         drop = [delta Ti x2]; dropout information for all subjects
%                delta: vector of indicators for informative dropout
%                       0: noninform/complete data
%                       1: informative dropout
%                Ti: dropout time (R_i) or censoring time (c_i), 
%                    depending on value of delta  (n x 1)
%                x2:  fixed effects for t_i model  (n X q)
%        maxit = maximum number of iterations to perform [default=100]
%          tol = convergence tolerance [default=1e-5]
%          ind = output indicator:
%                0: estimates at each iteration [default]
%                1: only iteration number is printed
%                2: same as 1 but also converged values at different tol values
%                   saved (for use with simulation studies).
%                3: no printed output
%
% OUTPUT:  out.* contains the parameter estimates and stderr.* the standard 
%          errors for the parameters listed below.
%             alpha1 = y(i) fixed effects parameters.
%            SigBeta = covariance matrix of subject random effects 
%              sig2u = variance of random pattern effects.
%              sig2e = variance of the errors for the y(i) model.
%             alpha2 = R(i) fixed effects parameters.
%                  b = R(i) pattern effects parameters.
%               iter = number of iterations performed.
%          resid = residuals from the fitted model, adjusted for censoring.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stdcal = 1;
u = zeros(max(patt),1);
const = [1e-13; 1e-25; 1e-31; 1e-20; 1e-31];
%% Define initial values
if exist('survinit.mat')==2
  load survinit;
else
  alpha1 = ones(size(X,2),1);
  [r,c] = size(drop);
  alpha2 = ones(c-2,1);
  SigBeta = eye(size(Z,2));
  sig2u = 1;
  sig2e = 1;
  sig2R = 1;
  b = 1;
end;
dimSB = size(SigBeta,1);
vecSB = dimSB^2;

%% Error Checks
if nargin<8;   ind = 0;    end;		% Estimates printed at each iteration
if nargin<7;   tol = 1e-5; end;		% Convergence Tolerance
if nargin<6; maxit = 100;  end;		% Maximum number of iterations
if (nargin<5) | (nargin>9)
  disp('ERROR: incorrect number of parameters');
  disp('       surv(data=[i Y], X, Z, patt, drop=[delta Ti w], maxit, tol, ind)');
  return;
end;


[sum_ni,cd] = size(data);
if cd~=2  % Not 2 columns
  disp('ERROR: data must contain 2 columns: [i Y]');
  return;
end;

[ni,n] = counti(data(:,1));
Y = data(:,2);

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
clear cd rx rz;

% Sort data so all subjects from same pattern are together
[patts,temp] = sortr(patt,[data X Z]);
Y = temp(:,2);          xend = 2+p;
X = temp(:,3:xend);       zend = xend+nz;
Z = temp(:,xend+1:zend);  
clear xend zend patts temp;
[patt,indpatt] = sort(patt);
drop = drop(indpatt,:);
clear indpatt;

% Get pattern counts and make sure patterns are 1,2,3,...
[mi,m] = counti(patt);
sum_mi = sum(mi);
patt = [];
for i=1:m
  patt = [patt; i*ones(mi(i),1)];
end;

[rd,cd] = size(drop);
if rd~=n
  disp('ERROR: dropout information not consistent with number of subjects');
  return;
end
if cd<3
  disp('ERROR: drop variable must contain at least 3 columns [delta Ti w]');
  return;
end;

delta = drop(:,1);
Ti = drop(:,2);
Tall = sort(unique(Ti));
x2 = drop(:,3:cd);
q = cd - 2;
clear rd cd;
nb = length(b);


%% Set up various subject matrices
sum_d = sum(delta);
nstart = 0;
sum_x2 = 0;
sum_x2x2t = 0;
sum_x21 = zeros(size(alpha2,1),max(patt));
sum_x2i = zeros(size(alpha2,1),max(patt));
nt = length(Tall);
indext = sum((Ti*ones(1,nt)) >= (ones(n,1)*Tall'),2);

for i=1:n
 sum_x2 = sum_x2 + delta(i)*x2(i,:);
 sum_x2x2t = sum_x2x2t + x2(i,:)'*x2(i,:);
 sum_x2i(:,patt(i)) = sum_x2i(:,patt(i)) + x2(i,:)';
 sum_x21(:,patt(i)) = sum_x21(:,patt(i)) + delta(i)*x2(i,:)';
end
inv_x2x2t = inv(sum_x2x2t);
inv_ni    = inv(sum_ni);
inv_n     = inv(n);
clear sum_x2x2t;

% Initial cummulative Hazard function
for i=1:nt;
  lam0(i) = sum((Ti==Tall(i)).*delta)/sum(Ti>=Tall(i));
end;

Lam0 = cumsum(lam0);

% Plot of initial Survivor Function
if ind==0
  figure(1);
  plot([0; Tall],[1 exp(-Lam0)],'k');
  hold on;
end;
  
%%%%%%%%%%%%%%%%%%%%%           ALGORITHM START            %%%%%%%%%%%%%%%%%%
conver = 0;
if ind==2;
  conver3 = 0;
  conver4 = 0;
end;
iter = 0;
if ind==0;
  disp('Iter  alpha1    SigBeta       sig2u   sig2e   alpha2    b');
  disp('  ');
  fprintf(' %1.0f',iter);
  fprintf('  %0.4f ',alpha1(1,1));
  fprintf(' %0.4f ',SigBeta(1,:));
  fprintf(' %0.4f  %0.4f  %0.4f  %0.4f \n ',sig2u,sig2e,alpha2(1,1),b);
  for i=2:max([size(alpha1,1) size(SigBeta,1) size(alpha2,1)]);
    if i<=size(alpha1,1);
      fprintf('    %0.4f ',alpha1(i,1));
    else
      fprintf('           ');
    end
    if i<=size(SigBeta,1);
      fprintf(' %0.4f ',SigBeta(i,:));
    else
      fprintf('                 ');
    end;
    if i<=size(alpha2,1);
      fprintf('                   %0.4f \n ',alpha2(i,1));
    else
      fprintf('\n ');
    end
  end
  fprintf('\n');
elseif ind==1 | ind==2
  fprintf('\nIteration: ');
  fprintf(' %1.0f',iter);
else
end;

while conver==0		% Do until convergence or max iterations
  iter = iter + 1;
  if ind<=2;
    fprintf(' %1.0f',iter);
  end;

%%%%%%%%%%%%%%%%%%%%% STEP 1:                               %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Parameter Estimation via EM algorithm %%%%%%%%%%%%%%%%%%
% Reset sums to zero for next iteration
%sum_plot = zeros(2,41);
sum_al1  = 0;
sum_XtX  = 0;
sum_al2  = 0;
sum_b1   = 0;
sum_b2   = 0;
sum_S2   = 0;
sum_lam  = 0;
ucu      = 0;
sum_ui   = 0;
sum_e    = 0;
sum_U2   = 0;
sum_vv   = 0;
sum_H1   = 0;
sum_H2   = 0; 
nstart   = 0;
subject  = 0;
subject2 = 0;
beta = [];
Eeu = [];
Eu = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E-STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:m			%% Group/pattern level stuff %%
 sum_X1ti = 0;
 sum_lami = 0;
 sum_temp1= 0;
 sum_temp2= 0;
 sum_temp3= 0;
 prod_ft = 1;
 sum_ft  = 0;
 prod_fy = 1;
 sum_fy1 = 0;
 sum_fy2 = 0;
 sum_fy3 = 0;
 sum_fy4 = 0;
 nid     = 0;
 sum_vv_1 = 0;
 sum_vv_2 = 0;
 sum_e_1  = 0;
 sum_e_2  = 0;
 sum_H2_1 = 0;
 sum_H2_2 = 0;

 groupstart = subject+1;

 for j=1:mi(i)			% Subject within group
   subject = subject + 1;
   nij = ni(subject);
   yij = Y(nstart+1:nstart+nij,1);
   Xij = X(nstart+1:nstart+nij,:);
   Zij = Z(nstart+1:nstart+nij,:);
   failij = delta(subject);
   Tij = Ti(subject);
   x2ij = x2(subject,:)';
   nstart = nstart + nij;
  
   mu_yij = Xij*alpha1;
   ex2ij = exp(x2ij'*alpha2);

   Qij = Zij*SigBeta*Zij' + sig2e*eye(nij);
   invQij = inv(Qij);

   % Standard Errors for fixed effects 
   SigXij = invQij*Xij;

   % Sums used in M-step  
   sum_XtX  = sum_XtX + Xij'*SigXij;
   sum_al1  = sum_al1 + SigXij'*yij;
   sum_X1ti = sum_X1ti + SigXij'*ones(nij,1);
   sum_lami = sum_lami + ((Tij >= Tall)*ex2ij); 

   % Sums used in Newton-Raphson for sig2e and SigBeta.
   inv_Vij = sig2e*invQij;
   temp1   = yij - mu_yij;
   sum_e   = sum_e + temp1'*inv_Vij*temp1;
   temp3   = inv_Vij*ones(nij,1)*temp1';
   sum_e_1 = sum_e_1 + trace(temp3) + trace(temp3');
   sum_e_2 = sum_e_2 + ones(1,nij)*inv_Vij*ones(nij,1);
   temp2  = Zij'*inv_Vij;
   vv1    = temp2*temp1*temp1'*temp2';
   sum_vv = sum_vv + vv1;
   vv2    = temp2*ones(nij,1)*temp1'*temp2';
   sum_vv_1 = sum_vv_1 + vv2 + vv2';
   vv3    = temp2*ones(nij,nij)*temp2';
   sum_vv_2 = sum_vv_2 + vv3;
   temp2  = temp2*Zij;
   temp   = kron(temp2,vv1);
   sum_H2 = sum_H2 + temp + temp';
   temp   = kron(temp2,vv2+vv2');
   sum_H2_1 = sum_H2_1 + temp + temp';
   temp   = kron(temp2,vv3);
   sum_H2_2 = sum_H2_2 + temp + temp';
   sum_U2 = sum_U2 + reshape(temp2,vecSB,1); 
   sum_H1 = sum_H1 + kron(temp2,temp2);

   % Products needed for expected values.
   if failij == 1
     nid = nid + 1;
     prod_ft = prod_ft*Lam0(indext(subject))*ex2ij;
   end
   sum_ft  = sum_ft + Lam0(indext(subject))*ex2ij;
   prod_fy = prod_fy*(2*pi)^(-nij/2)*(det(Qij)^(-1/2));
   sum_fy1 = sum_fy1 + temp1'*invQij*temp1;
   sum_fy2 = sum_fy2 + ones(1,nij)*invQij*temp1;
   sum_fy3 = sum_fy3 + temp1'*invQij*ones(nij,1);
   sum_fy4 = sum_fy4 + ones(1,nij)*invQij*ones(nij,1);   

 end				% End of subject loop for group calculations
 
 % Expected values for functions of u_i
 % Integral in f(u) evaluated at 5 points via quassian quadrature
 ui = u(patt(subject));
 limits = ui + [-2 2]*sqrt(sig2u); 
 const_fu = (2*pi*sig2u)^(-1/2);
 const_ft = exp(b*ui + b^2*sig2u/2);
 exp_bui = b*nid*ui + (b*nid)^2*sig2u/2;
 prod_fy = prod_fy*inv(-sum_fy1/2 + 0.5*sum_fy2 -0.125*sum_fy4);

 fui = survfui(1,limits,const_ft,0,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);
%fprintf(' fui = %0.16f \n',fui);
 ui = survfui(2,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(u_i | ...)
 uiui = survfui(3,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(u_i^2 | ...)
 Eeui = survfui(4,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(exp{b u_i}|...)
 Eeu = [Eeu; ones(mi(i),1)*Eeui];
 Eu = [Eu; ones(mi(i),1)*ui];

 % Sum used in M-step
 sum_al1 = sum_al1 - sum_X1ti*ui;
 sum_b1 = sum_b1 + mi(i)*uiui;
 sum_lam = sum_lam + sum_lami*Eeui;
 u(patt(subject)) = ui;

 % Sums used in calculations for gamma
 sum_S2 = sum_S2 + nid*ui;

 % Sums used in Newton-Raphson for Variance components.
 sum_vv = sum_vv - ui*sum_vv_1 + uiui*sum_vv_2;
 sum_e = sum_e - ui*sum_e_1 + uiui*sum_e_2;
 sum_H2 = sum_H2 - ui*sum_H2_1 + uiui*sum_H2_2;

end;				% End of group calculations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M-STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
al1_old = alpha1;
al2_old = alpha2;
b_old = b;
su_old = sig2u;
Lam_old = Lam0;
SB_old = SigBeta;
se_old = sig2e;

%fprintf('\n');

alpha1 = inv(sum_XtX) * sum_al1;
sig2u = inv_n * sum_b1;

% Esimate for cummulative baseline hazard
for t=1:nt
  Lam0(t) = sum((Ti<=Tall(t)).*delta./sum_lam(indext));
end

%% Parameter Estimation via Newton-Raphson %%
% Fixed effects estimates
minfn = '-(P1*x - sum(P2.*exp(P3*x)))';
opts  = optimset('TolFun',1e-5,'LevenbergMarquardt','off','Display','off'); 
gamma = [alpha2; b];
gamma = fminsearch(minfn,gamma,opts,[sum_x2 sum_S2],Lam0(indext)',[x2 Eu]);
alpha2 = gamma(1:q);
b = gamma(q+1:q+nb);


% Random effects covariance estimates
sum_vv = reshape(sum_vv,vecSB,1);
vecD = reshape(SigBeta/sig2e,vecSB,1);
sum_U2 = reshape(sum_U2,vecSB,1);
U = [(sum_ni - inv(sig2e)*sum_e); (sig2e*sum_U2 - sum_vv)];
temp4 = -sum_H1 + sig2e*sum_H2;
H = inv(sig2e)*[(-sum_ni + 2*inv(sig2e)*sum_e) sum_vv'; sum_vv temp4];
invHessian = inv(H);
para = [sig2e; vecD] - invHessian*U;
sig2e = para(1);
vecD = para(2:vecSB+1);
SigBeta = reshape(vecD,dimSB,dimSB);	      % Est. for D = SigBeta/sig2e.
SigBeta = (SigBeta + SigBeta')/2;	      % Ensuring matrix symmetric.
SigBeta = sig2e*SigBeta;		      % Adjust by sig2e to get SigBeta.


%%%%%%%%%%%%%%%%%%%% Display & Convergence Checking         %%%%%%%%%%%%%%%%%%

if ind==0;
  fprintf('  %0.4f ',alpha1(1,1));
  fprintf(' %0.4f ',SigBeta(1,:));
  fprintf(' %0.4f  %0.4f  %0.4f  %0.4f \n ',sig2u,sig2e,alpha2(1,1),b);
  for i=2:max([size(alpha1,1) size(SigBeta,1) size(alpha2,1)]);
    if i<=size(alpha1,1);
      fprintf('    %0.4f ',alpha1(i,1));
    else
      fprintf('           ');
    end
    if i<=size(SigBeta,1);
      fprintf(' %0.4f ',SigBeta(i,:));
    else
      fprintf('                 ');
    end;
    if i<=size(alpha2,1);
      fprintf('                   %0.4f \n ',alpha2(i,1));
    else
      fprintf('\n ');
    end
  end
  fprintf('\n');
end;


%%% Check for Convergence %%%
c_sum = 0;
if relerr(alpha1,al1_old,tol);  c_sum = c_sum+1; end;
if relerr(alpha2,al2_old,tol);  c_sum = c_sum+1; end;
if relerr(b,b_old,tol);         c_sum = c_sum+1; end;
if relerr(SigBeta,SB_old,tol);  c_sum = c_sum+1; end;
if relerr(sig2u,su_old,tol);    c_sum = c_sum+1; end;
if relerr(sig2e,se_old,tol);    c_sum = c_sum+1; end;
%if relerr(Lam0,Lam_old,tol);    c_sum = c_sum+1; end;

if c_sum==6;
  fprintf('\n');
  disp(['Convergence achieved after ' num2str(iter) ' iterations']);
  conver = 1;
end;
if iter==maxit;
  fprintf('\n');
  disp('Maximum number of iterations reached');
  conver = 1;
end;


end		% End of while loop

fprintf('\n');

out.alpha1 = alpha1;
out.alpha2 = alpha2;
out.SigBeta = SigBeta;
out.b = b;
out.sig2u = sig2u;
out.u = u;
out.sig2e = sig2e;
out.Lambda0 = Lam0;
out.iterations = iter;

subject=0;
nstart=0;
errY = [];
errS = [];
for e=1:n;
   subject = subject + 1;
   nij = ni(subject);
   yij = Y(nstart+1:nstart+nij,1);
   Xij = X(nstart+1:nstart+nij,:);
   Zij = Z(nstart+1:nstart+nij,:);
   Tij = Ti(subject);
   x2ij = x2(subject,:)';
   ui  = u(patt(subject));
   nstart = nstart + nij;
   errys = yij - Xij*alpha1 - ui;
   Vinv = inv(eye(nij)*sig2e + Zij*SigBeta*Zij');
   bij = SigBeta*Zij'*Vinv*errys;
   errY = [errY; errys - Zij*bij];
   lamind = sum(Tij>=Tall);
   errS = [errS; Lam0(lamind)*exp(x2ij'*alpha2 + b*ui) + (1-delta(subject))];
end;
resid.errY = errY;
resid.errS = errS;
resid.fail = delta;
resid.patt = patt;

%%%%%%%%%% Calculation of Standard Errors for the parameter estimates %%%%%%%%%%

% For parameters calculated using Newton-Raphson, std.errs can be estimated
% from the inversion Hession matrix.

nr = sqrt(diag(invHessian));
sqrtn = sqrt(n);
stderr.sig2e = nr(1)/sqrtn;
stderr.SigBeta = reshape(nr(2:vecSB+1),dimSB,dimSB)/sqrtn;


% Std. errors for the fixed effects
if stdcal==1;
  nu = 1/n;
  save survstdinit;
  disp('Standard errors calculations');
 % Alpha1 
  nal1 = length(alpha1);
  temp = nu*eye(nal1);
  c = 1;
  for ele = 1:nal1;
    stdest = surv_std(alpha1+temp(:,ele),alpha2,b);
    SS(:,c) = stdest;
    c = c+1;
  end
 % Alpha2
  nal2 = length(alpha2);
  temp = nu*eye(nal2);
  for ele=1:nal2;
    stdest = surv_std(alpha1,alpha2+temp(:,ele),b);
    SS(:,c) = stdest;
    c = c+1;
  end
 % b
  stdest = surv_std(alpha1,alpha2,b+nu);
  SS(:,c) = stdest;
  SS = SS/nu;
  stdest = inv(-SS);
  stdest = diag(stdest);
  stderr.alpha1 = stdest(1:nal1);
  stderr.alpha2 = stdest(nal1+1:nal1+nal2);
  stderr.b = stdest(nal1+nal2+1);
SS2 = SS/sum_ni;
SS2 = diag(inv(SS2));
std2.alpha1 = SS2(1:nal1);
std2.alpha2 = SS2(nal1+1:nal1+nal2);
std2.b = SS2(nal1+nal2+1);
end;
%stdall = sqrt(diag(inv(ucu)));
%stderr.alpha1 = stdall(1:length(alpha1))/sqrtn;
%stderr.alpha2 = stdall(length(alpha1)+1:length(stdall))/sqrtn;


if ind < 3
  figure(1)
  title('Baseline Survivor Function');
  xlabel('Time');
  hold on;
  plot([0; Tall],[1 exp(-Lam0)],'b');
  if ind==0
    legend('Initial Est.','Final Est.');
  end
  hold off
end
