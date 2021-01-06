function stdest = surv_std(alpha1fix,alpha2fix,bfix);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe  (copyright 2002)
%
% COMMAND: stdest = surv_std(stdind,ele,fixest);
%
% ACTION: Calculate profile likelihood standard errors for surv.m.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load survstdinit;

tol = 10^(log10(tol)/2);
maxit = 5;
iter = 0;
nal1 = length(alpha1);
nal2 = length(alpha2);
al2const = zeros(n,1);
%%%%%%%%%%%%%%%%%%%%%           ALGORITHM START            %%%%%%%%%%%%%%%%%%
conver = 0;
alpha1 = alpha1fix;
alpha2 = alpha2fix;
b = bfix;
fprintf('\nIteration: ');

while conver==0		% Do until convergence or max iterations
  iter = iter + 1;
  fprintf(' %1.0f',iter);

%%%%%%%%%%%%%%%%%%%%% STEP 1:                               %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Parameter Estimation via EM algorithm %%%%%%%%%%%%%%%%%%
% Reset sums to zero for next iteration
%sum_plot = zeros(2,41);
lam01 = delta./sum_lam(indext);
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
Eueu = [];

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
   ex2ij = exp(x2ij'*alpha2 + al2const(subject));

   Qij = Zij*SigBeta*Zij' + sig2e*eye(nij);
   invQij = inv(Qij);

   % Standard Errors for fixed effects 
   SigXij = invQij*Xij;

   % Sums used in M-step  
   sum_XtX  = sum_XtX + Xij'*SigXij;
   sum_al1  = sum_al1 + SigXij'*yij;
   sum_X1ti = sum_X1ti + SigXij'*ones(nij,1);
   sum_lami = sum_lami + ((Tij >= Tall)*ex2ij); 

   % Sums used in Newton-Raphson for sig2e and sigb.
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
 prod_fy = inv(-sum_fy1/2 + 0.5*sum_fy2 -0.125*sum_fy4)*prod_fy;

 fui = survfui(1,limits,const_ft,0,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);
%fprintf(' fui = %0.16f \n',fui);
 ui = survfui(2,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(u_i | ...)
 uiui = survfui(3,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(u_i^2 | ...)
 Eeui = survfui(4,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(exp{b u_i}|...)
 Eueui = survfui(5,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);	% E(u_i exp{b u_i}|...)
 Eeu = [Eeu; ones(mi(i),1)*Eeui];
 Eu = [Eu; ones(mi(i),1)*ui];
 Eueu = [Eueu; ones(mi(i),1)*Eueui];

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

%alpha1 = inv(sum_XtX) * sum_al1;
sig2u = inv_n * sum_b1;

% Esimate for cummulative baseline hazard
for t=1:nt
  Lam0(t) = sum((Ti<=Tall(t)).*delta./sum_lam(indext));
end

%% Parameter Estimation via Newton-Raphson %%
% Fixed effects estimates
%  minfn = '-(P1*x - sum(P2.*exp(P3*x)))';
%  opts  = optimset('TolFun',1e-5,'LevenbergMarquardt','off','Display','off'); 
%  gamma = [alpha2; b];
%  gamma = fminsearch(minfn,gamma,opts,[sum_x2 sum_S2],Lam0(indext)',[x2 Eu]);
%  alpha2 = gamma(1:q);
%  b = gamma(q+1:q+nb);

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
SigBeta = abs(SigBeta + SigBeta')/2;	      % Ensuring matrix symmetric.
SigBeta = sig2e*SigBeta;		      % Adjust by sig2e to get SigBeta.


%%%%%%%%%%%%%%%%%%%% Display & Convergence Checking         %%%%%%%%%%%%%%%%%%

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

nal2 = length(alpha2);

stda1 = sum_al1 - sum_XtX*alpha1;
stda2 = sum_x2 - sum(x2.*((Lam0(indext)'.*Eeu)*ones(1,nal2)));
stdb  = sum_S2 - sum(Lam0(indext)'.*Eueu);

stdest = [stda1; stda2'; b];
