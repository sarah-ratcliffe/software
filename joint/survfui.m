function [SS,const] = survfui(ind,limits,const_ft,fui,b,sig2u,nid,prod_ft,sum_ft,prod_fy,sum_fy1,sum_fy2,sum_fy3,sum_fy4,const_fu,exp_bui);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Function used by surv.m to calculate integrals via
% Gaussian Integration. (from lowlim to uplim)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

npts = 48;
load gauss48pts.dat;
X = gauss48pts(:,1);
W = gauss48pts(:,2);
clear gauss48pts;
lowlim = limits(1);
uplim = limits(2);

SS=0;
D1 = (uplim-lowlim)/2;
D2 = (uplim+lowlim)/2;

CC=prod_ft*prod_fy;
DD=sum_fy2+sum_fy3;
const = 1e-6;

xgrid = [];
ygrid1 = [];
ygrid2 = [];
ygrid3 = [];
ygrid4 = [];

for j=1:npts
 u1= D1*X(j)+D2;
 u2= -D1*X(j)+D2;

 ft = (exp(b*u1*nid)/const_ft)*exp(-sum_ft*(exp(b*u1)-const_fu));
 fy = exp(prod_fy*(-sum_fy1/2 + u1*sum_fy2 -1/2*(u1^2)*sum_fy4));
 fu = exp(-u1^2);
 F1 = ft*fy*fu;
%fprintf('\n%1.2f  %1.2f ',-sum_ft*(exp(b*u1)-const_ft),-sum_ft*(exp(b*u1)-const_fu));
%fprintf('\n%1.16f %1.16f %1.16f %1.16f',ft,fy,fu,F1);
%disp([ft,fy,fu,F1]);

 ft = (exp(b*u2*nid)/const_ft)*exp(-sum_ft*(exp(b*u2)-const_fu));
 fy = exp(prod_fy*(-sum_fy1/2 + u2*sum_fy2 -1/2*(u2^2)*sum_fy4));
 fu = exp(-u2^2);
 F2 = ft*fy*fu;
%fprintf('\n%1.2f  %1.2f ',-sum_ft*(exp(b*u2)-const_ft),-sum_ft*(exp(b*u1)-const_fu));
%fprintf('\n%1.16f %1.16f %1.16f %1.16f',ft,fy,fu,F2);
%disp([ft,fy,fu,F1]);

 xgrid = [u1; xgrid; u2];
 ygrid1 = [F1; ygrid1; F2];

 switch ind
  case 1		% Calculate integral in f(u_i)
    FF1 = F1;
    FF2 = F2;

  case 2		% Calculate E(u_i)
    gu1 = u1;
    FF1 = gu1*F1/fui;
    gu2 = u2;
    FF2 = gu2*F2/fui;
 ygrid2 = [FF1; ygrid2; FF2];

  case 3  		% Calculate E(u_i^T u_i)
    gu1 = u1^2;
    FF1 = gu1*F1/fui;
    gu2 = u2^2;
    FF2 = gu2*F2/fui;
 ygrid3 = [FF1; ygrid3; FF2];

  case 4		% Calculate E[e^(b u_i)]
    gu1 = exp(b*u1);
    FF1 = gu1*F1/fui;
    gu2 = exp(b*u2);
    FF2 = gu2*F2/fui;
 ygrid4 = [FF1; ygrid4; FF2];

  case 5		% Calculate E[u_i e^(b u_i)]
    gu1 = u1*exp(b*u1);
    FF1 = gu1*F1/fui;
    gu2 = u2*exp(b*u2);
    FF2 = gu2*F2/fui;

  case 6		% Calculate E[u_i^2 e^(b u_i)]
    gu1 = u1^2 *exp(b*u1);
    FF1 = gu1*F1/fui;
    gu2 = u2^2 *exp(b*u2);
    FF2 = gu2*F2/fui;

  otherwise	% Not a valid option
    disp('ERROR: invalid integral function to evaluate');
    return;
 end;

 SS=SS+W(j)*(FF1+FF2);
end;

SS = const*SS*D1;

return;
if ind==1
  figure;
  subplot(2,2,1);
  plot(xgrid,ygrid1);
  title('h(u_i)');
elseif ind==2
  subplot(2,2,2);
  plot(xgrid,ygrid2);
  title('u_i h(u_i)');
elseif ind==3
  subplot(2,2,3);
  plot(xgrid,ygrid3);
  title('u_i^T u_i h(u_i)');
elseif ind==4
  subplot(2,2,4);
  plot(xgrid,ygrid4);
  title('exp(b u_i) h(u_i)');
else
end

