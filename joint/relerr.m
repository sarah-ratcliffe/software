function x = relerr(new,old,tol);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Last Updated: 1 Feb 2001
%
% COMMAND  : x = relerr(new,old,tol); 
%  ACTION  : For two estimates of a value/vector of values (new
%             and old), relerr indicates if the relative error
%             of the two values (w.r.t new) is less than a given
%             tolerance level (tol). New and old must have the
%             same dimensions, and conparisons are on a point
%             by point basis.
%
%   INPUT  : new = new estimate.
%            old = old estimate.
%            tol = tolerance for error  [default = 1e-6].
%
%  OUTPUT  :   x = indictor for error less than tolerance
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin<3; tol = 1e-6; end;
if nargin<2;
  disp('ERROR: x = relerr(new,old,tol)');
  return;
end;

[r,c] = size(new);

if new==0
  err = abs(new - old);
else
  err = abs((new - old) ./ new);
end

compar = err < tol;

x = (sum(sum(compar)) == r*c);
