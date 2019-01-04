function [fluxvar, root] = hybrid_root (func, physcon, forcvar, surfvar, fluxvar, xa, xb, tol)

% Solve for the root of a function using the secant and Brent's methods given
% initial estimates xa and xb. The root is updated until its accuracy is tol. 
% func is the name of the function to solve. The variable root is returned as
% the root of the function. The function being evaluated has the definition statement: 
%
% function [fluxvar, fx] = func (physcon, forcvar, surfvar, fluxvar, x)
%
% The function func is exaluated at x and the returned value is fx. It uses variables
% in the physcon, forcvar, surfvar, and fluxvar structures. These are passed in as
% input arguments. It also calculates values for variables in the fluxvar structure
% so this must be returned in the function call as an output argument. The matlab
% function feval evaluates func.

% --- Evaluate func at xa and see if this is the root

x0 = xa;
[fluxvar, f0] = feval(func, physcon, forcvar, surfvar, fluxvar, x0);
if (f0 == 0)
   root = x0;
   return
end

% --- Evaluate func at xb and see if this is the root

x1 = xb;
[fluxvar, f1] = feval(func, physcon, forcvar, surfvar, fluxvar, x1);
if (f1 == 0)
   root = x1;
   return
end

% --- Order initial root estimates correctly

if (f1 < f0)
   minx = x1;
   minf = f1;
else
   minx = x0;
   minf = f0;
end

% --- Iterative root calculation. Use the secant method, with Brent's method as a backup

itmax = 40;
for iter = 1:itmax
   dx = -f1 * (x1 - x0) / (f1 - f0);
   x = x1 + dx;

   % Check if x is the root. If so, exit the iteration

   if (abs(dx) < tol)
      x0 = x;
      break
   end

   % Evaluate the function at x

   x0 = x1;
   f0 = f1;
   x1 = x;
   [fluxvar, f1] = feval(func, physcon, forcvar, surfvar, fluxvar, x1);
   if (f1 < minf)
      minx = x1;
      minf = f1;
   end

   % If a root zone is found, use Brent's method for a robust backup strategy
   % and exit the iteration

   if (f1 * f0 < 0)
      [fluxvar, x] = brent_root (func, physcon, forcvar, surfvar, fluxvar, x0, x1, tol);
      x0 = x;
      break
   end

   % In case of failing to converge within itmax iterations stop at the minimum function

   if (iter == itmax)
      [fluxvar, f1] = feval(func, physcon, forcvar, surfvar, fluxvar, minx);
      x0 = minx;
   end

end

root = x0;
