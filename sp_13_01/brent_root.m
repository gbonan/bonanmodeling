function [flux, root] = brent_root (func, physcon, atmos, leaf, flux, xa, xb, tol)

% Use Brent's method to find the root of a function, which is known to exist between
% xa and xb. The root is updated until its accuracy is tol. func is the name of the
% function to solve. The variable root is returned as the root of the function. The
% function being evaluated has the definition statement: 
%
% function [flux, fx] = func (physcon, atmos, leaf, flux, x)
%
% The function func is exaluated at x and the returned value is fx. It uses variables
% in the physcon, atmos, leaf, and flux structures. These are passed in as
% input arguments. It also calculates values for variables in the flux structure
% so this must be returned in the function call as an output argument. The matlab
% function feval evaluates func.

% --- Evaluate func at xa and xb and make sure the root is bracketed

a = xa;
b = xb;
[flux, fa] = feval(func, physcon, atmos, leaf, flux, a);
[flux, fb] = feval(func, physcon, atmos, leaf, flux, b);

if ((fa > 0 & fb > 0) | (fa < 0 & fb < 0))
   error('brent_root error: root must be bracketed')
end

% --- Initialize iteration

itmax = 50;      % Maximum number of iterations
eps1 = 1e-08;    % Relative error tolerance

c = b;
fc = fb;

% --- Iterative root calculation

for iter = 1:itmax
   if ((fb > 0 & fc > 0) | (fb < 0 & fc < 0))
      c = a;
      fc = fa;
      d = b - a;
      e = d;
   end
   if (abs(fc) < abs(fb))
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
   end
   tol1 = 2 * eps1 * abs(b) + 0.5 * tol;
   xm = 0.5 * (c - b);

   % Check to end iteration

   if (abs(xm) <= tol1 | fb == 0)
      break
   end

   if (abs(e) >= tol1 & abs(fa) > abs(fb))
      s = fb / fa;
      if (a == c)
         p = 2 * xm * s;
         q = 1 - s;
      else
         q = fa / fc;
         r = fb / fc;
         p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
         q = (q - 1) * (r - 1) * (s - 1);
      end
      if (p > 0)
         q = -q;
      end
      p = abs(p);
      if (2*p < min(3*xm*q-abs(tol1*q), abs(e*q)))
         e = d;
         d = p / q;
      else
         d = xm;
         e = d;
      end
   else
      d = xm;
      e = d;
   end
   a = b;
   fa = fb;
   if (abs(d) > tol1)
      b = b + d;
   else
      if (xm >= 0)
         b = b + abs(tol1);
      else
         b = b - abs(tol1);
      end
   end
   [flux, fb] = feval(func, physcon, atmos, leaf, flux, b);

   % Check to end iteration

   if (fb == 0)
      break
   end

   % Check to see if failed to converge

   if (iter == itmax)
      error('brent_root error: Maximum number of interations exceeded')
   end

end
root = b;
