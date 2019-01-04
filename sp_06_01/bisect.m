function [c] = bisect (func_name, a, b, delta, var)

% -----------------------------------------------------------------
% Use the bisection method to find the root of a function f
% between a and b. The root is refined until its accuracy is delta.
%
% Input:  func_name  ! Name of the function to solve
%         a          ! Low endpoint of the interval
%         b          ! High endpoint of the interval
%         delta      ! Tolerance/accuracy
%         var        ! Input variables for function
% Output: c          ! Root
% -----------------------------------------------------------------

% Evaluate function at a and b

fa = feval(func_name, a, var);
fb = feval(func_name, b, var);

% Error check: root must be bracketed

if (sign(fa) == sign(fb))
   error('bisect error: f must have different signs at the endpoints a and b')
end

% Iterate to find root

while (abs(b - a) > 2*delta)
   c = (b + a)/2;
   fc = feval(func_name, c, var);
   if (sign(fc) ~= sign(fb))
      a = c; fa = fc;
   else
      b = c; fb = fc;
   end
end
