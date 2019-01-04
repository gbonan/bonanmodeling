% Supplemental program 2.1

% ----------------------------------------------
% Calculate and graph leaf area density profiles
% ----------------------------------------------

% --- Parameters for beta distribution

p = [2.5, 3.5, 11.5];
q = [2.5, 2.0,  3.5];

% --- Canopy parameters

LAI = 5;              % leaf area index (m2/m2)
hc = 10;              % canopy height (m)

fprintf('Leaf area index = %6.2f\n',LAI)

% --- Create a vector of heights (z) with linearly spaced
% values from z_min to z_max with increment dz

z_min = 0;            % minimum
z_max = hc;           % maximum
dz = 0.1;             % increment
z = z_min:dz:z_max;   % z is linearly spaced between z_min and z_max with increment dz

% --- Calculate the leaf area density profile for each [p,q]

% Loop over each p,q

for i = 1:length(p)

   % Loop over each height

   sum = 0;
   for j = 1:length(z)
      x = z(j) / hc;
      lad = (LAI / hc) * (x^(p(i)-1) * (1 - x)^(q(i)-1)) / beta(p(i),q(i));

      % Numerically sum leaf area for each height

      sum = sum + lad * dz;

      % Save output for graphing

      if (i == 1)
         y1(j) = lad;
      elseif (i == 2)
         y2(j) = lad;
      elseif (i == 3)
         y3(j) = lad;
      end
   end

   fprintf('p, q = %6.2f %6.2f\n',p(i),q(i))
   fprintf('Leaf area index (numerical) = %8.4f\n',sum)

end

% --- Make graph for leaf area density in relation to relative height (z/hc)

z_rel = z ./ hc;

plot(y1,z_rel,'b-',y2,z_rel,'r-',y3,z_rel,'g-')
title('Profiles')
xlabel('Leaf area density (m^2 m^{-3})')
ylabel('Height (z/h_c)')
legend('p,q = 2.5,2.5','p,q = 3.5,2.0','p,q = 11.5,3.5','Location','southeast')
