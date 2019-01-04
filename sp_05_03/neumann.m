function [dummy] = neumann

% --------------------------------------------------------
% Analytical solution for Neumann problem (Lunardini 1981)
% --------------------------------------------------------

% soil depths (meters)

nsoi = 4;
depth(1) = 0.25;
depth(2) = 0.55;
depth(3) = 0.85;
depth(4) = 1.15;

% number of days to simulate

ndays = 60;

% time-invariant surface temperature (deg C)

ts = -10;

% initial soil temperature (deg C)

for i = 1:nsoi
   t0(i) = 2;
end

% freezing temperature (deg C)

tf = 0;

% soil dependent constants
% a1 - frozen diffusivity (cm**2/hr -> m**2/s)
% a2 - unfrozen diffusivity (cm**2/hr -> m**2/s)
% rm - "m" from Jumikis (1966) and adjust for above unit conversion
% rg - "gamma" from Lunardini (1981)

a1 = 42.55 / 3600 / (100*100);
a2 = 23.39 / 3600 / (100*100);
rm = 3.6 / (60 * 100);
rg = rm / (2 * sqrt(a1));

% save output for day zero (i.e., initial conditions)

m = 1;
iday_out(m) = 0;
xfr_out(m) = 0;
t1_out(m) = t0(1);
t2_out(m) = t0(2);
t3_out(m) = t0(3);
t4_out(m) = t0(4);
z1_out(m) = depth(1) * 100;
z2_out(m) = depth(2) * 100;
z3_out(m) = depth(3) * 100;
z4_out(m) = depth(4) * 100;

% begin time loop

for iday = 1:ndays

   % time (seconds)

   time = iday * 24 * 3600;

   % depth of frost penetration XF (meters) based on time (seconds)

   xf = 2 * rg * sqrt(a1*time);

   % calculate soil temperatures given XF and time

   for i = 1:nsoi
      if (depth(i) <= xf)
         x = depth(i) / (2*sqrt(a1*time));
         y = rg;
         t(i) = ts + (tf - ts) * erf(x) / erf(y);
      else
         x = depth(i) / (2*sqrt(a2*time));
         y = rg * sqrt(a1/a2);
         t(i) = t0(i) - (t0(i) - tf) * (1-erf(x)) / (1-erf(y));
      end
   end

   % save for output

   m = m + 1;
   iday_out(m) = iday;
   xf_out(m) = xf * 100;
   t1_out(m) = t(1);
   t2_out(m) = t(2);
   t3_out(m) = t(3);
   t4_out(m) = t(4);
   z1_out(m) = depth(1) * 100;
   z2_out(m) = depth(2) * 100;
   z3_out(m) = depth(3) * 100;
   z4_out(m) = depth(4) * 100;

end

A = [iday_out; xf_out; z1_out; t1_out; z2_out; t2_out; z3_out; t3_out; z4_out; t4_out];
fileID = fopen('data_analytical.txt','w');
fprintf(fileID,'%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n','day','z0C','z1','t1','z2','t2','z3','t3','z4','t4');
fprintf(fileID,'%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', A);
fclose(fileID);

dummy = 0;
