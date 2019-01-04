% Supplemental program 9.1

% -------------------------------------------------
% Use Green-Ampt equation to calculate infiltration
% -------------------------------------------------

% --- Input parameters

% sandy loam

  Ksat = 12.48 * 10 / 3600;   % cm/h -> mm/s
  psi_sat = -218;             % mm
  theta_sat = 0.435;

% loam

% Ksat = 2.5 * 10 / 3600;   % cm/h -> mm/s
% psi_sat = -478;           % mm
% theta_sat = 0.451;

% clay loam

% Ksat = 0.88 * 10 / 3600;   % cm/h -> mm/s
% psi_sat = -630;            % mm
% theta_sat = 0.476;

% Initial soil moisture

theta_dry = 0.3;
delta_theta = theta_sat - theta_dry;

% Matrix potential at wetting front

psi_w = 0.76 * psi_sat;

% --- Solve Green-Ampt equation for cumulative infiltration at time t

err = 1e-06;

% Time step (s)

dt = 0.1;

% Loop through one hour

for i = 1:36000
   t = i * dt;

   fval = @(x) x - (Ksat * t + abs(psi_w) * delta_theta * log(1 + x / (abs(psi_w)*delta_theta)));

   % Use bisection to sovle for I

   aval = 1e03;
   bval = 0;

   fa = fval(aval);
   fb = fval(bval);

   if (fa * fb > 0)
      error('bisection error: f(a) and f(b) do not have opposite signs')
   end

   while (abs(bval-aval) > err)
      cval = (aval + bval) / 2;
      fc = fval(cval);
      if (fa * fc < 0)
         bval = cval; fb = fc;
      else
         aval = cval; fa = fc;
      end
   end

   % Time (minutes)

   time(i) = t / 60;

   % Cumulative infiltration (mm)

   I_cum(i) = cval;

   % Infiltration rate (mm/s)

   i_rate(i) =  Ksat * (abs(psi_w) * delta_theta / I_cum(i) + 1);

end

fprintf('I = %12.8f\n',I_cum(36000))

% --- Write output file

B = [time; i_rate; I_cum];
fileID = fopen('data.txt','w');
fprintf(fileID,'%12s %12s %12s\n','time','i','I');
fprintf(fileID,'%12.6f %12.6f %12.3f\n', B);
fclose(fileID);

% --- Make graph

%figure
semilogy(time,i_rate)
%plot(time,i_rate)
title('Green-Ampt infiltration')
xlabel('Time (minutes)')
ylabel('Infiltration rate (mm s^{-1})')
