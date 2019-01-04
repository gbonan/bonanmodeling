% Supplemental program 14.1

% -----------------------------------------------------------------------------
% Use Eqs. (14.25) and (14.26) to numerically integrate the function G(Z) over
% nine leaf inclination angles (in increments of 10 degrees) to calculate the
% direct beam extinction coefficient Kb. Leaf azimuth angles are randomly
% distributed (uniform PDF). Leaf angle classes are (5, 15, ..., 75, 85).
% -----------------------------------------------------------------------------

% --- Solar elevation angle and zenith angle

solar_elevation = 60;                               % degrees
solar_zenith = (90 - solar_elevation) * (pi / 180); % radians

% --- The variable "leaf" defines the leaf angle distribution

  leaf = 'Planophile';
% leaf = 'Erectophile';
% leaf = 'Plagiophile';
% leaf = 'Uniform';
% leaf = 'Spherical';

% leaf = 'Ellipsoidal - spherical';
% leaf = 'Ellipsoidal - horizontal';
% leaf = 'Ellipsoidal - vertical';

% --- Parameters for ellipsoidal leaf angle distribution

xe = 0;
switch leaf
   case 'Ellipsoidal - spherical'
   xe = 1;
   case 'Ellipsoidal - horizontal'
   xe = 100;
   case 'Ellipsoidal - vertical'
   xe = 0.01 ;
end

if (xe == 1)
   le = 2;
elseif (xe < 1)
   ee = sqrt(1 - xe*xe);
   le = xe + asin(ee) / ee;
elseif (xe > 1)
   ee = sqrt(1 - 1/(xe*xe));
   le = xe + log((1 + ee)/(1 - ee)) / (2 * ee * xe);
end

% --- Create leaf angle probability density function (leaf_pdf) in 0.1 degree increments.
% This is needed to calculate the fractional abundance of leaves in each 10 degree bin.

dleaf_angle = 0.1 * (pi / 180);     % Leaf angle increment (0.1 deg -> radians)

for i = 1:900
   leaf_pdf(i) = 0;
   ang1 = (i-1) * dleaf_angle;
   ang2 =  i    * dleaf_angle;
   angle = (ang1 + ang2) / 2;

   switch leaf
      case 'Planophile'
      leaf_pdf(i) = 2 / pi * (1 + cos(2*angle));

      case 'Erectophile'
      leaf_pdf(i) = 2 / pi * (1 - cos(2*angle));

      case 'Plagiophile'
      leaf_pdf(i) = 2 / pi * (1 - cos(4*angle));

      case 'Uniform'
      leaf_pdf(i) = 2 / pi;

      case 'Spherical'
      leaf_pdf(i) = sin(angle);

      case {'Ellipsoidal - spherical', 'Ellipsoidal - horizontal', 'Ellipsoidal - vertical'}
      leaf_pdf(i) = 2 * xe^3 * sin(angle) / (le * (cos(angle)^2 + xe^2 * sin(angle)^2)^2);
   end
end

% Average leaf angle

sum = 0;
ave = 0;
for i = 1:900
   ang1 = (i-1) * dleaf_angle;
   ang2 =  i    * dleaf_angle;
   angle = (ang1 + ang2) / 2;
   sum = sum + leaf_pdf(i) * dleaf_angle;
   ave = ave + angle * leaf_pdf(i) * dleaf_angle;
end

fprintf(' \n')
fprintf('Leaf type = %30s\n', leaf)
fprintf('Sum of leaf angle distribution = %15.4f\n', sum)
fprintf('Mean leaf angle = %15.4f\n', ave*180/pi)

% Relative leaf angle distribution (lad) in 9 10-degree bins (5, 15, ..., 75, 85).
% This is the fraction of leaves in each 10-degree bin.

for j = 1:9                      % 10-degree bin
   lad(j) = 0;
   for i = (j-1)*100+1:j*100     % Leaf angles that fall in this bin
      lad(j) = lad(j) + leaf_pdf(i) * dleaf_angle;
   end
end

% --- Numerically integrate G(Z) over the nine leaf inclination angles

sum = 0;
for i = 1:9
   ang1 = (i-1) * (10 * pi / 180);
   ang2 =  i    * (10 * pi / 180);
   leaf_angle = (ang1 + ang2) / 2;

   a = cos(solar_zenith) * cos(leaf_angle);
   b = sin(solar_zenith) * sin(leaf_angle);

   if (leaf_angle <= (pi/2-solar_zenith))
      Gfunc_i = a;
   else
      c = sqrt( sin(leaf_angle)^2 - cos(solar_zenith)^2 );
      Gfunc_i = 2 / pi * (a * asin(a/b) + c);
   end

   sum = sum + Gfunc_i * lad(i);
end

% G(Z) - projected leaf area normal to solar beam

Gfunc = sum;

% --- Extinction coefficient for direct beam

Kb = Gfunc / cos(solar_zenith);

% --- Now calculate G(Z) using Ross index

F1 = lad(1) + lad(2) + lad(3);
F2 = lad(4) + lad(5) + lad(6);
F3 = lad(7) + lad(8) + lad(9);
xl = 0.5 * (abs(0.134-F1) + abs(0.366-F2) + abs(0.5-F3));
if ((0.5-F3) < 0)
   xl = -xl;
end

xl = min(max(xl, -0.4), 0.6); 
phi1 = 0.5 - 0.633 * xl - 0.330 * xl^2;
phi2 = 0.877 * (1 - 2 * phi1);
Gfunc_xl = phi1 + phi2 * cos(solar_zenith);
Kb_xl = Gfunc_xl / cos(solar_zenith);

% --- Print output

fprintf('Solar zenith = %15.4f\n', solar_zenith*180/pi)
fprintf('G(Z) = %15.4f\n', Gfunc)
fprintf('Kb = %15.4f\n', Kb)
fprintf('Ross index = %15.4f\n', xl)
fprintf('G(Z) = %15.4f\n', Gfunc_xl)
fprintf('Kb = %15.4f\n', Kb_xl)
