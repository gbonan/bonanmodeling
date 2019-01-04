function [esat, desat] = satvap (tc)

% Compute saturation vapor pressure and change in saturation vapor pressure
% with respect to temperature. Polynomial approximations are from:
% Flatau et al. (1992) Polynomial fits to saturation vapor pressure.
% Journal of Applied Meteorology 31:1507-1513. Input temperature is Celsius.

% --- For water vapor (temperature range is 0C to 100C)
 
a0 =  6.11213476;        b0 =  0.444017302;
a1 =  0.444007856;       b1 =  0.286064092e-01;
a2 =  0.143064234e-01;   b2 =  0.794683137e-03;
a3 =  0.264461437e-03;   b3 =  0.121211669e-04;
a4 =  0.305903558e-05;   b4 =  0.103354611e-06;
a5 =  0.196237241e-07;   b5 =  0.404125005e-09;
a6 =  0.892344772e-10;   b6 = -0.788037859e-12;
a7 = -0.373208410e-12;   b7 = -0.114596802e-13;
a8 =  0.209339997e-15;   b8 =  0.381294516e-16;
 
% --- For ice (temperature range is -75C to 0C)
 
c0 =  6.11123516;        d0 =  0.503277922;
c1 =  0.503109514;       d1 =  0.377289173e-01;
c2 =  0.188369801e-01;   d2 =  0.126801703e-02;
c3 =  0.420547422e-03;   d3 =  0.249468427e-04;
c4 =  0.614396778e-05;   d4 =  0.313703411e-06;
c5 =  0.602780717e-07;   d5 =  0.257180651e-08;
c6 =  0.387940929e-09;   d6 =  0.133268878e-10;
c7 =  0.149436277e-11;   d7 =  0.394116744e-13;
c8 =  0.262655803e-14;   d8 =  0.498070196e-16;

% --- Limit temperature to -75C to 100C
 
tc = min(tc, 100);
tc = max(tc, -75);

% --- Saturation vapor pressure (esat, mb) and derivative (desat, mb)

if (tc >= 0)
   esat  = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 ...
         + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8)))))));
   desat = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 ...
         + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8)))))));
else
   esat  = c0 + tc*(c1 + tc*(c2 + tc*(c3 + tc*(c4 ...
         + tc*(c5 + tc*(c6 + tc*(c7 + tc*c8)))))));
   desat = d0 + tc*(d1 + tc*(d2 + tc*(d3 + tc*(d4 ...
         + tc*(d5 + tc*(d6 + tc*(d7 + tc*d8)))))));
end

% --- Convert from mb to Pa

esat  = esat  * 100;
desat = desat * 100;
