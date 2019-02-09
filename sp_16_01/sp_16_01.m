% Supplemental program 16.1

% This code replicates the multilayer canopy model of Bonan et al. (2018)
% at US-UMB.2006 for the first time step

% --- Physical constants

physcon.vkc = 0.4;                  % von Karman constant
physcon.grav = 9.80665;             % Gravitational acceleration (m/s2)
physcon.tfrz = 273.15;              % Freezing point of water (K)
physcon.hvap = 2.501e6;             % Latent heat of evaporation (J/kg)
physcon.mmh2o = 18.02 / 1000;       % Molecular mass of water (kg/mol)
physcon.mmdry = 28.97 / 1000;       % Molecular mass of dry air (kg/mol)
physcon.cpd = 1005;                 % Specific heat of dry air at constant pressure (J/kg/K)
physcon.cpw = 1846;                 % Specific heat of water vapor at constant pressure (J/kg/K)
physcon.rgas = 8.31446;             % Universal gas constant (J/K/mol)

% --- Array indices for leaves

leafvar.nleaf = 2;                  % Two leaf types (sunlit and shaded)
leafvar.isun = 1;                   % Sunlit leaf
leafvar.isha = 2;                   % Shaded leaf

% --- Define canopy structure

npts = 1;                           % Number of grid points

for p = 1:npts

   % Plant area index of canopy (m2/m2)

   surfvar.pai(p) = 5.051612734794617;

   % Atmospheric forcing reference height (m)

   forcvar.zref(p) = 46;

   % Canopy height (m)

   surfvar.hc(p) = 21;

   % Determine number of within-canopy layers by specifying height increment (m)

   dht = 0.5;

   % Define the number of within-canopy layers and adjust dht if needed

   nveg = round(surfvar.hc(p) / dht);
   dht = surfvar.hc(p) / nveg;

   % Set array indices for within canopy layers

   surfvar.nsoi(p) = 1;                                     % First layer is soil
   surfvar.nbot(p) = surfvar.nsoi(p) + 1;                   % Bottom leaf layer
   surfvar.ntop(p) = surfvar.nbot(p) + nveg - 1;            % Top leaf layer

   % Calculate heights at layer interfaces (zw). These are the heights
   % for the conductance between two scalar concentrations. They are
   % defined for ic = nsoi (ground) to ic = ntop (top of the canopy).

   ic = surfvar.ntop(p);
   surfvar.zw(p,ic) = surfvar.hc(p);
   for ic = surfvar.ntop(p)-1: -1: surfvar.nsoi(p)
      surfvar.zw(p,ic) = surfvar.zw(p,ic+1) - dht;
   end

   ic = surfvar.nsoi(p);
   if (surfvar.zw(p,ic) > 1e-10 || surfvar.zw(p,ic) < 0)
      error('zw improperly defined at ground level')
   end

   % Now calculate the above-canopy layers and their heights

   dz_to_zref = forcvar.zref(p) - surfvar.hc(p);
   n_to_zref = round(dz_to_zref / dht);
   dht = dz_to_zref / n_to_zref;
   surfvar.nlev(p) = surfvar.ntop(p) + n_to_zref;

   ic = surfvar.nlev(p);
   surfvar.zw(p,ic) = forcvar.zref(p);
   for ic = surfvar.nlev(p)-1: -1: surfvar.ntop(p)+1
      surfvar.zw(p,ic) = surfvar.zw(p,ic+1) - dht;
   end

   % Determine heights of the scalar concentration and scalar source
   % (zs). These are physically centered between the conductance points
   % (i.e., in the middle of the layer).

   ic = surfvar.nsoi(p);
   surfvar.zs(p,ic) = 0;
   for ic = surfvar.nbot(p):surfvar.nlev(p)
      surfvar.zs(p,ic) = 0.5 * (surfvar.zw(p,ic) + surfvar.zw(p,ic-1));
   end

   % Determine plant area index increment for each layer by numerically
   % integrating the plant area density (beta distribution) between
   % the bottom and top heights for that layer

   pbeta = 3.5;     % Parameter for beta distribution
   qbeta = 2.0;     % Parameter for beta distribution

   for ic = surfvar.nbot(p):surfvar.ntop(p)
      zl = surfvar.zw(p,ic-1);
      zu = surfvar.zw(p,ic);

      surfvar.dpai(p,ic) = 0;

      % Numerical integration between zl and zu using 100 sublayers

      num_int = 100;
      dz_int = (zu - zl) / num_int;
      for ic_int = 1:num_int

         if (ic_int == 1)
            z_int = zl + 0.5 * dz_int;
         else
            z_int = z_int + dz_int;
         end

         % beta distribution probability density function

         zrel = min(z_int/surfvar.hc(p), 1);
         beta_pdf = (zrel^(pbeta-1) * (1 - zrel)^(qbeta-1)) / beta(pbeta,qbeta);

         % Plant area density (m2/m3)

         pad = (surfvar.pai(p) / surfvar.hc(p)) * beta_pdf;

         % Plant area index (m2/m2)

         surfvar.dpai(p,ic) = surfvar.dpai(p,ic) + pad * dz_int;

      end
   end

   % Check to make sure sum of numerical integration matches canopy plant area index

   pai_sum = 0;
   for ic = surfvar.nbot(p):surfvar.ntop(p)
      pai_sum = pai_sum + surfvar.dpai(p,ic);
   end

   if (abs(pai_sum - surfvar.pai(p)) > 1e-06)
      error('plant area index error')
   end

   % Set layers with small plant area index to zero

   pai_miss = 0;
   for ic = surfvar.nbot(p):surfvar.ntop(p)
      if (surfvar.dpai(p,ic) < 0.01)
         pai_miss = pai_miss + surfvar.dpai(p,ic);
         surfvar.dpai(p,ic) = 0;
      end
   end

   % Distribute the missing plant area across vegetation layers
   % in proportion to the plant area profile

   if (pai_miss > 0)
      pai_old = pai_sum;
      pai_new = pai_old - pai_miss;
      for ic = surfvar.nbot(p):surfvar.ntop(p)
         surfvar.dpai(p,ic) = surfvar.dpai(p,ic) + pai_miss * (surfvar.dpai(p,ic) / pai_new);
      end
   end

   % Zero out non-canopy layers

   for ic = surfvar.ntop(p)+1:surfvar.nlev(p)
      surfvar.dpai(p,ic) = 0;
   end

end

% --- Canopy layer variables

for p = 1:npts

   % Wet and dry fraction of each layer

   for ic = surfvar.nbot(p):surfvar.ntop(p)
      surfvar.fwet(p,ic) = 0;
      surfvar.fdry(p,ic) = 0.8218390792391702;
   end

   % Sunlit and shaded fraction of each layer

   Kb = 1.762817445019839;
   for ic = surfvar.ntop(p): -1: surfvar.nbot(p)
      if (ic == surfvar.ntop(p))
         sumpai = 0.5 * surfvar.dpai(p,ic);
      else
         sumpai = sumpai + 0.5 * (surfvar.dpai(p,ic+1) + surfvar.dpai(p,ic));
      end
      surfvar.fracsun(p,ic) = exp(-Kb * sumpai);
      surfvar.fracsha(p,ic) = 1 - surfvar.fracsun(p,ic);
   end

   % Leaf conductances of each layer - specify boundary layer and stomatal conductances rather than calculate

   for ic = surfvar.nbot(p):surfvar.ntop(p)
      leafvar.gbh(p,ic,leafvar.isun) = 2.268731551029694;
      leafvar.gbh(p,ic,leafvar.isha) = 2.268731551029694;
      leafvar.gbv(p,ic,leafvar.isun) = 2.496430918408511;
      leafvar.gbv(p,ic,leafvar.isha) = 2.4964309184085112;
   end

   ic =  2; leafvar.gs(p,ic,leafvar.isun) = 0.1055192299591170; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  3; leafvar.gs(p,ic,leafvar.isun) = 0.1055242801518130; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  4; leafvar.gs(p,ic,leafvar.isun) = 0.1055442917843868; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  5; leafvar.gs(p,ic,leafvar.isun) = 0.1055941138135227; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  6; leafvar.gs(p,ic,leafvar.isun) = 0.1056908793262697; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  7; leafvar.gs(p,ic,leafvar.isun) = 0.1058528110177235; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  8; leafvar.gs(p,ic,leafvar.isun) = 0.1060982785799807; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  9; leafvar.gs(p,ic,leafvar.isun) = 0.1064449769725067; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 10; leafvar.gs(p,ic,leafvar.isun) = 0.1069092429994712; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 11; leafvar.gs(p,ic,leafvar.isun) = 0.1075055939313528; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 12; leafvar.gs(p,ic,leafvar.isun) = 0.1082466006575611; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 13; leafvar.gs(p,ic,leafvar.isun) = 0.1091431941575199; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 14; leafvar.gs(p,ic,leafvar.isun) = 0.1102054367680129; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 15; leafvar.gs(p,ic,leafvar.isun) = 0.1114436673549660; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 16; leafvar.gs(p,ic,leafvar.isun) = 0.1128697739318455; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 17; leafvar.gs(p,ic,leafvar.isun) = 0.1139609997662955; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 18; leafvar.gs(p,ic,leafvar.isun) = 0.1171854769922459; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 19; leafvar.gs(p,ic,leafvar.isun) = 0.1187962807892496; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 20; leafvar.gs(p,ic,leafvar.isun) = 0.1208179425962061; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 21; leafvar.gs(p,ic,leafvar.isun) = 0.1228222278933264; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 22; leafvar.gs(p,ic,leafvar.isun) = 0.1265187201079511; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 23; leafvar.gs(p,ic,leafvar.isun) = 0.1302060211116664; leafvar.gs(p,ic,leafvar.isha) = 4.2012509631176821E-003;
   ic = 24; leafvar.gs(p,ic,leafvar.isun) = 0.1324769753414125; leafvar.gs(p,ic,leafvar.isha) = 5.2367461490856618E-003;
   ic = 25; leafvar.gs(p,ic,leafvar.isun) = 0.1368698950422385; leafvar.gs(p,ic,leafvar.isha) = 5.5301527204908814E-003;
   ic = 26; leafvar.gs(p,ic,leafvar.isun) = 0.1409781702869481; leafvar.gs(p,ic,leafvar.isha) = 9.3004358317120180E-003;
   ic = 27; leafvar.gs(p,ic,leafvar.isun) = 0.1453771563594110; leafvar.gs(p,ic,leafvar.isha) = 9.4127133298050371E-003;
   ic = 28; leafvar.gs(p,ic,leafvar.isun) = 0.1500692855259890; leafvar.gs(p,ic,leafvar.isha) = 1.2673855474610703E-002;
   ic = 29; leafvar.gs(p,ic,leafvar.isun) = 0.1550625267179841; leafvar.gs(p,ic,leafvar.isha) = 1.7001098609145948E-002;
   ic = 30; leafvar.gs(p,ic,leafvar.isun) = 0.1612573684054735; leafvar.gs(p,ic,leafvar.isha) = 2.3169937120171465E-002;
   ic = 31; leafvar.gs(p,ic,leafvar.isun) = 0.1670070755851650; leafvar.gs(p,ic,leafvar.isha) = 2.7922863298797198E-002;
   ic = 32; leafvar.gs(p,ic,leafvar.isun) = 0.1729085723882166; leafvar.gs(p,ic,leafvar.isha) = 3.7604357650450893E-002;
   ic = 33; leafvar.gs(p,ic,leafvar.isun) = 0.1789636092244607; leafvar.gs(p,ic,leafvar.isha) = 4.6973864571307165E-002;
   ic = 34; leafvar.gs(p,ic,leafvar.isun) = 0.1851674300515763; leafvar.gs(p,ic,leafvar.isha) = 5.9133058475084523E-002;
   ic = 35; leafvar.gs(p,ic,leafvar.isun) = 0.1934936004678613; leafvar.gs(p,ic,leafvar.isha) = 7.2069695598315137E-002;
   ic = 36; leafvar.gs(p,ic,leafvar.isun) = 0.1998800427890550; leafvar.gs(p,ic,leafvar.isha) = 8.7715286232307205E-002;
   ic = 37; leafvar.gs(p,ic,leafvar.isun) = 0.2062193012225244; leafvar.gs(p,ic,leafvar.isha) = 0.1061384302864682;
   ic = 38; leafvar.gs(p,ic,leafvar.isun) = 0.2125360681852630; leafvar.gs(p,ic,leafvar.isha) = 0.1228400116453672;
   ic = 39; leafvar.gs(p,ic,leafvar.isun) = 0.2173608304404607; leafvar.gs(p,ic,leafvar.isha) = 0.1417637263927490;
   ic = 40; leafvar.gs(p,ic,leafvar.isun) = 0.2229022164242044; leafvar.gs(p,ic,leafvar.isha) = 0.1584830781836984;
   ic = 41; leafvar.gs(p,ic,leafvar.isun) = 0.2272705787429249; leafvar.gs(p,ic,leafvar.isha) = 0.1712625620594075;
   ic = 42; leafvar.gs(p,ic,leafvar.isun) = 0.2303710181620154; leafvar.gs(p,ic,leafvar.isha) = 0.1801177896104632;
   ic = 43; leafvar.gs(p,ic,leafvar.isun) = 0.2315642727778051; leafvar.gs(p,ic,leafvar.isha) = 0.1844525262388218;

   % Net radiation of each leaf - specify Rn rather than calculate

%  for ic = surfvar.nbot(p):surfvar.ntop(p)
%     leafvar.rnleaf(p,ic,leafvar.isun) = ...;
%     leafvar.rnleaf(p,ic,leafvar.isha) = ...;
%  end

   ic =  2; leafvar.rnleaf(p,ic,leafvar.isun) = 274.1550972520949; leafvar.rnleaf(p,ic,leafvar.isha) = 135.5805463019297;
   ic =  3; leafvar.rnleaf(p,ic,leafvar.isun) = 152.6466432111793; leafvar.rnleaf(p,ic,leafvar.isha) = 14.07209226101406;
   ic =  4; leafvar.rnleaf(p,ic,leafvar.isun) = 143.3106465111127; leafvar.rnleaf(p,ic,leafvar.isha) = 4.736095560947426;
   ic =  5; leafvar.rnleaf(p,ic,leafvar.isun) = 141.1349601292557; leafvar.rnleaf(p,ic,leafvar.isha) = 2.560409179090552;
   ic =  6; leafvar.rnleaf(p,ic,leafvar.isun) = 140.3515782682666; leafvar.rnleaf(p,ic,leafvar.isha) = 1.777027318101450;
   ic =  7; leafvar.rnleaf(p,ic,leafvar.isun) = 139.9989618149844; leafvar.rnleaf(p,ic,leafvar.isha) = 1.424410864819143;
   ic =  8; leafvar.rnleaf(p,ic,leafvar.isun) = 139.8193779719216; leafvar.rnleaf(p,ic,leafvar.isha) = 1.244827021756379;
   ic =  9; leafvar.rnleaf(p,ic,leafvar.isun) = 139.7226856943809; leafvar.rnleaf(p,ic,leafvar.isha) = 1.148134744215676;
   ic = 10; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6715516492443; leafvar.rnleaf(p,ic,leafvar.isha) = 1.097000699079142;
   ic = 11; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6486777617586; leafvar.rnleaf(p,ic,leafvar.isha) = 1.074126811593461;
   ic = 12; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6454611443424; leafvar.rnleaf(p,ic,leafvar.isha) = 1.070910194177239;
   ic = 13; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6574886918659; leafvar.rnleaf(p,ic,leafvar.isha) = 1.082937741700723;
   ic = 14; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6825518740815; leafvar.rnleaf(p,ic,leafvar.isha) = 1.108000923916256;
   ic = 15; leafvar.rnleaf(p,ic,leafvar.isun) = 139.7197003420071; leafvar.rnleaf(p,ic,leafvar.isha) = 1.145149391841863;
   ic = 16; leafvar.rnleaf(p,ic,leafvar.isun) = 139.7687630478708; leafvar.rnleaf(p,ic,leafvar.isha) = 1.194212097705659;
   ic = 17; leafvar.rnleaf(p,ic,leafvar.isun) = 139.8300945195577; leafvar.rnleaf(p,ic,leafvar.isha) = 1.255543569392491;
   ic = 18; leafvar.rnleaf(p,ic,leafvar.isun) = 139.9044342012941; leafvar.rnleaf(p,ic,leafvar.isha) = 1.329883251128853;
   ic = 19; leafvar.rnleaf(p,ic,leafvar.isun) = 139.9928218561823; leafvar.rnleaf(p,ic,leafvar.isha) = 1.418270906017080;
   ic = 20; leafvar.rnleaf(p,ic,leafvar.isun) = 140.0965359142808; leafvar.rnleaf(p,ic,leafvar.isha) = 1.521984964115638;
   ic = 21; leafvar.rnleaf(p,ic,leafvar.isun) = 140.2170312214647; leafvar.rnleaf(p,ic,leafvar.isha) = 1.642480271299498;
   ic = 22; leafvar.rnleaf(p,ic,leafvar.isun) = 140.3558546400128; leafvar.rnleaf(p,ic,leafvar.isha) = 1.781303689847561;
   ic = 23; leafvar.rnleaf(p,ic,leafvar.isun) = 140.5145139618931; leafvar.rnleaf(p,ic,leafvar.isha) = 1.939963011727875;
   ic = 24; leafvar.rnleaf(p,ic,leafvar.isun) = 140.6942683355502; leafvar.rnleaf(p,ic,leafvar.isha) = 2.119717385385012;
   ic = 25; leafvar.rnleaf(p,ic,leafvar.isun) = 140.8957966907126; leafvar.rnleaf(p,ic,leafvar.isha) = 2.321245740547398;
   ic = 26; leafvar.rnleaf(p,ic,leafvar.isun) = 141.1186839776834; leafvar.rnleaf(p,ic,leafvar.isha) = 2.544133027518268;
   ic = 27; leafvar.rnleaf(p,ic,leafvar.isun) = 141.3606433537155; leafvar.rnleaf(p,ic,leafvar.isha) = 2.786092403550295;
   ic = 28; leafvar.rnleaf(p,ic,leafvar.isun) = 141.6163675209725; leafvar.rnleaf(p,ic,leafvar.isha) = 3.041816570807255;
   ic = 29; leafvar.rnleaf(p,ic,leafvar.isun) = 141.8758806459934; leafvar.rnleaf(p,ic,leafvar.isha) = 3.301329695828233;
   ic = 30; leafvar.rnleaf(p,ic,leafvar.isun) = 142.1222598264425; leafvar.rnleaf(p,ic,leafvar.isha) = 3.547708876277293;
   ic = 31; leafvar.rnleaf(p,ic,leafvar.isun) = 142.3286461575954; leafvar.rnleaf(p,ic,leafvar.isha) = 3.754095207430179;
   ic = 32; leafvar.rnleaf(p,ic,leafvar.isun) = 142.4546334499792; leafvar.rnleaf(p,ic,leafvar.isha) = 3.880082499813955;
   ic = 33; leafvar.rnleaf(p,ic,leafvar.isun) = 142.4425091001131; leafvar.rnleaf(p,ic,leafvar.isha) = 3.867958149947935;
   ic = 34; leafvar.rnleaf(p,ic,leafvar.isun) = 142.2145595429649; leafvar.rnleaf(p,ic,leafvar.isha) = 3.640008592799720;
   ic = 35; leafvar.rnleaf(p,ic,leafvar.isun) = 141.6738568074410; leafvar.rnleaf(p,ic,leafvar.isha) = 3.099305857275830;
   ic = 36; leafvar.rnleaf(p,ic,leafvar.isun) = 140.7125813068553; leafvar.rnleaf(p,ic,leafvar.isha) = 2.138030356690095;
   ic = 37; leafvar.rnleaf(p,ic,leafvar.isun) = 139.2336249812494; leafvar.rnleaf(p,ic,leafvar.isha) = 0.6590740310842307;
   ic = 38; leafvar.rnleaf(p,ic,leafvar.isun) = 137.1921092835705; leafvar.rnleaf(p,ic,leafvar.isha) = -1.382441666594701;
   ic = 39; leafvar.rnleaf(p,ic,leafvar.isun) = 134.6629482058528; leafvar.rnleaf(p,ic,leafvar.isha) = -3.911602744312447;
   ic = 40; leafvar.rnleaf(p,ic,leafvar.isun) = 131.9427836668230; leafvar.rnleaf(p,ic,leafvar.isha) = -6.631767283342214;
   ic = 41; leafvar.rnleaf(p,ic,leafvar.isun) = 129.7331041990678; leafvar.rnleaf(p,ic,leafvar.isha) = -8.841446751097436;
   ic = 42; leafvar.rnleaf(p,ic,leafvar.isun) = 129.8107940035027; leafvar.rnleaf(p,ic,leafvar.isha) = -8.763756946662468;
   ic = 43; leafvar.rnleaf(p,ic,leafvar.isun) = 143.7593131399208; leafvar.rnleaf(p,ic,leafvar.isha) = 5.184762189755574;

   % Leaf heat capacity

   for ic = surfvar.nbot(p):surfvar.ntop(p)
      leafvar.cpleaf(p,ic) = 744.5333333333334;
   end

end

% --- Soil variables

for p = 1:npts
   soilvar.tk(p) = 1.261326601469150;          % Soil thermal conductivity (W/m/K)
   soilvar.dz(p) = 7.1006354171935350E-003;    % Soil layer depth (m)
   soilvar.tsoi(p) = 294.8492736816406;        % Soil temperature (K)
   soilvar.resis(p) = 3361.509423807650;       % Soil evaporative resistance (s/m)
   soilvar.rhg(p) = 0.9984057411945876;        % Relative humidity of airspace at soil surface (fraction)
   fluxvar.rnsoi(p) = 1.867608703739659;       % Net radiation at ground (W/m2)
end

% --- Atmospheric forcing at a reference height

for p = 1:npts

   forcvar.pref(p) = 98620;               % Atmospheric pressure (Pa)
   forcvar.uref(p) = 5.169;               % Wind speed at reference height (m/s)
   forcvar.tref(p) = 295.9349938964844;   % Air temperature at reference height (K)

   % Water vapor (mol/mol) using constant relative humidity

   relhum = 53.871;
   [esat, desat] = satvap (forcvar.tref(p)-physcon.tfrz);
   eref = (relhum / 100) * esat;
   forcvar.qref(p) = eref / forcvar.pref(p);

   % Specific humidity (kg/kg)

   qref_kg_kg = physcon.mmh2o / physcon.mmdry * eref / (forcvar.pref(p) - (1 - physcon.mmh2o/physcon.mmdry) * eref);

   % Potential temperature and virtual potential temperature (K)

   forcvar.thref(p) = forcvar.tref(p) + 0.0098 * forcvar.zref(p);
   forcvar.thvref(p) = forcvar.thref(p) * (1 + 0.61 * qref_kg_kg);

   % Molar density (mol/m3)

   forcvar.rhomol(p) = forcvar.pref(p) / (physcon.rgas * forcvar.tref(p));

   % Air density (kg/m3)

   rhoair = forcvar.rhomol(p) * physcon.mmdry * (1 - (1 - physcon.mmh2o/physcon.mmdry) * eref / forcvar.pref(p));

   % Molecular mass of air (kg/mol)

   forcvar.mmair(p) = rhoair / forcvar.rhomol(p);

   % Specific heat of air at constant pressure (J/mol/K)

   forcvar.cpair(p) = physcon.cpd * (1 + (physcon.cpw/physcon.cpd - 1) * qref_kg_kg) * forcvar.mmair(p);

end

% --- Time step (s)

dt = 5 * 60;

% --- Initialize profiles

for p = 1:npts
   for ic = surfvar.nsoi(p):surfvar.nlev(p)
      fluxvar.tair_old(p,ic) = forcvar.tref(p);
      fluxvar.qair_old(p,ic) = forcvar.qref(p);
      fluxvar.tveg_old(p,ic,leafvar.isun) = forcvar.tref(p);
      fluxvar.tveg_old(p,ic,leafvar.isha) = forcvar.tref(p);
   end
  fluxvar.taf(p) = fluxvar.tair_old(p,surfvar.ntop(p));
  fluxvar.qaf(p) = fluxvar.qair_old(p,surfvar.ntop(p));
end

% Initialize profile output variables

for p = 1:npts

   for ic = surfvar.nsoi(p):surfvar.nlev(p)
      fluxvar.wind(p,ic) = 0;
      fluxvar.tair(p,ic) = 0;
      fluxvar.qair(p,ic) = 0;
      fluxvar.tveg(p,ic,leafvar.isun) = 0;
      fluxvar.tveg(p,ic,leafvar.isha) = 0;
      fluxvar.ga_prof(p,ic) = 0;
   end

   for ic = surfvar.nbot(p):surfvar.ntop(p)
      fluxvar.shveg(p,ic) = 0;
      fluxvar.etveg(p,ic) = 0;
      fluxvar.stveg(p,ic) = 0;
   end

   for ic = surfvar.nbot(p):surfvar.nlev(p)
      fluxvar.shair(p,ic) = 0;
      fluxvar.etair(p,ic) = 0;
      fluxvar.stair(p,ic) = 0;
   end

end

% --- Calculate fluxes and scalar profiles

for p = 1:npts
   surfvar.p = p;
   [fluxvar] = canopy_turbulence (dt, physcon, forcvar, surfvar, leafvar, soilvar, fluxvar);
end

% --- Write profile output data

p = 1;
for ic = surfvar.nsoi(p):surfvar.nlev(p)
   fprintf('%6.0f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',ic-1,surfvar.zs(p,ic),surfvar.dpai(p,ic), ...
   fluxvar.wind(p,ic),fluxvar.tair(p,ic),fluxvar.qair(p,ic)*forcvar.pref(p), ...
   fluxvar.tveg(p,ic,leafvar.isun),fluxvar.tveg(p,ic,leafvar.isha))
end

% --- Make graph profile output data

p = 1;
z = surfvar.zs(p,:) / surfvar.hc(p);
x = fluxvar.tair(p,:) - 273.15;

plot(x,z)
title('Profiles')
xlabel('Air temperature (^oC)')
ylabel('Relative height')
