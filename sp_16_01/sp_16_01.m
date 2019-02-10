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

   % Find the lowest vegetation layer

   for ic = surfvar.ntop(p): -1: surfvar.nbot(p)
      if (surfvar.dpai(p,ic) > 0)
         ic_bot = ic;
      end
   end
   surfvar.nbot(p) = ic_bot;

   % Zero out non-vegetation layers

   ic = surfvar.nsoi(p);
   surfvar.dpai(p,ic) = 0;

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
      leafvar.gbv(p,ic,leafvar.isha) = 2.496430918408511;
   end

   ic =  7; leafvar.gs(p,ic,leafvar.isun) = 0.1056193510550169; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  8; leafvar.gs(p,ic,leafvar.isun) = 0.1058669704208841; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic =  9; leafvar.gs(p,ic,leafvar.isun) = 0.1062166035088956; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 10; leafvar.gs(p,ic,leafvar.isun) = 0.1066846074875817; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 11; leafvar.gs(p,ic,leafvar.isun) = 0.1072854387286280; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 12; leafvar.gs(p,ic,leafvar.isun) = 0.1080315168674592; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 13; leafvar.gs(p,ic,leafvar.isun) = 0.1089335362366439; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 14; leafvar.gs(p,ic,leafvar.isun) = 0.1100012607812562; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 15; leafvar.gs(p,ic,leafvar.isun) = 0.1112447128077408; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 16; leafvar.gs(p,ic,leafvar.isun) = 0.1126755044648808; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 17; leafvar.gs(p,ic,leafvar.isun) = 0.1138467165585616; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 18; leafvar.gs(p,ic,leafvar.isun) = 0.1170524695200598; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 19; leafvar.gs(p,ic,leafvar.isun) = 0.1186451281076514; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 20; leafvar.gs(p,ic,leafvar.isun) = 0.1206859738130298; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 21; leafvar.gs(p,ic,leafvar.isun) = 0.1228219389652392; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 22; leafvar.gs(p,ic,leafvar.isun) = 0.1263235652964973; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 23; leafvar.gs(p,ic,leafvar.isun) = 0.1300019677357508; leafvar.gs(p,ic,leafvar.isha) = 2.0000000000000000E-003;
   ic = 24; leafvar.gs(p,ic,leafvar.isun) = 0.1322680545506565; leafvar.gs(p,ic,leafvar.isha) = 5.2146013029975334E-003;
   ic = 25; leafvar.gs(p,ic,leafvar.isun) = 0.1367071935229807; leafvar.gs(p,ic,leafvar.isha) = 5.5227688387169205E-003;
   ic = 26; leafvar.gs(p,ic,leafvar.isun) = 0.1408216759258680; leafvar.gs(p,ic,leafvar.isha) = 9.2945439124555301E-003;
   ic = 27; leafvar.gs(p,ic,leafvar.isun) = 0.1452273039039047; leafvar.gs(p,ic,leafvar.isha) = 9.4101275089266447E-003;
   ic = 28; leafvar.gs(p,ic,leafvar.isun) = 0.1499262843535941; leafvar.gs(p,ic,leafvar.isha) = 1.2582674218550544E-002;
   ic = 29; leafvar.gs(p,ic,leafvar.isun) = 0.1549264640058029; leafvar.gs(p,ic,leafvar.isha) = 1.6999874421743270E-002;
   ic = 30; leafvar.gs(p,ic,leafvar.isun) = 0.1611234013632947; leafvar.gs(p,ic,leafvar.isha) = 2.3036435105984941E-002;
   ic = 31; leafvar.gs(p,ic,leafvar.isun) = 0.1668845999057947; leafvar.gs(p,ic,leafvar.isha) = 2.7903866816023401E-002;
   ic = 32; leafvar.gs(p,ic,leafvar.isun) = 0.1727971327085968; leafvar.gs(p,ic,leafvar.isha) = 3.7385308959493969E-002;
   ic = 33; leafvar.gs(p,ic,leafvar.isun) = 0.1788628079180081; leafvar.gs(p,ic,leafvar.isha) = 4.6808450662473224E-002;
   ic = 34; leafvar.gs(p,ic,leafvar.isun) = 0.1850771375553107; leafvar.gs(p,ic,leafvar.isha) = 5.9036977283335762E-002;
   ic = 35; leafvar.gs(p,ic,leafvar.isun) = 0.1934140277837149; leafvar.gs(p,ic,leafvar.isha) = 7.1890808689035093E-002;
   ic = 36; leafvar.gs(p,ic,leafvar.isun) = 0.1998116684650200; leafvar.gs(p,ic,leafvar.isha) = 8.7547774703355410E-002;
   ic = 37; leafvar.gs(p,ic,leafvar.isun) = 0.2061626747018590; leafvar.gs(p,ic,leafvar.isha) = 0.1059444058487105;
   ic = 38; leafvar.gs(p,ic,leafvar.isun) = 0.2124795008223110; leafvar.gs(p,ic,leafvar.isha) = 0.1228398700721039;
   ic = 39; leafvar.gs(p,ic,leafvar.isun) = 0.2173241738995193; leafvar.gs(p,ic,leafvar.isha) = 0.1416660859387607;
   ic = 40; leafvar.gs(p,ic,leafvar.isun) = 0.2228796106202699; leafvar.gs(p,ic,leafvar.isha) = 0.1584170776550386;
   ic = 41; leafvar.gs(p,ic,leafvar.isun) = 0.2272584280787935; leafvar.gs(p,ic,leafvar.isha) = 0.1712280540285039;
   ic = 42; leafvar.gs(p,ic,leafvar.isun) = 0.2303662043528620; leafvar.gs(p,ic,leafvar.isha) = 0.1801048624092090;
   ic = 43; leafvar.gs(p,ic,leafvar.isun) = 0.2315636153119537; leafvar.gs(p,ic,leafvar.isha) = 0.1844507421254655;

   % Net radiation of each leaf - specify Rn rather than calculate

%  for ic = surfvar.nbot(p):surfvar.ntop(p)
%     leafvar.rnleaf(p,ic,leafvar.isun) = ...;
%     leafvar.rnleaf(p,ic,leafvar.isha) = ...;
%  end

   ic =  7; leafvar.rnleaf(p,ic,leafvar.isun) = 139.9869857739781; leafvar.rnleaf(p,ic,leafvar.isha) =  1.411488333307743;
   ic =  8; leafvar.rnleaf(p,ic,leafvar.isun) = 139.8100113537029; leafvar.rnleaf(p,ic,leafvar.isha) =  1.234513913032590;
   ic =  9; leafvar.rnleaf(p,ic,leafvar.isun) = 139.7147998761629; leafvar.rnleaf(p,ic,leafvar.isha) =  1.139302435492522;
   ic = 10; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6645467566822; leafvar.rnleaf(p,ic,leafvar.isha) =  1.089049316011852;
   ic = 11; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6422035725484; leafvar.rnleaf(p,ic,leafvar.isha) =  1.066706131878055;
   ic = 12; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6392966303582; leafvar.rnleaf(p,ic,leafvar.isha) =  1.063799189687813;
   ic = 13; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6514847604817; leafvar.rnleaf(p,ic,leafvar.isha) =  1.075987319811279;
   ic = 14; leafvar.rnleaf(p,ic,leafvar.isun) = 139.6766021357984; leafvar.rnleaf(p,ic,leafvar.isha) =  1.101104695127997;
   ic = 15; leafvar.rnleaf(p,ic,leafvar.isun) = 139.7137254254163; leafvar.rnleaf(p,ic,leafvar.isha) =  1.138227984745972;
   ic = 16; leafvar.rnleaf(p,ic,leafvar.isun) = 139.7627019640728; leafvar.rnleaf(p,ic,leafvar.isha) =  1.187204523402388;
   ic = 17; leafvar.rnleaf(p,ic,leafvar.isun) = 139.8238999626867; leafvar.rnleaf(p,ic,leafvar.isha) =  1.248402522016342;
   ic = 18; leafvar.rnleaf(p,ic,leafvar.isun) = 139.8980702313243; leafvar.rnleaf(p,ic,leafvar.isha) =  1.322572790653995;
   ic = 19; leafvar.rnleaf(p,ic,leafvar.isun) = 139.9862631887909; leafvar.rnleaf(p,ic,leafvar.isha) =  1.410765748120540;
   ic = 20; leafvar.rnleaf(p,ic,leafvar.isun) = 140.0897684653183; leafvar.rnleaf(p,ic,leafvar.isha) =  1.514271024647946;
   ic = 21; leafvar.rnleaf(p,ic,leafvar.isun) = 140.2100538053315; leafvar.rnleaf(p,ic,leafvar.isha) =  1.634556364661143;
   ic = 22; leafvar.rnleaf(p,ic,leafvar.isun) = 140.3486818847138; leafvar.rnleaf(p,ic,leafvar.isha) =  1.773184444043382;
   ic = 23; leafvar.rnleaf(p,ic,leafvar.isun) = 140.5071806149416; leafvar.rnleaf(p,ic,leafvar.isha) =  1.931683174271291;
   ic = 24; leafvar.rnleaf(p,ic,leafvar.isun) = 140.6868352048059; leafvar.rnleaf(p,ic,leafvar.isha) =  2.111337764135555;
   ic = 25; leafvar.rnleaf(p,ic,leafvar.isun) = 140.8883584829672; leafvar.rnleaf(p,ic,leafvar.isha) =  2.312861042296822;
   ic = 26; leafvar.rnleaf(p,ic,leafvar.isun) = 141.1113792315132; leafvar.rnleaf(p,ic,leafvar.isha) =  2.535881790842842;
   ic = 27; leafvar.rnleaf(p,ic,leafvar.isun) = 141.3536664423189; leafvar.rnleaf(p,ic,leafvar.isha) =  2.778169001648555;
   ic = 28; leafvar.rnleaf(p,ic,leafvar.isun) = 141.6099822011559; leafvar.rnleaf(p,ic,leafvar.isha) =  3.034484760485563;
   ic = 29; leafvar.rnleaf(p,ic,leafvar.isun) = 141.8704336551236; leafvar.rnleaf(p,ic,leafvar.isha) =  3.294936214453243;
   ic = 30; leafvar.rnleaf(p,ic,leafvar.isun) = 142.1181913093900; leafvar.rnleaf(p,ic,leafvar.isha) =  3.542693868719610;
   ic = 31; leafvar.rnleaf(p,ic,leafvar.isun) = 142.3264909566734; leafvar.rnleaf(p,ic,leafvar.isha) =  3.750993516003037;
   ic = 32; leafvar.rnleaf(p,ic,leafvar.isun) = 142.4550034158019; leafvar.rnleaf(p,ic,leafvar.isha) =  3.879505975131512;
   ic = 33; leafvar.rnleaf(p,ic,leafvar.isun) = 142.4460421185886; leafvar.rnleaf(p,ic,leafvar.isha) =  3.870544677918208;
   ic = 34; leafvar.rnleaf(p,ic,leafvar.isun) = 142.2218178601452; leafvar.rnleaf(p,ic,leafvar.isha) =  3.646320419474797;
   ic = 35; leafvar.rnleaf(p,ic,leafvar.isun) = 141.6851596824207; leafvar.rnleaf(p,ic,leafvar.isha) =  3.109662241750371;
   ic = 36; leafvar.rnleaf(p,ic,leafvar.isun) = 140.7277716843982; leafvar.rnleaf(p,ic,leafvar.isha) =  2.152274243727867;
   ic = 37; leafvar.rnleaf(p,ic,leafvar.isun) = 139.2518034108234; leafvar.rnleaf(p,ic,leafvar.isha) =  0.676305970153017;
   ic = 38; leafvar.rnleaf(p,ic,leafvar.isun) = 137.2114197261891; leafvar.rnleaf(p,ic,leafvar.isha) = -1.364077714481233;
   ic = 39; leafvar.rnleaf(p,ic,leafvar.isun) = 134.6805463548995; leafvar.rnleaf(p,ic,leafvar.isha) = -3.894951085770870;
   ic = 40; leafvar.rnleaf(p,ic,leafvar.isun) = 131.9550915266485; leafvar.rnleaf(p,ic,leafvar.isha) = -6.620405914021802;
   ic = 41; leafvar.rnleaf(p,ic,leafvar.isun) = 129.7361873094630; leafvar.rnleaf(p,ic,leafvar.isha) = -8.839310131207370;
   ic = 42; leafvar.rnleaf(p,ic,leafvar.isun) = 129.7993862020948; leafvar.rnleaf(p,ic,leafvar.isha) = -8.776111238575538;
   ic = 43; leafvar.rnleaf(p,ic,leafvar.isun) = 143.7045065806239; leafvar.rnleaf(p,ic,leafvar.isha) =  5.129009139953610;

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
   fluxvar.rnsoi(p) = 1.896127799819662;       % Net radiation at ground (W/m2)
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
