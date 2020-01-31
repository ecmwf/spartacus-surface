% Matlab script to plot the Russell Square urban case
base = 'russell_square';

scenarios = {'1stream','2stream','4stream'};
leg = {'1 stream','2 stream','4 stream'}

for iscen = 1:length(scenarios)
  d{iscen} = loadnc([base '_' scenarios{iscen} '_out.nc']);
  d{iscen}.net_lost = d{iscen}.top_flux_dn_sw ...
		      - d{iscen}.top_flux_net_sw ...
		      + sum(d{iscen}.veg_absorption_sw,1) ...
		      + sum(d{iscen}.roof_net_sw,1) ...
		      + sum(d{iscen}.wall_net_sw,1);
end

in = loadnc([base '.nc']);

mu0 = in.cos_solar_zenith_angle;
sza = acosd(mu0);

stys = {'-','--','-.'};

clf
for iscen = 1:length(scenarios)
  plot(sza,d{iscen}.ground_flux_dn_sw,['b' stys{iscen}])
  hold on
  plot(sza,sum(d{iscen}.veg_absorption_sw,1), ['g' stys{iscen}]);
  plot(sza,d{iscen}.top_flux_dn_sw - d{iscen}.top_flux_net_sw, ['r' stys{iscen}])
  plot(sza,d{iscen}.net_lost,['k' stys{iscen}])
  plot(sza,sum(d{iscen}.wall_net_sw,1), ['m' stys{iscen}])
  plot(sza,sum(d{iscen}.roof_net_sw,1), ['c' stys{iscen}])
end
