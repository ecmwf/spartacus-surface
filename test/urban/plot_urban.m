% Matlab script to plot the Russell Square urban case
base = 'russell_square';

scenarios = {'1stream','2stream','4stream'};
leg = {'1 stream','2 stream','4 stream'}


in = loadnc([base '_sza.nc']);

mu0 = in.cos_solar_zenith_angle;
sza = acosd(mu0);

if ~exist('d','var')
for iscen = 1:length(scenarios)
  d{iscen} = loadnc([base '_' scenarios{iscen} '_out.nc']);
  d{iscen}.net_lost = d{iscen}.top_flux_dn_sw ...
		      - d{iscen}.top_flux_net_sw ...
		      + sum(d{iscen}.veg_absorption_sw,1)' ...
		      + sum(d{iscen}.roof_flux_net_sw,1)'...
		      + sum(d{iscen}.wall_flux_net_sw,1)' ...
		      + d{iscen}.ground_flux_net_sw;
end
end

stys = {'-','--','-.'};

clf
set(gcf,'defaultlinelinewidth',1.5);
for iscen = 1:length(scenarios)
  plot(sza,d{iscen}.ground_flux_dn_sw,['b' stys{iscen}])
  hold on
  plot(sza,sum(d{iscen}.veg_absorption_sw,1), ['g' stys{iscen}]);
  plot(sza,d{iscen}.top_flux_dn_sw - d{iscen}.top_flux_net_sw, ['r' stys{iscen}])
  plot(sza,d{iscen}.net_lost,['k' stys{iscen}])
  plot(sza,sum(d{iscen}.wall_flux_net_sw,1), ['m' stys{iscen}])
  plot(sza,sum(d{iscen}.roof_flux_net_sw,1), ['c' stys{iscen}])
end
legend('Ground','Vegetation','Top','Net','Wall','Roof');
