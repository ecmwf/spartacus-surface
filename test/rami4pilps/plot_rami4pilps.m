% Matlab script to plot the RAMI4PILPS case with different
% SPARTACUS-Surface settings

base = 'vis-snw-0.3'; thetitle = 'Visible, snow, vegetation cover 0.3';
%base = 'nir-med-0.3'; thetitle = 'Near-IR, bare-soil, vegetation cover 0.3';

scenarios={[base '-1-1'], [base '-1-2'], [base '-2-1'], [base '-2-2']};
leg = {'1 region 2 stream','1 region 4 stream','2 region 2 stream','2 region 4 stream'};

for iscen = 1:length(scenarios)
  d{iscen}=loadnc(['rami4pilps_' scenarios{iscen} '_out.nc']);

  d{iscen}.net_lost = d{iscen}.top_flux_dn_sw-d{iscen}.top_flux_net_sw ...
		      +d{iscen}.veg_absorption_sw(2,:)' ...
		      +d{iscen}.ground_flux_net_sw;
end
in = loadnc('rami4pilps.nc');

mu0 = in.cos_solar_zenith_angle;
sza = acosd(mu0);
stys = {'-','--','-.',':'};

clf
plot(-1,-1,'b',-1,-1,'r',-1,-1,'g');
hold on

for iscen = 1:length(scenarios)
  plot(-1,-1,['k' stys{iscen}]);
  hold on
end

for iscen = 1:length(scenarios)
  plot(sza,d{iscen}.top_flux_dn_sw-d{iscen}.top_flux_net_sw,['r' stys{iscen}]);
  hold on
  plot(sza,d{iscen}.ground_flux_dn_sw,['b' stys{iscen}]);
  plot(sza,d{iscen}.veg_absorption_sw(2,:),['g' stys{iscen}]);
  %  plot(sza,d{iscen}.net_lost, ['k' stys{iscen}]);
end
legend([{'Reflectance','Transmittance','Absorptance'} leg])
axis([0 90 0 1]);
xlabel('Solar zenith angle (\circ)');
ylabel('Normalized flux')
title(thetitle)
