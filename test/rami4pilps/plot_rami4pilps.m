base = 'vis-snw-0.3'
scenarios={[base '-1-1'], [base '-1-2'], [base '-2-1'],[base '-2-2']};
%,...
%	   [base '-1-4'], [base '-2-4']};
%scenarios={'vis-med-0.1','nir-snw-0.5'};

for iscen = 1:length(scenarios)
  d{iscen}=loadnc(['rami4pilps_' scenarios{iscen} '_out.nc']);
end
in = loadnc('rami4pilps.nc');
mu0 = in.cos_solar_zenith_angle;
sza = acosd(mu0);
stys = {'-','--','-.',':'};

clf
for iscen = 1:length(scenarios)
  plot(sza,d{iscen}.ground_flux_dn_sw.*mu0,['b' stys{iscen}]);
  hold on
  plot(sza,d{iscen}.veg_absorption_sw(2,:),['g' stys{iscen}]);
end
