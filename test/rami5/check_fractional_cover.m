% This Matlab script should be run after typing "make overhead": it
% uses the direct flux at the surface
% (ground_spectral_flux_dn_direct_sw) to diagnose the fractional scene
% coverage, and compare it with the values provided on the RAMI-V web
% site.

cases = {'HET07_JPS_SUM',...
	 'HET08_OPS_WIN',...
	 'HET09_JBS_SUM',...
	 'HET14_WCO_UND',...
	 'HET15_JBS_WIN'};

% "True" values from the RAMI-V website
fractional_scene_coverage = [0.406 0.1248 0.5044 0.392 0.2510];

% Load SPARTACUS overhead results
for ic = 1:length(cases)
  d=loadnc(['out_nc/rami5_' cases{ic} '-00_out.nc']);
  fsc_spartacus(ic) = 1.0 - d.ground_spectral_flux_dn_direct_sw(end);
end

% Compare to "truth"
for ic = 1:length(cases)
  disp([cases{ic} ' FSCtrue=' num2str(fractional_scene_coverage(ic)) ...
	     ', FSCspartacus=' num2str(fsc_spartacus(ic)) ' (' ...
	     num2str(100.*(fsc_spartacus(ic)-fractional_scene_coverage(ic)) ...
		     ./ fractional_scene_coverage(ic)) '%)']);
end
