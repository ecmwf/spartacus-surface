% This Matlab script processes the output of SPARTACUS-Surface run on
% a single RAMI-V scene, producing ASCII files suitable for submission
% to the intercomparison in the "mes" directory.  You can call this
% script directly (optionally specifying the scene by predefining the
% "iscene" variable) or run the process_spartacus_scene.m script to
% run them all.

% Scene number in range 0 to 4
if ~exist('iscene','var')
  iscene = 1;
end

% Optionally use the RAMI-IV solar zenith angles
if ~exist('is_rami4','var')
  is_rami4 = 0;
end

% Obtain scene ID string and other information from the "iscene"
% variable
if iscene == 0
  scene_id = 'HET15_JBS_WIN'; zmax = 30.5130; solar_ids = {'diffuse', '76','56'};
  if is_rami4
    solar_ids = {'diffuse','54'};
  end
elseif iscene == 1
  scene_id = 'HET09_JBS_SUM'; zmax = 30.5130; solar_ids = {'diffuse', '56','41'};
  if is_rami4
    solar_ids = {'diffuse','37'};
  end
elseif iscene == 2
  scene_id = 'HET07_JPS_SUM'; zmax = 18.56;   solar_ids = {'diffuse', '56','41'};
  if is_rami4
    solar_ids = {'diffuse','37'};
  end
elseif iscene == 3
  scene_id = 'HET14_WCO_UND'; zmax = 4.12;    solar_ids = {'diffuse', '42','60','67'};
  if is_rami4
    solar_ids = {'diffuse','00','20','50'};
  end
elseif iscene == 4
  scene_id = 'HET08_OPS_WIN'; zmax = 15.0213; solar_ids = {'diffuse', '76','56'};
  if is_rami4
    solar_ids = {'diffuse','47'};
  end
else
  error(['iscene=' num2str(iscene) ' not recognised']);
end

% Band names
bands = {'O03','O04','O06','O08','O10','O11','O12','M08','O17','MD5','M11','MD7','M12'};

% Some outputs are only required in the photosynthetically active
% range
bands_par = bands(1:5);
iband_black = length(bands)+1; % 14th band is for black surfaces

% Filenames are tagged with the model name
model = 'spartacus';

% Make output directory
mkdir mes

% Loop over the solar configurations for this scene
for sid = 1:length(solar_ids)
  solar_id = solar_ids{sid};

  % Load the netCDF files
  in = loadnc(['scene_nc/rami5_' scene_id '_scene.nc']);
  out= loadnc(['out_nc/rami5_' scene_id '-' solar_id '_out.nc']);
  bs = loadnc(['out_nc/rami5_' scene_id '-' solar_id '_blacksoil_out.nc']);
  
  % If wavelengths are provided then write_ascii will attempt to plot
  % the results
  %wav = in.wavelength(1:13);
  wav = [];

  residual = out.ground_spectral_flux_net_sw ...
	     + sum(out.wall_spectral_flux_net_sw,2) ...
	     + sum(out.roof_spectral_flux_net_sw,2) ...
	     + sum(out.veg_spectral_absorption_sw,2) ...
	     - out.top_spectral_flux_net_sw;

  if strcmp(solar_id, 'diffuse')
    % White-sky albedo
    write_ascii(scene_id, bands, solar_id, 'bhr', model, out.top_spectral_flux_dn_sw - out.top_spectral_flux_net_sw, ...
		'White-sky albedo', wav);
  else
    % Black-sky albedo
    write_ascii(scene_id, bands, solar_id, 'dhr', model, out.top_spectral_flux_dn_sw - out.top_spectral_flux_net_sw, ...
		'Black-sky albedo', wav);
  end

  % Absorption by all vegetation (foliage and wood)
  write_ascii(scene_id, bands_par, solar_id, 'fabs_tot', model, ...
	      sum(out.veg_spectral_absorption_sw+out.wall_spectral_flux_net_sw+out.roof_spectral_flux_net_sw, 2), ...
	      'Vegetation absorption', wav);

  % Absorption by foliage only
  foliage_ratio = (ones(14,1)*in.foliage_extinction').*in.foliage_sw_ssa ...
		  ./ ((ones(14,1)*in.veg_extinction').*in.veg_sw_ssa);
  foliage_ratio(find(isnan(foliage_ratio))) = 0;
  write_ascii(scene_id, bands_par, solar_id, 'fabs_fol', model, ...
	      sum(out.veg_spectral_absorption_sw.*foliage_ratio, 2),...
	      'Foliage absorption', wav);

  % Transmission scattered one or more timesby vegetation but not by
  % soil
  if strcmp(solar_id, 'diffuse')
    write_ascii(scene_id, bands_par, solar_id, 'ftran_coco', model, ...
	       bs.ground_spectral_flux_dn_sw - bs.ground_spectral_flux_dn_sw(iband_black), ...
	       'Canopy-only collided transmission', wav);
  else
    write_ascii(scene_id, bands_par, solar_id, 'ftran_coco', model, ...
		bs.ground_spectral_flux_dn_sw - bs.ground_spectral_flux_dn_direct_sw, ...
		'Canopy-only collided transmission', wav);
  end

  % Unscattered transmission
  write_ascii(scene_id, bands_par, solar_id, 'ftran_uc', model, ...
	      bs.ground_spectral_flux_dn_sw(iband_black).*ones(size(bs.ground_spectral_flux_dn_sw)), ...
	      'Uncollided transmission', wav);

  % Total transmission
  write_ascii(scene_id, bands_par, solar_id, 'ftran_tot', model, out.ground_spectral_flux_dn_sw, ...
	      'Total transmission', wav);

  % Flux profile is average of layer top and layer base values, which
  % are slightly different at the interfaces because of the stepped
  % nature of the trunk description
  flux_up = [out.spectral_flux_up_layer_base_sw zeros(14,1)] ...
	    + [zeros(14,1) out.spectral_flux_up_layer_top_sw];
  flux_up(:,2:end-1) = 0.5.*flux_up(:,2:end-1);
  flux_dn = [out.spectral_flux_dn_layer_base_sw zeros(14,1)] ...
	    + [zeros(14,1) out.spectral_flux_dn_layer_top_sw];
  flux_dn(:,2:end-1) = 0.5.*flux_dn(:,2:end-1);

  write_ascii_profile(scene_id, bands_par, solar_id, 'ftran_tot_vprof', model, ...
		      zmax, out.height, flux_up, flux_dn);
end
