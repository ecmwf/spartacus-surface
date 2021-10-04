function write_ascii_profile(scene, bands, zenith, meas, model, zmax, z, flux_up, flux_dn);
  for iband = 1:length(bands)
    if strcmp(zenith,'diffuse')
      illum = 'DIFFUSE';
    else
      zen = str2num(zenith);
      azim = zeros(1,91);
      azim([56 41 76 42 60 67]+1) = [153 147 155 76 45 41];
      illum = sprintf('z%02da%03d', zen, azim(zen+1));
    end

    zint = linspace(zmax, 0, 11);
    fup  = interp1(z, flux_up(iband,:), zint);
    fdn  = interp1(z, flux_dn(iband,:), zint);

    filename = ['mes/' scene '_' bands{iband} '_' illum '-' meas '_' model '.mes'];
    disp(['Writing ' filename]);
    fid = fopen(filename,'w');
    fprintf(fid,'%4d %4d\t%.6f\n', 11, 3, zmax/10);
    for iz = 1:11
      fprintf(fid,'%.6f\t%.6f\t%.6f\n', zint(iz), fup(iz), fdn(iz));
    end
    fclose(fid);
  end
