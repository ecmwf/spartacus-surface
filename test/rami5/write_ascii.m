function write_ascii(scene, bands, zenith, meas, model, values, varname, wavelength)
  if nargin < 7
    wavelength = [];
  end

  nband = length(bands);

  if strcmp(zenith,'diffuse')
    illum = 'DIFFUSE';
  else
    zen = str2num(zenith);
    azim = zeros(1,91);
    azim([56 41 76 42 60 67]+1) = [153 147 155 76 45 41];
    illum = sprintf('z%02da%03d', zen, azim(zen+1));
  end

  for iband = 1:nband
    filename = ['mes/' scene '_' bands{iband} '_' illum '-' meas '_' model '.mes'];
    disp(['Writing ' filename]);
    fid = fopen(filename,'w');
    fprintf(fid,'%.6f\t%.6f\n', values(iband), -1);
    fclose(fid);
  end

  if ~isempty(wavelength) & length(values(:)) > 1
    for iw = 1:length(wavelength)
      wl{iw} = num2str(round(wavelength(iw)));
    end
    clf
    set(gcf,'paperposition',[0.5 0.5 18 12]);
    plot(values,'k')
    set(gca,'xticklabel',wl,'xtick',[1:length(wl)]);
    xlabel('Wavelength (nm)');
    h=ylabel([varname ' (' meas ')']);
    set(h,'interpreter','none');
    h=title([scene ' ' illum ' ' model]);
    set(h,'interpreter','none');
    drawnow
    filename = ['mes/' scene '_' illum '-' meas '_' model '.png'];
    print_png(filename, '150');
  end
