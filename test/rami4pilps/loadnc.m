function [data, attribute, dimensions] = load_nc_struct(nc_file, names);
% load_nc_struct -- Load NetCDF variables and attributes.
%
% [data, attribute] = load_nc_struct('nc_file') loads all variables of
%   'nc_file' into structure 'data' and all attributes into structure
%   'attribute', so variable 'X' could then be accessed using 'data.X'
%   and attribute 'long_name' of 'X' could be accessed with
%   'attribute.X.long_name'.  
 
if nargin < 1, help(mfilename), return, end

result = [];
if nargout > 0, data = []; attribute = []; end

rootid = netcdf.open(nc_file, 'nowrite');
if isempty(rootid), return, end
disp(['Loading ' nc_file]);

ver_str = version;
if str2num(ver_str(1)) >= 8
  ncid = get_science_data_group(rootid);
else
  ncid = rootid;
end

if nargout > 2
  dimids = netcdf.inqDimIDs(ncid);
  for id = dimids
    [dimname,dimlen] = netcdf.inqDim(ncid,id);
    dimnames{id+1} = dimname;
    dimensions.(dimname) = dimlen;
  end
end

if nargin < 2
  names = var_names(ncid);
end

allnames = var_names(ncid);

disp(['Variables:']);
for ii = 1:length(names)
  if any(strcmp(names{ii},allnames))
    newname = names{ii};
    newname(find(newname == '-')) = '_';
    newname(find(newname == '.')) = '_';
    newname(find(newname == '@')) = '_';
    varid = netcdf.inqVarID(ncid, names{ii});

    [varname, xtype,vdimids,natts] = netcdf.inqVar(ncid, varid);
    if xtype == 12 
       continue;
    end
    [add_offset, scale_factor] = get_scaling(ncid, varid);

    if xtype == netcdf.getConstant('NC_FLOAT') | ~isempty(add_offset) | ~isempty(scale_factor)
      missing_value = double(get_missing_value(ncid, varid));
      data.(newname) = double(netcdf.getVar(ncid, varid));
    else
      missing_value = get_missing_value(ncid, varid);
      data.(newname) = netcdf.getVar(ncid, varid);
    end

    if ~isempty(missing_value)
      data.(newname)(find(data.(newname) == missing_value)) = NaN;
    end

    if ~isempty(scale_factor)
       data.(newname) = data.(newname) .* scale_factor;
    else
       scale_factor = 1.0;
    end
    if ~isempty(add_offset)
       data.(newname) = data.(newname) + add_offset;
    else
       add_offset = 0.0;
    end
    
    if nargout > 1
      % Do attributes
      for jj = 1:natts
	% Check for underscore starting attribute name
	attname = netcdf.inqAttName(ncid, varid, jj-1);
	attname_mod = attname;
	if strcmp(attname,'_FillValue')
	  attname_mod = 'FillValue_';
	elseif attname(1) == '_'
	  warning([newname '.' attname ' changed to ' names{ii} ':X' attname]);
	  attname_mod = ['X' attname];
	end
	eval(['attribute.' names{ii} '.' attname_mod ' = netcdf.getAtt(ncid, varid, attname);']);
	if nargout > 2
	  attribute.(names{ii}).dimensions = dimnames(vdimids+1);
	end
      end
    end
  
    the_size = size(data.(newname));
    if the_size(end) == 1;
       the_size = the_size(1:end-1);
    end

    long_name = '';
    try
      long_name = netcdf.getAtt(ncid, varid, 'long_name');
      long_name = clean_up_string(long_name);
      long_name = ['"' long_name '"'];
    end
    units = '';
    try
      units = netcdf.getAtt(ncid, varid, 'units');
      units = clean_up_string(units);
      units = [' (' units ')'];
    end
    
    names{ii} = newname;
    

    size_str = num2str(the_size(1));
    for ii = 2:length(the_size);
      size_str = [size_str ',' num2str(the_size(ii))];
    end

    newname = [newname ' (' size_str ')'];
    
    namefill = blanks(max(0,30-length(newname)));


    if add_offset ~= 0.0 | scale_factor ~= 1.0
       scale_str = ' scaled';
    else
       scale_str = '';
    end
    disp([namefill newname ': ' long_name units scale_str]);

  end
end

if nargout > 1
  % Do global attributes
  disp('Global attributes:');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for ii = 1:ngatts
    attname = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),ii-1);

    newattname = attname;
    if attname(1) == '_';
      newattname = [attname(2:end) '_'];
    end
    attribute.global.(newattname) = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);

% if isempty(find(attnames{ii} == '/' | attnames{ii} == '-' | attnames{ii} == '+' | attnames{ii} == ' ') > 0)
%   eval(['attr = f.' attnames{ii} '(:);']);
%   if ischar(attr)
%     attr = clean_up_string(attr);
%   end
%   eval(['attribute.global.' attnames{ii} ' = attr;']);
    namefill = blanks(max(0,14-length(attname)));
    disp([namefill newattname ': ' num2str(attribute.global.(newattname))]);
% else
%   disp(['Attribute name ' attnames{ii} ' incompatible with Matlab - not loaded.'])
% end
  end
end
netcdf.close(rootid)

function newstr = clean_up_string(oldstr)
newstr = num2str(oldstr);
if length(newstr) > 1
  if newstr(end-1) == '\' & newstr(end) == '0'
    newstr = deblank(newstr(1:end-2));
  end
end

function names = var_names(ncid)
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for ii = 0:nvars-1
    names{ii+1} = netcdf.inqVar(ncid, ii);
  end

function missing_value = get_missing_value(ncid, varid)
  missing_value = [];
  try
    missing_value = netcdf.getAtt(ncid, varid, 'missing_value');
  catch exception
    try
      missing_value = netcdf.getAtt(ncid, varid, '_FillValue');
    end
  end

function [add_offset,scale_factor] = get_scaling(ncid, varid)
  add_offset = [];
  scale_factor = [];
  try
    add_offset   = netcdf.getAtt(ncid, varid, 'add_offset');
  end
  try
    scale_factor = netcdf.getAtt(ncid, varid, 'scale_factor');
  end

function sciid = get_science_data_group(ncid)
  group_id = netcdf.inqGrps(ncid);
  sciid = ncid;
  if ~isempty(group_id)
    for ii = 1:length(group_id)
      if strcmp(netcdf.inqGrpName(group_id(ii)),'ScienceData')
	 sciid = group_id(ii);
	 return;
      end
    end
  end

