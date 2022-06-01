function observations = load_data(spatial_sampling, data_dir)
    %%%
    %
    % This function loads in Q and Vs observations and associated data (see
    % Return description for required fields). As written, it relies on
    % finding the GLAD25-related model files, 'Vs_Model.mat' and 'QL6_Model.mat'.
    % These files are not included in the repository (see the top level
    % readme for information on how to get them).
    %
    % Parameters
    % ----------
    % spatial_sampling : struct
    %     a structure with information on sampling the observations. Fields
    %     include:
    %       .elim : float
    %           eastern limit [decimal degrees; multiple of 0.5]
    %       .wlim : float
    %           western limit [decimal degrees; multiple of 0.5]
    %       .nlim : float
    %           nothern limit [decimal degrees; multiple of 0.5]
    %       .slim : float
    %           southern limit [decimal degrees; multiple of 0.5]
    %       .xy_res : int
    %           spatial resolution [decimal degrees; multiple of 0.5]
    %       .z_start_index : int
    %           index of the top depth [between 1 and 40].
    %       .z_end_index : int
    %           index of the bottom depth [between 1 and 40].
    %       .z_res : int
    %           depth sampling factor [each depth slice is 10 km].
    %
    % data_dir : string
    %    the path where the GLAD25 data files reside
    %
    % Returns
    % -------
    % structure with the following fields:
    %
    %   .vs(lat, lon, depth)  3D shear wave velocity in km/s
    %   .Q(lat, lon, depth)   3D quality factor
    %   .vs_err(lat, lon, depth)  error for each vs measurement
    %   .Q_err(lat, lon, depth)   error for each Q measurement 
    %
    %   .lats   1D latitude array, degrees 
    %   .lons   1D longitude array, degrees
    %   .zs    1D depth array, in km
    %
    %   .nz, .nlat, .nlon integer values with the length of each lat, lon, zs array
    %   .npts  integer values with total points (nz * nlat * nlon)
    %%%
    
    % The first set is the tomography model (vs and Q)
    % TO DO: move these files up a level or consider a plugin-framework
    GLAD25_Vs = load([data_dir, '/Vs_Model.mat']);
    QL6_Q = load([data_dir, '/QL6_Model.mat']);
    lats = GLAD25_Vs.Vs_Model.Latitude;
    lons = GLAD25_Vs.Vs_Model.Longitude;

    % Set the vs and Q models as "vs" and "Q"
    % trim model to correct area and resolution
    idxlat = find((lats>=spatial_sampling.slim)&(lats<=spatial_sampling.nlim));
    idxlon = find((lons>=spatial_sampling.wlim)&(lons<=spatial_sampling.elim));
    idxlat = idxlat(1:spatial_sampling.xy_res*2:end);
    idxlon = idxlon(1:spatial_sampling.xy_res*2:end);
    lats = lats(idxlat);
    lons = lons(idxlon);
    
    z0 = spatial_sampling.z_start_index;
    z1 = spatial_sampling.z_end_index;
    z_res = spatial_sampling.z_res;
    idz = z0:z_res:z1;
    vs = GLAD25_Vs.Vs_Model.Vs(idxlat,idxlon,idz); % km/s
    Q = QL6_Q.QL6_Model.Q(idxlat,idxlon,idz);
    zs = GLAD25_Vs.Vs_Model.Depth(idz); % km
    
    % store in the observations structure
    observations.nlat = length(lats);
    observations.nlon = length(lons);
    observations.nz = length(zs);
    observations.npts = observations.nlat*observations.nlon*observations.nz;
    observations.vs = vs;
    observations.Q = Q;
    observations.Q_err = ones(size(Q))*10.;
    observations.vs_err = ones(size(vs))*0.05;
    observations.lats = lats;
    observations.lons = lons;
    observations.zs = zs;

end
