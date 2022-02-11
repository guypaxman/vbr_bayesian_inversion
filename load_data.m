function observations = load_data(spatial_sampling)
    %%%
    % must return a structure with the following fields 
    %   .vs(lat, lon, depth)  3D shear wave velocity in km/s
    %   .Q(lat, lon, depth)   3D quality factor
    %   .vs_err(lat, lon, depth)  error for each vs measurement
    %   .Q_err(lat, lon, depth)   error for each Q measurent 
    %
    %   .lat   1D latitdue array, degrees 
    %   .lon   1D longitude array, degrees
    %   .zs    1D depth array, in km (???)
    %%%
    
    % The first set is the tomography model (vs and Q)
    GLAD25_Vs = load('Vs_Model.mat');
    QL6_Q = load('QL6_Model.mat');
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
    zs = GLAD25_Vs.Vs_Model.Depth(idz); % remove top 3 slices (10, 20, 30 km)
    
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
