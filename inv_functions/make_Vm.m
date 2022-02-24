%% INVERSION FUNCTIONS 4
function Vm = make_Vm(std_T,std_phi,std_g,lscale,lats,lons,zs,npts,nmod)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % make_Vm = model covariance matrix
    %
    % the final covariance matrix will have a size of (nmod, nmod) with
    % covariances for each parameter (T, phi, g) arranged in internal blocks.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get the cartesian grid at ever lat, lon zs permuation
    [X, Y, Z] = get_cartesian_grid(lats, lons, zs);

    % build the covariance matrix for each parameter
    co_T = zeros(npts);
    co_phi = zeros(npts);
    co_g = zeros(npts);

    for isite = 1:npts
        d = get_dist(X(isite), Y(isite), Z(isite), X, Y, Z);
        fac = transpose(exp(-abs(d)/lscale)(:));
        std_facd =  std_T(isite, :) .* fac;
        co_T(isite, :) = std_T(isite, isite) .* std_facd;
        co_phi(isite, :) = std_phi(isite, isite) .* (std_phi(isite, :) .* fac);
        co_g(isite, :) = std_g(isite, isite) .* (std_g(isite, :) .* fac);
    end

    % now assemble into expected form.
    Vm = zeros(nmod); % nmod = npts * npar (npar == 3)

    % note that if npar is ever not 3, this will fail.
    npt_range = 1:npts;
    icols = (npt_range - 1) * 3 + 1;
    irows = (npt_range - 1) * 3 + 1;
    Vm(irows, icols) = co_T;

    icols = icols + 1;
    irows = irows + 1;
    Vm(irows, icols) = co_phi;

    icols = icols + 1;
    irows = irows + 1;
    Vm(irows, icols) = co_g;

end

function [X, Y, Z] = get_cartesian_grid(lats, lons, zs)

    % make vectors of repeated lats and lons and zs
    [vlats, vlons, vzs] = ndgrid(lats, lons, zs);
    vlats = permute(vlats, [3, 2, 1]);
    vlons = permute(vlons, [3, 2, 1]);
    vzs = permute(vzs, [3, 2, 1]);

    % get the cartesian position of the full grid
    [X, Y, Z] = geodesic_to_cart(vlats, vlons, vzs);

end

function [x, y, z] = geodesic_to_cart(lat, lon, depth)
    % depth in km, lat/lon in deg
    r = 6371. - depth;
    x = r .* sind(90. - lat) .* cosd(lon);
    y = r .* sind(90. - lat) .* sind(lon);
    z = r .* cosd(90. - lat);
end


function d = get_dist(x1, y1, zz1, x2, y2, zz2)

    % get distance between the two points
    dx = x2 - x1;
    dy = y2 - y1;
    dz = zz2 - zz1;
    d = sqrt(dx.*dx+dy.*dy+dz.*dz);

end