%% INVERSION FUNCTIONS 4
function Vm = make_Vm(std_T,std_phi,std_g,lscale,lats,lons,zs,npts,nmod)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % make_Vm = model covariance matrix
    %
    % the final covariance matrix will have a size of (nmod, nmod) with
    % covariances for (T, phi, g) arranged in internal 3x3 subblocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get the cartesian grid at ever lat, lon zs permuation
    [X, Y, Z] = get_cartesian_grid(lats, lons, zs);

    % build the covariance matrix for each parameter
    Vm = zeros(nmod); % nmod = npts * npar (npar == 3)

    % get the columns for T, phi and g given the 3x3 subblocks
    npt_range = 1:npts;
    icols_T = (npt_range - 1) * 3 + 1;
    icols_phi = icols_T + 1;
    icols_g = icols_T + 2;

    % for every site, calculate the covariance with all others, weighted by
    % the distance between sites.
    for isite = 1:npts
        d = get_dist(X(isite), Y(isite), Z(isite), X, Y, Z);
        fac = transpose(exp(-abs(d)/lscale)(:));  % the distance weighting

        irow = (isite - 1) * 3 + 1;
        Vm(irow, icols_T) = std_T(isite) .* std_T(isite, :) .* fac;
        Vm(irow+1, icols_phi) = std_phi(isite) .* (std_phi(isite, :) .* fac);
        Vm(irow+2, icols_g) = std_g(isite) .* (std_g(isite, :) .* fac);
    end
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