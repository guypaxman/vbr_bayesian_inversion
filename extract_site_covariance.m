function vpo_site = extract_site_covariance(bayes_result, observations, latitude, longitude, exact_only)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % extract the covariance from the bayesian inversion results for a specified
    % latitude and longitude.
    %
    % Parameters
    % ----------
    % bayes_result : struct
    %    the bayesian results structure output from bayesian_inversion()
    % observations : struct
    %    the observations structure used in bayesian_inversion()
    % latitude : float
    %    the latitude to extract at
    % longitude : float
    %    the longitude to extract at
    % exact_only: integer
    %    if 1, then will ony find exact lat/lon matches. if 0, will return results
    %    at the closest lat/lon.
    %
    % Returns
    % -------
    % vpo_site : array
    %    the covariance as a function of depth at the specified lat/lon.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lats = observations.lats;
    lons = observations.lons;

    npar = bayes_result.npar;
    nz = observations.nz;
    nlon = observations.nlon;

    if exact_only == 1
        ilat = find(lats==latitude);
        ilon = find(lons==longitude);
    else
        dlat = abs(lats - latitude);
        dlon = abs(lons - longitude);
        ilat = find(dlat == min(dlat));
        ilon = find(dlon == min(dlon));
    end

    ilat = npar*nz*nlon*(ilat-1);
    ilon = npar*nz* ilon;

    start_idx = ilat + ilon - npar * nz;
    end_idx = ilat + ilon;

    vpo_site = bayes_result.Vpo(start_idx:end_idx, :);
end
