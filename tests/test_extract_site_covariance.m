function [current_vals, diff_tolerances] = test_extract_site_covariance()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % retrieve results from bayesian_inversion()
    %
    % Returns
    % -------
    % current_vals : struct
    %    the results of calling the test function, must contain the following
    %    fields:
    %
    %    .elapsed_time : scalar
    %        the elapsed time of the function call in seconds
    %    .arrays_to_compare  : cell array
    %        a cell array with the nested field names pointing to the arrays to
    %        run comparisons on. For example, with the following cell array:
    %
    %        arrays_to_compare = { ...
    %                              {'output'; 'array_1'} ...
    %                              {'output'; 'array_2'} ...
    %                              {'output'; 'out_1'; 'array_3'} ...
    %                            };
    %
    %        the test comparison will expect to the following arrays :
    %
    %               current_vals.output.array_1
    %               current_vals.output.array_2
    %               current_vals.out_1.array_3
    %
    %    .{array_fields}  : arrays
    %         any of the array fields specified by .arrays_to_compare should
    %         also be fields themselves.
    %
    % diff_tolerances : struct
    %     the allowed tolerances for all array comparisons for this test, must
    %     have the following fields:
    %
    %     .max_abs_diff : scalar
    %         the max allowed difference in absolute value:
    %             abs(expected - current)
    %     .max_frac_diff : scalar
    %         the max allowed fractional difference:
    %             abs(expected - current)/expected
    %         (values with zero in the denominator will be ignored).
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set the tolerances for this test
    diff_tolerances = struct();
    diff_tolerances.max_abs_diff = 1e-8;
    diff_tolerances.max_frac_diff = 1e-10;

    % go execute the code!
    [results, elapsed_time] = call_extract_site_covariance();

    % store the results, including required fields (elapsed time, the
    % array_comparisons cell array)
    current_vals = struct();
    current_vals.elapsed_time = elapsed_time;

    % the array_comparisons cell array is used to access nested fields in the
    % structure. Each row refers to the nested order of access for the
    % current_vals structure
    current_vals.arrays_to_compare= { ...
                                      {'results'; 'vpo_exact'} ...
                                      {'results'; 'vpo_closest'} ...
                                    };
    current_vals.results = results;

end


function [results_to_save, end_time] = call_extract_site_covariance()

    addpath(genpath('./functions'))
    nz = 5;
    nlon = 7;
    nlat = 6;
    npts = nz * nlon * nlat;
    npar = 3; % T, phi, g = 3 parameters
    vpo_size = npts * npar;

    vpo_test = linspace(0, 1, vpo_size);
    [vpog1, vpog2] =  ndgrid(vpo_test, vpo_test);

    bayes_result.Vpo = vpog1;
    bayes_result.npar = npar;
    observations = struct();
    observations.lats = linspace(-90, 90, nlat);
    observations.lons = linspace(-180, 180, nlon);
    observations.zs = linspace(0, 500, npts);
    observations.nz = nz;
    observations.nlon = nlon;
    observations.nlat = nlat;

    tic();

    latitude = observations.lats(2);
    longitude = observations.lons(2);
    vpo_exact = extract_site_covariance(bayes_result, observations, latitude, longitude, 1);

    latitude = 42;
    longitude = -110;
    vpo_closest = extract_site_covariance(bayes_result, observations, latitude, longitude, 0);
    end_time = toc();

    results_to_save = struct();
    results_to_save.vpo_exact = vpo_exact;
    results_to_save.vpo_closest = vpo_closest;

end
