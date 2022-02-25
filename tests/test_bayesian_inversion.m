function [current_vals, diff_tolerances] = test_bayesian_inversion()
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
    [results, elapsed_time] = call_bayesian_inversion();

    % store the results, including required fields (elapsed time, the
    % array_comparisons cell array)
    current_vals = struct();
    current_vals.elapsed_time = elapsed_time;

    % the array_comparisons cell array is used to access nested fields in the
    % structure. Each row refers to the nested order of access for the
    % current_vals structure
    current_vals.arrays_to_compare= { ...
                                      {'results'; 'Vpo'} ...
                                      {'results'; 'Xk1_temp'} ...
                                      {'results'; 'Xk1_phi'} ...
                                      {'results'; 'Xk1_g'} ...
                                    };
    current_vals.results = results;

end


function [results_to_save, end_time] = call_bayesian_inversion()

    % tests are run from top-level, so this is relative to top.
    addpath(genpath('./functions'))
    addpath('./GLAD25')


    % load in the observations, limited by spatial extent
    spatial_sampling.elim = -10;
    spatial_sampling.wlim = -80.;
    spatial_sampling.nlim = 86.;
    spatial_sampling.slim = 56.;
    spatial_sampling.xy_res = 8;  % spatial resolution [decimal degrees; multiple of 0.5]
    spatial_sampling.z_start_index = 4; % remove top 3 slices (10, 20, 30 km)
    spatial_sampling.z_end_index = 40;
    spatial_sampling.z_res = 10; % depth sampling factor
    observations = load_data(spatial_sampling);
    % load in the vbr predictions for selected method
    qmethod = 'xfit_premelt';     % Choose from xfit_premelt, xfit_mxw, eburgers_psp, andrade_psp
    vbr_predictions = load_vbr_box('sweep_box.mat', qmethod, observations);

    % run the inversion
    % set the starting models for melt fraction and grain size
    bayesian_settings.phi0 = 0.01; % melt fraction
    bayesian_settings.g0 = log10(1000.); % 1 mm (1000 microns)
    bayesian_settings.std_T = 50.; % degree C
    bayesian_settings.std_phi = 0.0025; % melt fraction
    bayesian_settings.std_g = 0.2; % grain size in log units
    bayesian_settings.lscale = 200.; % distance scale, km
    bayesian_settings.output_dir = 0;

    tic();
    results = bayesian_inversion(bayesian_settings, observations, vbr_predictions);
    end_time = toc();

    results_to_save = struct();
    fields_to_save = {'Vpo'; 'Xk1_temp'; 'Xk1_phi'; 'Xk1_g'};
    for ifi = 1:numel(fields_to_save)
        fld = fields_to_save{ifi};
        results_to_save.(fld) = results.(fld);
    end


end
