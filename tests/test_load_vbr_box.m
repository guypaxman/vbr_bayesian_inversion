function [current_vals, diff_tolerances] = test_load_vbr_box()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % retrieve results from load_vbr_box()
    %
    % note that this requires an existing sweep on disk, which will be
    % generated on the first run and saved in ./tests/testdata/. Subsequent
    % runs will re-load that sweep, so only the `load_vbr_box` is timed here.
    %
    % see test_make_vm for description of the output structures here.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set the tolerances for this test
    diff_tolerances = struct();
    diff_tolerances.max_abs_diff = 1e-8;
    diff_tolerances.max_frac_diff = 1e-10;

    % go execute the code!
    box_file = './tests/testdata/test_box.mat';
    initialize_sweep(box_file);
    [output, elapsed_time] = call_load_vbr(box_file);

    % store the results, including required fields (elapsed time, the
    % array_comparisons cell array)
    current_vals = struct();
    current_vals.elapsed_time = elapsed_time;

    % the array_comparisons cell array is used to access nested fields in the
    % structure. Each row refers to the nested order of access for the
    % current_vals structure
    current_vals.arrays_to_compare= { ...
                                      {'output'; 'Q_vbr'} ...
                                      {'output'; 'vs_vbr'} ...
                                      {'output'; 'g'} ...
                                      {'output'; 'T'} ...
                                      {'output'; 'phi'} ...
                                    };
    current_vals.output = output;  % ({'output'; 'Vm'} will be used to access output.Vm

end


function [output, end_time] = call_load_vbr(box_file)
    addpath(genpath('./functions'))
    qmethod = "xfit_premelt";
    observations.nz = 6; % match the sweep_params from initialize_sweep()
    tic()
    output = load_vbr_box(box_file, qmethod, observations);
    end_time = toc();
end


function initialize_sweep(box_file)
    % generate the sweep box if it does not exist locally. This is tested
    % separately from testing generate_parameter_sweep_HL
    if exist(box_file) ~= 2
        addpath(genpath('./functions'))
        initialize_vbr()

        sweep_params.per_bw_max = 0.8; % max period (s)
        sweep_params.per_bw_min = 1.2; % min period (s)
        sweep_params.z = linspace(50,400,6)*1000.; % depth (m)
        sweep_params.T = linspace(1200, 1800, 4); %[degrees C]
        sweep_params.phi = (0.0:0.01:0.02); % melt fraction
        sweep_params.gs = linspace(0.001,0.01,4)*1e6; % grain size [micrometres]
        sweep_params.verbose = 0;
        sweep = generate_parameter_sweep_HL(sweep_params);

        save(box_file, 'sweep')
    end
end
