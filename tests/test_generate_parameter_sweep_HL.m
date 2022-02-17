function [current_vals, diff_tolerances] = test_generate_parameter_sweep_HL()

    % set the tolerances for this test
    diff_tolerances = struct();
    diff_tolerances.max_abs_diff = 1e-8;
    diff_tolerances.max_frac_diff = 1e-10;

    % go execute the code!
    [sweep, elapsed_time] = make_a_sweep();

    % store the results, including required fields (elapsed time, the
    % array_comparisons cell array)
    current_vals = struct();
    current_vals.elapsed_time = elapsed_time;

    % the array_comparisons cell array is used to access nested fields in the
    % structure. Each row refers to the nested order of access for the
    % current_vals structure
    current_vals.arrays_to_compare= { ...
                                      {'sweep'; 'gs'} ...
                                      {'sweep'; 'phi'} ...
                                      {'sweep'; 'T'} ...
                                      {'sweep'; 'P_GPa'} ...
                                    };

    % programmatically add the {'sweep'; 'box_i'; 'qmethod'; 'Vs_or_q'} entries
    i_array = numel(current_vals.arrays_to_compare) + 1;
    box_ops = {'box_1'; 'box_n'};
    method_ops = {'andrade_psp'; 'eburgers_psp'; 'xfit_premelt'; 'xfit_mxw'};
    var_ops = {'Vs'; 'Q'};

    for ibox = 1:numel(box_ops)
        for imeth = 1:numel(method_ops)
            for ivar = 1:numel(var_ops)
                new_entry = {'sweep'; box_ops{ibox}; method_ops{imeth}; var_ops{ivar}};
                current_vals.arrays_to_compare{i_array} = new_entry;
                i_array = i_array + 1;
            end
        end
    end

    current_vals.sweep = sweep;  % ({'output'; 'Vm'} will be used to access output.Vm

end


function [sweep, end_time] = make_a_sweep()

    % tests are run from top-level, so this is relative to top.
    addpath('./functions')
    initialize_vbr()

    tic();

    sweep_params.per_bw_max = 0.8; % max period (s)
    sweep_params.per_bw_min = 1.2; % min period (s)
    sweep_params.z = linspace(50,400,6)*1000.; % depth (m)
    sweep_params.T = linspace(1200, 1800, 4); %[degrees C]
    sweep_params.phi = (0.0:0.01:0.02); % melt fraction
    sweep_params.gs = linspace(0.001,0.01,4)*1e6; % grain size [micrometres]
    sweep_params.verbose = 0;
    sweep = generate_parameter_sweep_HL(sweep_params);

    end_time = toc();

    % pop up some of the Box indices to make comparisons easier
    sweep.box_1 = sweep.Box(1);
    sweep.box_n = sweep.Box(numel(sweep.Box));
end
