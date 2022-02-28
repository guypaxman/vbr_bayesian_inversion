function [current_vals, diff_tolerances] = test_make_vm()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % retrieve results from make_Vm()
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
    [Vm, elapsed_time] = call_make_vm();

    % store the results, including required fields (elapsed time, the
    % array_comparisons cell array)
    current_vals = struct();
    current_vals.elapsed_time = elapsed_time;

    % the array_comparisons cell array is used to access nested fields in the
    % structure. Each row refers to the nested order of access for the
    % current_vals structure
    current_vals.arrays_to_compare= { ...
                                      {'output'; 'Vm'} ...
                                    };
    current_vals.output.Vm = Vm;  % ({'output'; 'Vm'} will be used to access output.Vm
    
end 


function [Vm, end_time] = call_make_vm()

    % tests are run from top-level, so this is relative to top.
    addpath(genpath('./functions'))

    lons = linspace(-10, -80,4);
    lats = linspace(56, 86, 5);
    zs = linspace(50, 400, 6);

    npts = length(lons) * length(lats) * length(zs);

    std_T = ones(npts)*50.;
    std_phi = ones(npts)*0.0025;
    std_g = ones(npts)*0.2; % log units
    lscale = 200.; % km

    pars = {'T'; 'phi'; 'g'};
    npar = numel(pars); 
    nmod = npts * npar;

    tic();
    Vm = make_Vm(std_T, std_phi, std_g, lscale, lats, lons, zs, npts, nmod);
    end_time = toc();
end 
