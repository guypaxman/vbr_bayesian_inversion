function [current_vals, diff_tolerances] = test_make_vm()

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
    addpath('./inv_functions')

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
