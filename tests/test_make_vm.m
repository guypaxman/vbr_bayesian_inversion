function test_make_vm(store_answer)
    % 
    answer_file = 'test_make_vm.mat';
    bad_file = 'test_make_vm_bad.mat';  % only saved if store_answer == 0 and test fails
    diff_val_frac_tol = 1e-14;
    
    [Vm, end_time] = call_make_vm();  % go execute the code!

    time_diff_tol = -0.05;  % negative is slower
    diff_val_tol = 1e-14;

    if store_answer==1
        expected = struct();
        expected.Vm = Vm;
        expected.execution_time = end_time;
        store_answer_file(answer_file, expected)
    else
        % reload the stored result 
        expected = load_answer(answer_file);

        % calculate differences
        diffvals = abs(Vm - expected.Vm)./Vm;
        tdiff_abs = expected.execution_time - end_time;
        tdiff_frac = tdiff_abs / expected.execution_time;

        if max(diffvals(:)) < diff_val_frac_tol
            tpct = num2str(tdiff_frac * 100);
            tda = num2str(tdiff_abs);
            disp(['test_make_vm passed with dt of ', tda,' s (', tpct, ' percent)'])
        else            
            answer_file_bad = 'test_make_vm_bad.mat';
            bad_answer.Vm = Vm;
            bad_answer.execution_time = end_time;
            save(answer_file_bad, '-struct', 'bad_answer')
            error(['test_make_vm failed, check ', answer_file_bad]);
        end
    end
    
end 


function [Vm, end_time] = call_make_vm()

    addpath('./inv_functions')  % tests are run from top-level, so this is relative to top.

    lons = linspace(-10, -80,5);
    lats = linspace(56, 86, 6);
    zs = linspace(50, 400, 7);

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
