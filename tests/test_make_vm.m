function test_make_vm(store_answer)
    addpath('./inv_functions')
    
    answer_dir = './tests/testdata';
    answer_file = [answer_dir, '/test_make_vm.mat'];

    % call the functions 

    lons = linspace(-10, -80, 7);
    lats = linspace(56, 86, 6);
    zs = linspace(50, 400, 10);

    npts = length(lons) * length(lats) * length(zs);

    std_T = ones(npts)*50.;
    std_phi = ones(npts)*0.0025;
    std_g = ones(npts)*0.2; % log units
    % std_g = ones(npts)*2000; % linear units (microns)
    lscale = 200.; % km

    pars = {'T'; 'phi'; 'g'};
    npar = numel(pars); 
    nmod = npts * npar;
    disp(nmod)
    disp(npar)
    disp(npts)
    tic();
    Vm = make_Vm(std_T, std_phi, std_g, lscale, lats, lons, zs, npts, nmod);
    end_time = toc();
    
    time_diff_tol = -0.05;  % negative is slower
    diff_val_tol = 1e-14;
    
    if store_answer==1
        expected.Vm = Vm;
        expected.execution_time = end_time;
        save(answer_file, '-struct', 'expected')
    else 
        % reload the stored result 
        expected = load(answer_file);
        
        % calculate differences
        diffvals = abs(Vm - expected.Vm)./Vm;
        timediff = (expected.execution_time - end_time) / expected.execution_time;
                        
        if max(diffvals(:)) < diff_val_tol
            tpct = num2str(timediff * 100);
            disp(['test_make_vm passed with a time difference of ', tpct, ' percent'])
        else            
            answer_file_bad = [answer_dir, '/test_make_vm_bad.mat'];
            bad_answer.Vm = Vm;
            bad_answer.execution_time = end_time;
            save(answer_file_bad, '-struct', 'bad_answer')
            error(['test_make_vm failed, check ', answer_file_bad]);
        end
    end
    
end 
