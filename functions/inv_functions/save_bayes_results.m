function save_bayes_results(output_dir, observations, results)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save output from bayesian inversion
    %
    % Parameters
    % ----------
    % output_dir : string
    %   the directory to save to (will try to mkdir if not found)
    % observations : struct
    %   the observations structure including lats, lons and zs arrays
    % results : struct
    %   the results structure from bayesian_inversion()
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lats = observations.lats;
    lons = observations.lons;
    zs = observations.zs;

    % also save some output
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save_max_likelihood(output_dir, lats, lons, zs, results)

    Vpo_var = results.Vpo_var;
    vpofile = [output_dir, '/Vpo.mat'];
    disp(["saving Vpo_var to ", vpofile])
    save(vpofile,'Vpo_var');

    fullresults = [output_dir, '/bayes_results.mat'];
    disp(["saving full bayesian inversion results to ", fullresults])
    save(fullresults, '-struct', 'results')
end


function save_max_likelihood(output_dir, lats, lons, depth, results)

    fields_to_save = {'Xk1_temp'; 'Xk1_phi'; 'Xk1_g'};

    % build a struct with field aliases for each field
    % each one is {'fieldname'; 'structure_name'};
    save_info.Xk1_temp = {'temperature'; 'Temp'};
    save_info.Xk1_phi = {'meltfraction'; 'Phi'};
    save_info.Xk1_g = {'grainsize'; 'GrainSize'};

    for ifield = 1:numel(fields_to_save)
        field = fields_to_save{ifield};

        % (re)initialize the struct to save
        OuterStruct = struct();
        s_name = save_info.(field){2};
        OuterStruct.(s_name) = struct();
        OuterStruct.(s_name).latitude = lats;
        OuterStruct.(s_name).longitude = lons;
        OuterStruct.(s_name).depth = depth;

        % add the field to save
        fieldnametosave = save_info.(field){1};
        OuterStruct.(s_name).(fieldnametosave) = results.(field);

        filename = [output_dir, '/', s_name, '.mat'];
        disp(['Saving max likelihood of ', fieldnametosave, ' to ', filename])
        save(filename, '-struct', 'OuterStruct', s_name)
    end
end
