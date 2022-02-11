% should exactly reproduce run_invert_xfit_premelt.m, with the caveat that the results
% end up in results/state_variables_xfit_premelt/

% load in the observations, limited by spatial extent
spatial_sampling.elim = -10;
spatial_sampling.wlim = -80.;
spatial_sampling.nlim = 86.;
spatial_sampling.slim = 56.;
spatial_sampling.xy_res = 2;  % spatial resolution [decimal degrees; multiple of 0.5]
spatial_sampling.z_start_index = 4; % remove top 3 slices (10, 20, 30 km)
spatial_sampling.z_end_index = 40;
spatial_sampling.z_res = 1; % depth sampling factor
observations = load_data(spatial_sampling);

% load in the vbr predictions for selected method
qmethod = 'xfit_premelt';
vbr_predictions = load_vbr_box('sweep_box.mat', qmethod, observations);

% run the inversion
bayesian_inversion(observations, vbr_predictions);
