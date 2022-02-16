% should exactly reproduce run_invert_xfit_premelt.m, with the caveat that the results
% end up in results/state_variables_xfit_premelt/
addpath('./functions')
addpath('./inv_functions')
addpath('./GLAD25')

% load in the observations, limited by spatial extent
spatial_sampling.elim = -10.;     % eastern limit [decimal degrees; multiple of 0.5]
spatial_sampling.wlim = -80.;     % western limit [decimal degrees; multiple of 0.5]
spatial_sampling.nlim = 86.;      % nothern limit [decimal degrees; multiple of 0.5]
spatial_sampling.slim = 56.;      % southern limit [decimal degrees; multiple of 0.5]
spatial_sampling.xy_res = 2;  % spatial resolution [decimal degrees; multiple of 0.5]
spatial_sampling.z_start_index = 4; % index of the top depth [between 1 and 40].
spatial_sampling.z_end_index = 40;  % index of the bottom depth [between 1 and 40].
spatial_sampling.z_res = 1;         % depth sampling factor [each depth slice is 10 km].
observations = load_data(spatial_sampling);

% load in the vbr predictions for selected method
qmethod = 'xfit_premelt';     % Choose from xfit_premelt, xfit_mxw, eburgers_psp, andrade_psp
vbr_predictions = load_vbr_box('sweep_box.mat', qmethod, observations);

% run the inversion
bayesian_inversion(observations, vbr_predictions);
