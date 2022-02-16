% a lower resolution for quick calculation
addpath('./functions')
addpath('./inv_functions')
addpath('./GLAD25')

% load in the observations, limited by spatial extent
spatial_sampling.elim = -10;
spatial_sampling.wlim = -80.;
spatial_sampling.nlim = 86.;
spatial_sampling.slim = 56.;
spatial_sampling.xy_res = 4;  % spatial resolution [decimal degrees; multiple of 0.5]
spatial_sampling.z_start_index = 4; % remove top 3 slices (10, 20, 30 km)
spatial_sampling.z_end_index = 40;
spatial_sampling.z_res = 10; % depth sampling factor

observations = load_data(spatial_sampling);

% load in the vbr data for selected method
qmethod = 'xfit_premelt';
vbr_predictions = load_vbr_box('sweep_box.mat', qmethod, observations);

% run the inversion
% To Do: move adjustable parameters out of bayesian_inversion up to this level
bayesian_inversion(observations, vbr_predictions);


% spatial_sampling.elim = -10.;     % eastern limit [decimal degrees; multiple of 0.5]
% spatial_sampling.wlim = -80.;     % western limit [decimal degrees; multiple of 0.5]
% spatial_sampling.nlim = 86.;      % nothern limit [decimal degrees; multiple of 0.5]
% spatial_sampling.slim = 56.;      % southern limit [decimal degrees; multiple of 0.5]
% spatial_sampling.xy_res = 2;  % spatial resolution [decimal degrees; multiple of 0.5]
% spatial_sampling.z_start_index = 4; % index of the top depth [between 1 and 40].
% spatial_sampling.z_end_index = 40;  % index of the bottom depth [between 1 and 40].
% spatial_sampling.z_res = 1;         % depth sampling factor [each depth slice is 10 km].
% observations = load_data(spatial_sampling);
