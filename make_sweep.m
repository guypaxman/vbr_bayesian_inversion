close all; clc

% % Determine where your m-file's folder is.
% folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath('./functions')
initialize_vbr();
%% THIS SCRIPT WILL GENERATE LOOK UP TABLES

fileout = 'sweep_box.mat';
sweep_params.per_bw_max = 0.8; % max period (s)
sweep_params.per_bw_min = 1.2; % min period (s)
sweep_params.z = (40:10:400)*1000.; % depth (m)


% OPTION 1: LONG RUN (~HOURS)
% things to choice:
sweep_params.T = 1100:20:1800; %[degrees C]
sweep_params.phi = (0.0:0.0025:0.05); % melt fraction

% discretize grain size:
% linear:
% sweep_params.gs = linspace(0.001,0.01,10)*1e6; % grain size [micrometres]
% sweep_params.gs_params = struct('type','linear');
% log:
gsmin = 0.0001*1e6; gsmax = 0.01*1e6; gsref = 0.001*1e6; 
sweep_params.gs = gsref * exp(linspace(log(gsmin/gsref),log(gsmax/gsref),21));
sweep_params.gs_params = struct('type','log','gsmin',gsmin,'gsmax',gsmax,'gsref',gsref);

% OPTION 2: QUICK TEST (10-ish min)
%sweep_params.T = 1200:50:1800; %[degrees C]
%sweep_params.phi = (0.0:0.005:0.05); % melt fraction
% sweep_params.gs = linspace(0.001,0.01,10)*1e6; % grain size [micrometres]


% % OPTION 3: QUICK QUICK TEST
% sweep_params.T = 1200:250:1800; %[degrees C]
% sweep_params.phi = (0.0:0.01:0.05); % melt fraction
% sweep_params.gs = linspace(0.001,0.01,4)*1e6; % grain size [micrometres]

sweep = generate_parameter_sweep_HL(sweep_params);

save(fileout,'sweep');

clear sweep_params