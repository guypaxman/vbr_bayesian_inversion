function sweep = generate_parameter_sweep_HL(sweep_params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sweep = generate_parameter_sweep(sweep_params);
%
% Calculates the mean Vs and Q at a range of depths for a combination of
% VBR input parameters (temperature, melt fraction, and grain size) across
% the values specified in sweep_params.  This will be used as a
% look-up-table to find the best fitting parameter combination.
%
% Calculated Vs and Q is dependent on the values of the state variables
% assumed, as well as the frequency range assumed for the observed data
% and the calibration assumed for anelastic behaviour.
%
% Parameters:
% -----------
%        sweep_params       structure with the following required fields
%               T               vector of temperature values [deg C]
%               phi             vector of melt fractions [vol fraction]
%               gs              vector of grain sizes [micrometres]
%               per_bw_max      maximum period (min. freq.) considered [s]
%               per_bw_min      minimum period (max. freq.) considered [s]
%               z               vector of depths [m] %%%% Added by HL
%
% Output:
% -------
%        sweep              structure with the following fields
%               NO: z               vector of depths [m] --- removed by HL!
%               VBR             structure of fixed values for the VBR
%                               input, including information on assumed
%                               pressure, stress, water content, density
%               Box             structure, described separately below
%               All the fields on the input structure, sweep_params
%
%       sweep.Box          (numel(sweep_params.T) x numel(sweep_params.phi)
%                           x numel(sweep_params.gs)) structure.  Each 
%                          element contains a field for each of the 
%                          anelastic methods in given in
%                          VBR.in.anelastic.methods_list
%
%       sweep.Box.[anelastic method name]
%                           structure with the following fields
%               meanVs          vector of calculated Vs (mean within the
%                               given frequency range) [km/s]
%               meanQinv        vector of calculated attentuation (mean
%                               within the given frequency range) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% construct state variable fields
%%z = linspace(0,300,150)*1e3;  % change z into input choice
z= sweep_params.z'; %% switched by HL to allow user to choose z array
VBR.in.z = z;
VBR.in.SV.P_GPa = z * 3300 * 9.8 /1e9;
VBR.in.SV.sig_MPa = 0.1*ones(size(z));
VBR.in.SV.Ch2o = zeros(size(z)); % in PPM!
VBR.in.SV.rho = 3300 * ones(size(z)); % [Pa]
VBR.in.SV.chi = ones(size(z));
VBR.in.SV.f = logspace(-2.2,-1.3,10);
VBR.in.SV.T_K = (1450 + 273) * ones(size(z));
VBR.in.SV.phi = zeros(size(z));
VBR.in.SV.dg_um = 1000 * ones(size(z));
solidus_C = SoLiquidus(VBR.in.SV.P_GPa*1e9, zeros(size(z)),zeros(size(z)),...
    'hirschmann');
VBR.in.SV.Tsolidus_K = solidus_C.Tsol + 273;
% Note, the parameters that we are sweeping through will be 
% overwritten in calculate_sweep()!

% write method list (these are the things to calculate)
% Use all available methods except xfit_premelt
elastic = feval(fetchParamFunction('elastic'), '');
VBR.in.elastic.methods_list = elastic.possible_methods; 
viscous = feval(fetchParamFunction('viscous'), '');
VBR.in.viscous.methods_list = viscous.possible_methods;
anelastic = feval(fetchParamFunction('anelastic'), '');
VBR.in.anelastic.methods_list = anelastic.possible_methods; 

% load in settings that you might want to overwrite (optional)
%  (each will be called internally if you don't call them here)
% VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
% VBR.in.elastic.poro_Takei=Params_Elastic('poro_Takei'); % unrelaxed poro-elasticity
% VBR.in.viscous.HK2003 = Params_Viscous('HK2003'); % viscous parameters
% VBR.in.viscous.LH2012 = Params_Viscous('LH2012'); % viscous parameters
% VBR.in.anelastic.eBurgers = Params_Anelastic('eBurgers'); % anelasticity
% VBR.in.anelastic.AndradePsP = Params_Anelastic('AndradePsP'); % anelasticity

% Generate parameter sweep and calculate VBR at each combination
sweepBox = calculate_sweep(VBR, sweep_params);

sweep = sweep_params;
sweep.z = VBR.in.z;  %% HL
sweep.Box = sweepBox;
sweep.VBR = VBR;
sweep.P_GPa = VBR.in.SV.P_GPa;
sweep.cH2O = VBR.in.SV.Ch2o;
sweep.state_names = {'T', 'phi', 'gs'};




end


function [sweepBox] = calculate_sweep(VBR_init, sweep_params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sweepBox = calculate_sweep(sweep_params);
%
% Calls extract_meanVs_Q() to calculate mean Vs and Q values given the 
% fixed parameters given in VBR_init and across all combinations of the 
% parameters given in sweep_params (T, phi, gs) for the frequency range
% given in sweep_params.
%
% Parameters:
% -----------
%       sweep_params       structure with the following required fields
%               T               vector of temperature values [deg C]
%               phi             vector of melt fractions [vol (?) fraction]
%               gs              vector of grain sizes [micrometres]
%               per_bw_max      maximum period (min. freq.) considered [s]
%               per_bw_min      minimum period (max. freq.) considered [s]
%               z               array of z, added by HL
%
%       VBR                 structure of fixed values for the VBR input
%                           including information on assumed pressure, 
%                           stress, water content, density
%
% Output:
% -------
%       sweepBox          (numel(sweep_params.T) x numel(sweep_params.phi)
%                           x numel(sweep_params.gs)) structure.  Each 
%                          element contains a field for each of the 
%                          anelastic methods in given in
%                          VBR.in.anelastic.methods_list
%       sweep.Box.[anelastic method name]
%                           structure with the following fields
%               meanVs          vector of calculated Vs (mean within the
%                               given frequency range) [km/s]
%               meanQinv        vector of calculated attentuation (mean
%                               within the given frequency range) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define values for masking
freq_range = sort([1./sweep_params.per_bw_max 1./sweep_params.per_bw_min]);
min_freq = freq_range(1);
max_freq = freq_range(2);

n_T   = numel(sweep_params.T);
n_phi    = numel(sweep_params.phi);
n_gs     = numel(sweep_params.gs);
n_z = numel(sweep_params.z);
N_TOT = n_T*n_phi*n_gs;

if isfield(sweep_params, 'verbose') == 0
    sweep_params.verbose = 1;
end


i_state=1;
for i_T = n_T:-1:1  % run backwards so structure is preallocated
    for i_phi = n_phi:-1:1
        for i_gs = n_gs:-1:1
            if sweep_params.verbose
                disp([num2str(i_state),' of ',num2str(N_TOT)])
            end
            % copy initial VBR structure and overwrite input state
            % variables for asthenospheric depth
            VBR = VBR_init;
            VBR.in.SV.phi = sweep_params.phi(i_phi) ...
                * ones(size(VBR.in.z));
            VBR.in.SV.dg_um = sweep_params.gs(i_gs) ...
                * ones(size(VBR.in.z));
            VBR.in.SV.T_K = sweep_params.T(i_T) + 273 ...
                * ones(size(VBR.in.z));
            
            % run VBR for asthenospheric depths (and current phi, gs, T)
            VBR = VBR_spine(VBR);
            anelastic_methods = fieldnames(VBR.out.anelastic);

            for i_an=1:length(anelastic_methods)
                Vs_zf = mean(VBR.out.anelastic.(anelastic_methods{i_an}).V,2);
                Q_zf = mean(VBR.out.anelastic.(anelastic_methods{i_an}).Q,2);
                calculated_vals.(anelastic_methods{i_an}).Vs = Vs_zf;
                calculated_vals.(anelastic_methods{i_an}).Q = Q_zf;

            end
            sweepBox(i_T, i_phi, i_gs) = calculated_vals;
            

            i_state=i_state+1;
        end
    end
end

sweep = sweep_params;
sweep.z = VBR.in.z;
sweep.Box = sweepBox; 
sweep.VBR = VBR;
sweep.P_GPa = VBR.in.SV.P_GPa;
sweep.cH2O = VBR.in.SV.Ch2o;
sweep.state_names = {'T', 'phi', 'gs'};


end

function calculated_vals = extract_meanVs_Q(VBR, min_freq, max_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sweepBox = extract_meanVs_Q(VBR, min_freq, max_freq);
%
% Calculates the mean Vs and Q at a range of depths for the initial VBR
% values given
%
% Parameters:
% -----------
%       VBR                 structure of fixed values for the VBR input
%                           including information on assumed pressure, 
%                           stress, water content, density, temperature,
%                           grain size, and melt fraction.
%
% Output:
% -------
%       calculated_vals    Structure contains a field for each of the 
%                          anelastic methods in given in
%                          VBR.in.anelastic.methods_list
%       calculated_vals.[anelastic method name]
%                           structure with the following fields
%               meanVs          vector of calculated Vs (mean within the
%                               given frequency range) [km/s]
%               meanQinv        vector of calculated attentuation (mean
%                               within the given frequency range) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anelastic_methods = fieldnames(VBR.out.anelastic);
freq = VBR.in.SV.f;
nn_pts = 1000; % number of points to interpolate frequency axis

for i_an = 1:length(anelastic_methods)
    Q_zf = VBR.out.anelastic.(anelastic_methods{i_an}).Q;
    Vs_zf = VBR.out.anelastic.(anelastic_methods{i_an}).V/1000; % m/s to km/s
    
    % 2D interplations of Qinv(z,f) and Vs(z,f) over frequency only
    [Vs_zf, ~, ~] = interp_FreqZ(Vs_zf, freq, nn_pts,...
        VBR.in.z, length(VBR.in.z));
    [Q_zf, freq_interp, ~] = interp_FreqZ(Q_zf, freq, nn_pts,...
        VBR.in.z, length(VBR.in.z));
    
    % Mask Vs and Q in frequency and depth
    f_mask = ((freq_interp >= min_freq) & (freq_interp <= max_freq));
    Vs_zf_mask = Vs_zf(:, f_mask);
    Q_zf_mask = Q_zf(:, f_mask);
    
    % Calculate mean Vs and Q of mask
    meanVs = mean(Vs_zf_mask, 2);
    meanQ = mean(Q_zf_mask, 2);
    calculated_vals.(anelastic_methods{i_an}).meanVs = meanVs;
    calculated_vals.(anelastic_methods{i_an}).meanQ = meanQ;
    
end

end