function vbr_predictions = load_vbr_box(box_file, qmethod, observations)
    %%%
    % Parameters
    % ----------
    %
    % box_file 
    %   file path to box to load, e.g., 'sweep_box.mat'
    % qmethod
    %   the qmethod to use, e.g., 'xfit_premelt'
    %
    % Returns
    % -------
    % vbr_predictions structure with following fields:
    %
    %   .vs_vbr     isotropic shear wave velocity, km/s
    %   .Q_vbr      quality factor
    %   .T          temperature in degrees celsius
    %   .phi        volumetric melt fraction 
    %   .g          log10 of grain size in micrometers 
    %   .qmethod    the anelastic method, string  
    %%%
    
    load(box_file); % e.g., 'sweep_box.mat', containing a sweep structure
    params = make_param_grid(sweep.state_names, sweep);
    
    %  Extract the VBR box
    T = unique(params.T);
    phi = unique(params.phi);
    g = log10(unique(params.gs));
    nT = length(T);
    nphi = length(phi);
    ng = length(g);
    nz = observations.nz;
    vs_vbr = zeros(nz,nT,nphi,ng);
    Q_vbr = zeros(nz,nT,nphi,ng);

    for iT=1:nT
        for iphi=1:nphi
            for ig=1:ng
                for iz=1:nz
                    vs_vbr(iz,iT,iphi,ig) = sweep.Box(iT,iphi,ig).(qmethod).Vs(iz)*1.e-3; % km/s
                    Q_vbr(iz,iT,iphi,ig) = sweep.Box(iT,iphi,ig).(qmethod).Q(iz);
                end
            end
        end
    end
    
    vbr_predictions.vs_vbr = vs_vbr;
    vbr_predictions.Q_vbr = Q_vbr;
    vbr_predictions.phi = phi;
    vbr_predictions.T = T;
    vbr_predictions.g = g;
    vbr_predictions.qmethod = qmethod;
end 
