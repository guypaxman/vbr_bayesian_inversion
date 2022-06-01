function results = bayesian_inversion(bayesian_settings, observations, vbr_predictions)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for the best T, phi, g for a given  tomography model using
    % non-linear Bayesian inference.
    %
    % Parameters
    % ----------
    % bayesian_settings : struct
    %   settings for the bayesian inversion. Required fields:
    %   .phi0 : float
    %       starting melt fraction
    %   .g0 : float
    %       starting log10(grain size in microns)
    %   .std_T : float
    %       temperature prior model standard deviation
    %   .std_phi : float
    %       melt fraction prior model standard deviation
    %   .std_g : float
    %       log10(grain size) prior model standard deviation
    %   .lscale : float
    %       distance scale (km) for covariance matrix
    %
    % observations : struct
    %   tomography model observations of Vs and Q, see load_data() for
    %   expected fields
    %
    % vbr_predictions : struct
    %   the forward model for Vs and Q, see load_vbr_box() for expected
    %   fields
    %
    % Returns
    % -------
    % results : struct
    %    a structure containing all the bayesian results, with the
    %    following fields:
    %
    %    .Vpo : array
    %       the full covariance array
    %    .Vpo_var : array
    %       the diagonal of the full covariance array
    %    .npar : int
    %       the number of parameters in the inversion (2)
    %    .Xk1_temp : array
    %       the max likelihood for temperature, function of lat, lon and depth
    %    .Xk1_phi : array
    %       the max likelihood for melt fraction, function of lat, lon and depth
    %    .Xk1_g : array
    %       the max likelihood for log10(grain size), function of lat, lon and depth
    %
    % Notes
    % -----
    % grain size uses a log-normal distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % extract variables from observations structure
    vs = observations.vs;
    Q = observations.Q;
    nlat = observations.nlat;
    nlon = observations.nlon;
    nz = observations.nz;
    npts = observations.npts;
    vs_err = observations.vs_err;
    Q_err = observations.Q_err;
    vs_vbr = vbr_predictions.vs_vbr;
    Q_vbr = vbr_predictions.Q_vbr;
    phi = vbr_predictions.phi;
    g = vbr_predictions.g;
    T = vbr_predictions.T;
    nT = length(T);
    ng = length(g);
    nphi = length(phi);
    
    
    % calculate log10(Q) quantities
    
    lQ = log10(Q);
    lQ_err = log10(Q_err);
    lQ_vbr= log10(Q_vbr);

    %% BEGIN MAKING INITIAL VECTORS AND MATRICES

    %  For all spatial dimensions, the loop order is same as tomography:
    %      for ilat=1:nlat ... Latitude
    %         for ilon=1:nlon ... Longitude
    %            for iz=1:nz ... Depth
    %  The order for the model parameters is T1, phi1, g1, T2, phi2, g2, ...
    %  Our total model space is (npts * number of physical parameters):
    pars = {'T'; 'phi'; 'g'};    % set the free thermodynamic parameters in the inversion
    npar = numel(pars); 
    nmod = npts*npar;
    ntype = 2; % number of observational datasets (Vs and Q)
    ndata = npts*ntype; % total data

    %  (1) Starting (X0) and Prior (Xpr) Models
    %  Let's assume X0 and Xpr are the same.
    phi0 = bayesian_settings.phi0; % melt fraction
    g0 = bayesian_settings.g0; % grain size (in microns)

    [X0,iT,iphi,ig] = make_X0(vs,nz,nlat,nlon,nmod,npts,vs_vbr,phi0,g0,phi,g,T);
    Xpr = X0;

    %  (2) Frechet kernel (F) of the starting model
    F = make_F(vs_vbr,lQ_vbr,nlat,nlon,nz,ndata,nmod,T,phi,g,iT,iphi,ig,nT,nphi,ng);

    %  (3) Make Covariance Matrices
    %      Data covariance matrix (Vd):
    Vdi = get_Vd_inv(ndata,nlat,nlon,nz,lQ_err,vs_err); % calls make_Vd

    %      Model covariance matrix (Vm):
    Vmi = get_VM(bayesian_settings, observations, npts, nmod);

    %  (4) Residual (y - f(X))
    yres = calc_yres(ndata,nlat,nlon,nz,iT,iphi,ig,vs_vbr,lQ_vbr,vs,lQ);
    
    %% BEGIN INVERSION

    %  Only Xk, F, and yres are updated in this scheme...

    irun = 1;
    Xk = X0;
    iter=1;
    iphi_idx = (2:3:nmod);  % needed for the loop, so extract here

    while irun == 1
        
        % form inversion
        % Xk1 = Xk + matA * matB
        Xk1 = update_Xk1(Xk, F, Vdi, Vmi, yres, Xpr, nmod);
        Xk1(Xk1(iphi_idx)<0)=0.; % can't have negative melt

        % calculate misfit
        [iTs,iphis,igs] = find_idx(Xk1,npts,nlat,nlon,nz,phi,g,T);
        yres = calc_yres(ndata,nlat,nlon,nz,iTs,iphis,igs,vs_vbr,lQ_vbr,vs,lQ);
        chi2 = yres'*Vdi*yres/ndata;
        display(['Chi-Squared ',num2str(chi2),' for iteration ',num2str(iter)])

        max_iter = 1;
        chi_limit = 1;
        if ((chi2 > chi_limit)&&(iter<max_iter))
            irun=1;
            
            % update F
            F = make_F(vs_vbr,lQ_vbr,nlat,nlon,nz,ndata,nmod,...
                T,phi,g,iTs,iphis,igs,nT,nphi,ng);
            % update model Xk
            Xk = Xk1;
            iter=iter+1;
        else 
            irun=0;
            if ((iter>=max_iter) && (chi2 > chi_limit))
                disp(['Did not converge after ' num2str(iter) ' iterations']);
            else
                disp(['Converged to chi-squared < ' num2str(chi_limit) ' after ' num2str(iter) ' iterations']);
            end
        end
    end

    %% EXTRACT OUTPUT

    % Calculate the posterior distribution.
    Vpo = get_Vpo(F, F', Vdi, Vmi, nmod);
    
    % Extract the diagonal (variances)
    Vpo_var = diag(Vpo);

    % Extract maximum likelihood state variables as functions of lat, lon, z
    iT_idx = (1:3:nmod);
    ig_idx = (3:3:nmod);
    Xk1_temp = Xk1(iT_idx);
    Xk1_phi = Xk1(iphi_idx); % already pulled out iphi_idx above for the loop.
    Xk1_g = Xk1(ig_idx);

    Xk1_temp = permute(reshape(Xk1_temp,[nz,nlon,nlat]),[3 2 1]);
    Xk1_phi = permute(reshape(Xk1_phi,[nz,nlon,nlat]),[3 2 1]);
    Xk1_g = permute(reshape(Xk1_g,[nz,nlon,nlat]),[3 2 1]);

    % store in the final structure to return
    results.Vpo = Vpo;
    results.Vpo_var = Vpo_var;
    results.npar = npar;
    results.Xk1_temp =  Xk1_temp;
    results.Xk1_phi =  Xk1_phi;
    results.Xk1_g =  Xk1_g;

end

% the below functions split up all of the inversion steps into isolated functions
% to help with memory usage. This ensures intermediate matrices are released from
% memory when they are no longer needed.

function Xk1 = update_Xk1(Xk, F, Vdi, Vmi, yres, Xpr, nmod)
    matA_x_matB = get_matA_x_matB(F, Vdi, Vmi, yres, Xk, Xpr, nmod);
    Xk1 = Xk + matA_x_matB;
end


function matA_x_matB = get_matA_x_matB(F, Vdi, Vmi, yres, Xk, Xpr, nmod)

    F_t = F'; % store the transpose since we need it twice
    matA =  get_Vpo(F, F_t, Vdi, Vmi, nmod);

    mat1 = F_t*Vdi*yres;
    mat2 = Vmi*(Xk-Xpr);
    matB = mat1-mat2;

    matA_x_matB = matA*matB;
end


function Vpo = get_Vpo(F, F_t, Vdi, Vmi, nmod)
    % note: F_t = F'
    matA = get_matA(F, F_t, Vdi, Vmi);
    Vpo = matA\eye(nmod);  % this it the killer inversion
end


function matA = get_matA(F, F_t, Vdi, Vmi)
    mat1 = F_t*Vdi*F;
    mat2 = Vmi;
    matA = mat1+mat2;
end


function Vmi = get_VM(bayesian_settings, observations, npts, nmod)
    % isolate this here, so that we do not keep std_T, phi, g in memory more
    % longer than needed
    std_T = ones(npts) * bayesian_settings.std_T;
    std_phi = ones(npts) * bayesian_settings.std_phi;
    std_g = ones(npts) * bayesian_settings.std_g; % log units!!
    lscale = bayesian_settings.lscale; % km
    Vm = make_Vm(std_T, std_phi, std_g, lscale, ...
                 observations.lats, observations.lons, observations.zs, ...
                 npts, nmod);
    % only the inverse is used, only return that to release the above from memory
    Vmi = Vm\eye(nmod);
end

function Vdi = get_Vd_inv(ndata, nlat,nlon,nz,lQ_err,vs_err)
    Vd = make_Vd(ndata,nlat,nlon,nz,lQ_err,vs_err);
    Vdi = Vd\eye(ndata);
end
