%% RUN A BAYESIAN INVERSION  
%
%  Here we will solve for the best T, phi, g for a given 
%  tomography model using non-linear Bayesian inference


function bayesian_inversion(observations, vbr_predictions)
    
    addpath('./functions')
    addpath('./inv_functions')

    % extract variables from observations structure
    vs = observations.vs;
    Q = observations.Q;
    lons = observations.lons;
    lats = observations.lats;
    zs = observations.zs;
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
    pars = ['T','phi','g'];    % set the free thermodynamic parameters in the inversion
    npar = numel(pars); 
    nmod = npts*npar;
    ntype = 2; % number of observational datasets (Vs and Q)
    ndata = npts*ntype; % total data

    % TO DO: move up to user-defined
    %  (1) Starting (X0) and Prior (Xpr) Models
    %  Let's assume X0 and Xpr are the same.
    phi0 = 0.01; % melt fraction
    g0 = log10(1000.); % 1 mm (1000 microns)

    % TO DO: switch following function to just make_X0
    disp('building X0')
    [X0,iT,iphi,ig] = make_X0_xfit_premelt(vs,nz,nlat,nlon,nmod,npts,vs_vbr,phi0,g0,phi,g,T);
    Xpr = X0;

    %  (2) Frechet kernel (F) of the starting model
    disp('building Frechet kernel')
    F = make_F(vs_vbr,lQ_vbr,nlat,nlon,nz,ndata,nmod,T,phi,g,iT,iphi,ig,nT,nphi,ng);

    %  (3) Make Covariance Matrices
    %      Data covariance matrix (Vd):
    disp('building data covariance matrix')
    Vd = make_Vd(ndata,nlat,nlon,nz,lQ_err,vs_err);
    
    % TO DO: move up std, lscale values to user-defined
    %      Model covariance matrix (Vm):
    std_T = ones(npts)*50.;
    std_phi = ones(npts)*0.0025;
    std_g = ones(npts)*0.2; % log units
    % std_g = ones(npts)*2000; % linear units (microns)
    lscale = 200.; % km
    disp('building model covariance matrix')
    Vm = make_Vm(std_T,std_phi,std_g,lscale,lats,lons,zs,npts,nmod);

    %  (4) Residual (y - f(X))
    disp('calculating initial residual')
    yres = calc_yres(ndata,nlat,nlon,nz,iT,iphi,ig,vs_vbr,lQ_vbr,vs,lQ);

    %% BEGIN INVERSION

    %  Only Xk, F, and yres are updated in this scheme...

    irun = 1;
    % invert some matrices:
    Vdi = Vd\eye(ndata);
    Vmi = Vm\eye(nmod);
    Xk = X0;

    iter=1;

    iphi_idx = (2:3:nmod);
    iT_idx = (1:3:nmod);
    ig_idx = (3:3:nmod);
    disp('Beginning inversion')
    while irun == 1
        disp(['Itration ', num2str(iter)])
        % form inversion

        mat1 = F'*Vdi*F;
        mat2 = Vmi;
        matA = mat1+mat2;
        matA = matA\eye(nmod);
        
        mat1 = F'*Vdi*yres;
        mat2 = Vmi*(Xk-Xpr);
        matB = mat1-mat2;
        
        Xk1 = Xk + matA*matB;

        iphi0s = find(Xk1(iphi_idx)<0); % can't have negative melt
        Xk1(iphi_idx(iphi0s))=0.;

        % calculate misfit    
        [iTs,iphis,igs] = find_idx(Xk1,npts,nlat,nlon,nz,phi,g,T);
        yres = calc_yres(ndata,nlat,nlon,nz,iTs,iphis,igs,vs_vbr,lQ_vbr,vs,lQ);
        chi2 = yres'*Vdi*yres/ndata;
        display(['    Chi-Squared ',num2str(chi2),' for iteration ',num2str(iter)])

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
            if (iter>=max_iter)
                disp(['Did not converge after ' num2str(iter) ' iterations']);
            else
                disp(['Converged to chi-squared < ' num2str(chi_limit) ' after ' num2str(iter) ' iterations']);
            end
        end
    end

    %% EXTRACT OUTPUT
    % TO DO: decide what in the following should be in the repository
    % Calculate the posterior distribution.
    mat1 = F'*Vdi*F;
    mat2 = Vmi;
    matA = mat1+mat2;
    Vpo = matA\eye(nmod);
    Vpo_var = diag(Vpo);

    % Extract covariances for chosen sites
    lat_site1 = 78;
    lon_site1 = -42;
    % TO DO: make a function out of this
    start_idx1 = (npar*nz*nlon*(find(lats==lat_site1)-1)) + (npar*nz*(find(lons==lon_site1)-1)) + 1;
    end_idx1 = (npar*nz*nlon*(find(lats==lat_site1)-1)) + (npar*nz*find(lons==lon_site1));
    Vpo_site1 = Vpo(start_idx1:end_idx1,:);

    lat_site2 = 68;
    lon_site2 = -34;
    start_idx2 = (npar*nz*nlon*(find(lats==lat_site2)-1)) + (npar*nz*(find(lons==lon_site2)-1)) + 1;
    end_idx2 = (npar*nz*nlon*(find(lats==lat_site2)-1)) + (npar*nz*find(lons==lon_site2));
    Vpo_site2 = Vpo(start_idx2:end_idx2,:);

    lat_site3 = 66;
    lon_site3 = -48;
    start_idx3 = (npar*nz*nlon*(find(lats==lat_site3)-1)) + (npar*nz*(find(lons==lon_site3)-1)) + 1;
    end_idx3 = (npar*nz*nlon*(find(lats==lat_site3)-1)) + (npar*nz*find(lons==lon_site3));
    Vpo_site3 = Vpo(start_idx3:end_idx3,:);

    % Extract maximum likelihood state variables as functions of lat, lon, z
    Xk1_temp = Xk1(iT_idx);
    Xk1_phi = Xk1(iphi_idx);
    Xk1_g = Xk1(ig_idx);

    Xk1_temp = permute(reshape(Xk1_temp,[nz,nlon,nlat]),[3 2 1]);
    Xk1_phi = permute(reshape(Xk1_phi,[nz,nlon,nlat]),[3 2 1]);
    Xk1_g = permute(reshape(Xk1_g,[nz,nlon,nlat]),[3 2 1]);
            
    Temp.latitude = lats;
    Temp.longitude = lons;
    Temp.depth = zs;
    Temp.temperature = Xk1_temp;

    Phi.latitude = lats;
    Phi.longitude = lons;
    Phi.depth = zs;
    Phi.meltfraction = Xk1_phi;

    GrainSize.latitude = lats;
    GrainSize.longitude = lons;
    GrainSize.depth = zs;
    GrainSize.grainsize = Xk1_g;

    save_dir = ['results/', 'state_variables_', vbr_predictions.qmethod, '/'];
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    
    save([save_dir, 'Temp.mat'],'Temp');
    save([save_dir, 'Phi.mat'],'Phi');
    save([save_dir, 'GrainSize.mat'],'GrainSize');
    save([save_dir, 'Vpo.mat'],'Vpo_var');
    save([save_dir, 'Vpo_site1.mat'],'Vpo_site1');
    save([save_dir, 'Vpo_site2.mat'],'Vpo_site2');
    save([save_dir, 'Vpo_site3.mat'],'Vpo_site3');
  
end
