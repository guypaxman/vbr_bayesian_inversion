%% RUN A BAYESIAN INVERSION  
%
%  Here we will solve for the best T, phi, g for a given 
%  tomography model using non-linear Bayesian inference

% TO DO: delete this file?

clc
clear
addpath(genpath('./functions'))
addpath('./GLAD25')

%% INPUT FILES

% regional limits
elim = -10;
wlim = -80.;
nlim = 86.;
slim = 56.;
xy_res = 2;    % spatial resolution [decimal degrees; multiple of 0.5]

% The first set is the tomography model (vs and Q)
GLAD25_Vs = load('Vs_Model.mat');
QL6_Q = load('QL6_Model.mat');
lats = GLAD25_Vs.Vs_Model.Latitude;
lons = GLAD25_Vs.Vs_Model.Longitude;

% Set the vs and Q models as "vs" and "Q"
% trim model to correct area and resolution
idxlat = find((lats>=slim)&(lats<=nlim));
idxlon = find((lons>=wlim)&(lons<=elim));
idxlat = idxlat(1:xy_res*2:end);
idxlon = idxlon(1:xy_res*2:end);
lats = lats(idxlat);
lons = lons(idxlon);
vs = GLAD25_Vs.Vs_Model.Vs(idxlat,idxlon,4:40); % km/s
Q = QL6_Q.QL6_Model.Q(idxlat,idxlon,4:40);
zs = GLAD25_Vs.Vs_Model.Depth(4:40); % remove top 3 slices (10, 20, 30 km)
nlat = length(lats);
nlon = length(lons);
nz = length(zs);
npts = nlat*nlon*nz;

% Set errors in Vs and Q
lQ = log10(Q); % try log10 of Q
vs_err = ones(size(vs))*0.05;
Q_err = ones(size(Q))*10.;
lQ_err = log10(Q_err);

%  Now we import the parameter sweep 
%  This should be made at a wide range of T, phi, g at fine increments
%  and should predict vs and Q at the right frequency band
%  Use make_sweep.m to do this and import that file here
%  Make sure the depth of sweep matches the depth in the tomography model
sweep_in = 'sweep_box.mat';
load(sweep_in);
params = make_param_grid(sweep.state_names, sweep);
qmethod = 'xfit_premelt';
%  Extract the VBR box
T = unique(params.T);
phi = unique(params.phi);
g = log10(unique(params.gs)); %% log grain size
% g = unique(params.gs); %% linear grain size
nT = length(T);
nphi = length(phi);
ng = length(g);
vs_vbr = zeros(nz,nT,nphi,ng);
Q_vbr = zeros(nz,nT,nphi,ng);
lQ_vbr = zeros(nz,nT,nphi,ng);
for iT=1:nT
    for iphi=1:nphi
        for ig=1:ng
            for iz=1:nz
                vs_vbr(iz,iT,iphi,ig) = sweep.Box(iT,iphi,ig).(qmethod).Vs(iz)*1.e-3; % km/s
                Q_vbr(iz,iT,iphi,ig) = sweep.Box(iT,iphi,ig).(qmethod).Q(iz);
                lQ_vbr(iz,iT,iphi,ig) = log10(Q_vbr(iz,iT,iphi,ig));
            end
        end
    end
end

clear GLAD25_Vs QL6_Q sweep; 

%% BEGIN MAKING INITIAL VECTORS AND MATRICES

%  For all spatial dimensions, the loop order is same as tomography:
%      for ilat=1:nlat ... Latitude
%         for ilon=1:nlon ... Longitude
%            for iz=1:nz ... Depth
%  The order for the model parameters is T1, phi1, g1, T2, phi2, g2, ...
%  Our total model space is (npts * number of physical parameters):
npar = 3;
pars = ['T','phi','g'];
nmod = npts*npar;
ntype = 2; % type of data (vs and Q)
ndata = npts*ntype; % total data

%  (1) Starting (X0) and Prior (Xpr) Models
%  Let's assume X0 and Xpr are the same.
phi0 = 0.01; % melt fraction
g0 = log10(1000.); % 1 mm (1000 microns)
% g0 = 4000.; % 5 mm

[X0,iT,iphi,ig] = make_X0_xfit_premelt(vs,nz,nlat,nlon,nmod,npts,vs_vbr,phi0,g0,phi,g,T);
Xpr = X0;

%  (2) Frechet kernel (F) of the starting model
F = make_F(vs_vbr,lQ_vbr,nlat,nlon,nz,ndata,nmod,...
    T,phi,g,iT,iphi,ig,nT,nphi,ng);

%  (3) Make Covariance Matrices
%      Data covariance matrix (Vd):
Vd = make_Vd(ndata,nlat,nlon,nz,lQ_err,vs_err);
%      Model covariance matrix (Vm):
std_T = ones(npts)*50.;
std_phi = ones(npts)*0.0025;
std_g = ones(npts)*0.2; % log units
% std_g = ones(npts)*2000; % linear units (microns)
lscale = 200.; % km
Vm = make_Vm(std_T,std_phi,std_g,lscale,lats,lons,zs,npts,nmod);

%  (4) Residual (y - f(X))
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

while irun == 1
    
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
    
%     ig0s = find(Xk1(ig_idx)<0.1e3); % can't have small grain size
%     Xk1(ig_idx(ig0s))=0.1e3;
%    
    % calculate misfit    
    [iTs,iphis,igs] = find_idx(Xk1,npts,nlat,nlon,nz,phi,g,T);
    yres = calc_yres(ndata,nlat,nlon,nz,iTs,iphis,igs,vs_vbr,lQ_vbr,vs,lQ);
    chi2 = yres'*Vdi*yres/ndata;
    display(['Chi-Squared',string(chi2),'for iteration',string(iter)])

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

save('Temp.mat','Temp');
movefile('./Temp.mat', './state_variables_xfit_premelt/');    
save('Phi.mat','Phi');
movefile('./Phi.mat', './state_variables_xfit_premelt/');    
save('GrainSize.mat','GrainSize');
movefile('./GrainSize.mat', './state_variables_xfit_premelt/');   

save('Vpo.mat','Vpo_var');
movefile('./Vpo.mat', './state_variables_xfit_premelt/');    
save('Vpo_site1.mat','Vpo_site1');
movefile('./Vpo_site1.mat', './state_variables_xfit_premelt/');    
save('Vpo_site2.mat','Vpo_site2');
movefile('./Vpo_site2.mat', './state_variables_xfit_premelt/');    
save('Vpo_site3.mat','Vpo_site3');
movefile('./Vpo_site3.mat', './state_variables_xfit_premelt/');        
