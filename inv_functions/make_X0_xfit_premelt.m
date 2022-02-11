%% INVERSION FUNCTIONS 1
function [X0,iTs,iphis,igs] = ...
    make_X0_xfit_premelt(vs,nz,nlat,nlon,nmod,npts,vs_vbr,phi0,g0,phis,gs,Ts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% make_X0 = assume all changes due to temperature and assume prior gs and 
%           phi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0 = zeros(nmod,1);
iTs =  zeros(npts,1);
iphis = zeros(npts,1);
igs = zeros(npts,1);

[~,iphi] = min(abs(phis-phi0));
[~,ig] = min(abs(gs-g0));

% run through entire vs model:

imod=1;
ipts=1;
for ilat=1:nlat
    for ilon=1:nlon
        for iz=1:nz
            val = vs(ilat,ilon,iz);
            [~,iT] = min(abs(vs_vbr(iz,:,iphi,ig)-val));
            X0(imod) = Ts(iT)-273;
            imod=imod+1;
            % now do phi:
            X0(imod) = phi0;
            imod=imod+1;
            % now do g:
            X0(imod) = g0;
            imod=imod+1;
            
            % record indices:
            iTs(ipts) = iT;
            iphis(ipts) = iphi;
            igs(ipts) = ig;
            ipts=ipts+1;            
        end
    end
end

end