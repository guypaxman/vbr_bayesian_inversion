%% INVERSION FUNCTIONS 6
function [iTs,iphis,igs] = find_idx(X,npts,nlat,nlon,nz,phis,gs,Ts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% find_idx = find the indexes of T, phi, g in VBR output for model X
%           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iTs =  zeros(npts,1);
iphis = zeros(npts,1);
igs = zeros(npts,1);

% run through entire vs model:

imod=1;
ipts=1;
for ilat=1:nlat
    for ilon=1:nlon
        for iz=1:nz
            % T, phi, g
            [~,iT] = min(abs(Ts-X(imod)));
            iTs(ipts) = iT;
            imod=imod+1;
            [~,iphi] = min(abs(phis-X(imod)));
            iphis(ipts) = iphi;
            imod=imod+1;
            [~,ig] = min(abs(gs-X(imod)));
            igs(ipts) = ig;
            imod=imod+1;
            ipts=ipts+1;
        end
    end
end

end