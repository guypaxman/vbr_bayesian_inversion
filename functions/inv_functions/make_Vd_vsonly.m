%% INVERSION FUNCTIONS 3
function Vd = make_Vd_vsonly(ndata,nlat,nlon,nz,vs_err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% make_Vd = data covariance matrix, assume uncorrelated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vd = zeros(ndata);

% vs first
ipts=1;
for i=1:nlat
    for j=1:nlon
        for k=1:nz
            Vd(ipts,ipts) = vs_err(i,j,k)*vs_err(i,j,k);
            ipts=ipts+1;
        end
    end
end

% % Q next
% 
% for i=1:nlat
%     for j=1:nlon
%         for k=1:nz
%             Vd(ipts,ipts) = Q_err(i,j,k)*Q_err(i,j,k);
%             ipts=ipts+1;
%         end
%     end
% end

end

