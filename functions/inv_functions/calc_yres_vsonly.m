%% INVERSION FUNCTIONS 5
function yres = calc_yres_vsonly(ndata,nlat,nlon,nz,iT,iphi,ig,vs_vbr,vs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calc_res = calculate residual matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yres = zeros(ndata,1);

% vs first
ipts = 1;
irow = 1;
for ilat=1:nlat
    for ilon=1:nlon
        for iz=1:nz
            val = vs_vbr(iz,iT(ipts),iphi(ipts),ig(ipts));
            yres(irow) = vs(ilat,ilon,iz) - val;
            ipts=ipts+1;
            irow=irow+1;
        end
    end
end

% % Q second - do not reset irow
% ipts = 1;
% for ilat=1:nlat
%     for ilon=1:nlon
%         for iz=1:nz
%             val = lQ_vbr(iz,iT(ipts),iphi(ipts),ig(ipts));
%             yres(irow) = lQ(ilat,ilon,iz) - val;
%             ipts=ipts+1;
%             irow=irow+1;
%         end
%     end
% end

end
