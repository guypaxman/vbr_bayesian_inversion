%% INVERSION FUNCTIONS 2
function F = make_F_vsonly(vs_vbr,nlat,nlon,nz,ndata,nmod,...
    T,phi,g,iTs,iphis,igs,nT,nphi,ng)

F = zeros(ndata,nmod);

% run through vs perturbations first
irow=1;
icol=1;
isite=1;
for ilat=1:nlat
    for ilon=1:nlon
        for iz=1:nz
            % dvs/dT at isite
            iphi = iphis(isite);
            ig = igs(isite);
            iTp = iTs(isite)+1;
            iTm = iTs(isite)-1;
            if (iTp>nT)
                iTp = nT;
            end
            if (iTm<1)
                iTm=1;
            end
            dvs = (vs_vbr(iz,iTp,iphi,ig)-vs_vbr(iz,iTm,iphi,ig))/...
                (T(iTp)-T(iTm));
            F(irow,icol) = dvs;
            % dvs/dphi at isite
            icol=icol+1;
            iT = iTs(isite);
            ig = igs(isite);
            iphip = iphis(isite)+1;
            iphim = iphis(isite)-1;
            if (iphip>nphi)
                iphip = nphi;
            end
            if (iphim<1)
                iphim=1;
            end
            dvs = (vs_vbr(iz,iT,iphip,ig)-vs_vbr(iz,iT,iphim,ig))/...
                (phi(iphip)-phi(iphim));
            F(irow,icol) = dvs;
            % dvs/dg at isite
            icol=icol+1;
            iT = iTs(isite);
            iphi = iphis(isite);
            igp = igs(isite)+1;
            igm = igs(isite)-1;
            if (igp>ng)
                igp = ng;
            end
            if (igm<1)
                igm=1;
            end
            dvs = (vs_vbr(iz,iT,iphi,igp)-vs_vbr(iz,iT,iphi,igm))/...
                (g(igp)-g(igm));
            F(irow,icol) = dvs;
            icol=icol+1;
            isite=isite+1;
            irow=irow+1;
        end
    end
end

% % run through Q perturbations second
% % do not reset irow 
% icol=1;
% isite=1;
% for ilat=1:nlat
%     for ilon=1:nlon
%         for iz=1:nz
%             % dQ/dT at isite
%             iphi = iphis(isite);
%             ig = igs(isite);
%             iTp = iTs(isite)+1;
%             iTm = iTs(isite)-1;
%             if (iTp>nT)
%                 iTp = nT;
%             end
%             if (iTm<1)
%                 iTm=1;
%             end
%             dQ = (Q_vbr(iz,iTp,iphi,ig)-Q_vbr(iz,iTm,iphi,ig))/...
%                 (T(iTp)-T(iTm));
%             F(irow,icol) = dQ;
%             % dvs/dphi at isite
%             icol=icol+1;
%             iT = iTs(isite);
%             ig = igs(isite);
%             iphip = iphis(isite)+1;
%             iphim = iphis(isite)-1;
%             if (iphip>nphi)
%                 iphip = nphi;
%             end
%             if (iphim<1)
%                 iphim=1;
%             end
%             dQ = (Q_vbr(iz,iT,iphip,ig)-Q_vbr(iz,iT,iphim,ig))/...
%                 (phi(iphip)-phi(iphim));
%             F(irow,icol) = dQ;
%             % dvs/dg at isite
%             icol=icol+1;
%             iT = iTs(isite);
%             iphi = iphis(isite);
%             igp = igs(isite)+1;
%             igm = igs(isite)-1;
%             if (igp>ng)
%                 igp = ng;
%             end
%             if (igm<1)
%                 igm=1;
%             end
%             dQ = (Q_vbr(iz,iT,iphi,igp)-Q_vbr(iz,iT,iphi,igm))/...
%                 (g(igp)-g(igm));
%             F(irow,icol) = dQ;
%             icol=icol+1;
%             isite=isite+1;
%             irow=irow+1;
%         end
%     end
% end


end
