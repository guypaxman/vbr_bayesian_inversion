%% INVERSION FUNCTIONS 4
function Vm = make_Vm(std_T,std_phi,std_g,lscale,lats,lons,zs,npts,nmod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% make_Vm = model covariance matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% make vectors of repeated lats and lons and zs
vlats = zeros(npts,1);
vlons = zeros(npts,1);
vzs = zeros(npts,1);
irow=1;
for i=1:length(lats)
    for j=1:length(lons)
        for k=1:length(zs)
            vlats(irow) = lats(i);
            vlons(irow) = lons(j);
            vzs(irow) = zs(k);
            irow=irow+1;
        end
    end
end                  
            
Vm = zeros(nmod);

for isite=1:npts
    icol = 1;
    for jsite=1:npts
        irow = (isite-1)*3+1;
        % each is a 3x3 subblock
        d = get_dist(vlats(isite),vlats(jsite),vlons(isite),vlons(jsite),...
            vzs(isite),vzs(jsite));
        fac = exp(-abs(d)/lscale);
        % T first:
        val = std_T(isite)*std_T(jsite)*fac;
        Vm(irow,icol) = val;
        irow=irow+1;
        icol=icol+1;
        % phi second;
        val = std_phi(isite)*std_phi(jsite)*fac;
        Vm(irow,icol) = val;
        irow=irow+1;
        icol=icol+1;
        % g third;
        val = std_g(isite)*std_g(jsite)*fac;
        Vm(irow,icol) = val;
        icol=icol+1; 
    end
end

end



function d = get_dist(lat1,lat2,lon1,lon2,z1,z2)

% z in km, lat/lon in deg
r1 = 6371.-z1;
r2 = 6371.-z2;
x1 = r1*sind(90.-lat1)*cosd(lon1);
y1 = r1*sind(90.-lat1)*sind(lon1);
zz1 = r1*cosd(90.-lat1);
x2 = r2*sind(90.-lat2)*cosd(lon2);
y2 = r2*sind(90.-lat2)*sind(lon2);
zz2 = r2*cosd(90.-lat2);

dx = x2-x1;
dy = y2-y1;
dz = zz2-zz1;

d = sqrt(dx*dx+dy*dy+dz*dz);

end