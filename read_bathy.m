function [lon,lat,elev]=read_bathy(fname,lonv,latv)
%fname = '../GEBCO_data/gebco_2020_n3.457031249999993_s-86.89453125_w124.57031250000001_e234.9609375.nc';

lon_all=ncread(fname,'lon');
lat_all = ncread(fname,'lat');
idx = find(lon_all>=lonv(1) & lon_all<=lonv(2));
idy = find(lat_all>=latv(1) & lat_all<=latv(2));

lon = lon_all(idx);
lat = lat_all(idy);
elev = ncread(fname,'elevation',[idx(1) idy(1)],[idx(end)-idx(1)+1 idy(end)-idy(1)+1]);

end
