function  [lon,lat,sst]= read_SST_dt(filename,lonv,latv)
mtime = datenum([1981 1 1])+ncread(filename,'time')./60./60./24;
lons = ncread(filename,'lon');
lats = ncread(filename,'lat');

if lonv(1)<=180 
indxx = find(lons>lonv(1) & lons<lonv(2));
indxy = find(lats>latv(1) & lats<latv(2));

lon_map = ncread(filename,'lon',indxx(1),indxx(end)-indxx(1));
lat_map = ncread(filename,'lat',indxy(1),indxy(end)-indxy(1));

sst_map = ncread(filename,'analysed_sst',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]);

    

if lonv(1)<180 & lonv(2)>180
    lonv2 = lonv(2)-360;
    indxx2 = find(lons<lonv2 & lons>-180);
    lon_map2 = ncread(filename,'lon',indxx2(1),indxx2(end)-indxx2(1));
    sst_map2= ncread(filename,'analysed_sst',[indxx2(1) indxy(1) 1],[indxx2(end)-indxx2(1),indxy(end)-indxy(1) 1]);

    lon_map = [lon_map; lon_map2+360];
    sst_map_ext = [sst_map; sst_map2];
else
    sst_map_ext = sst_map;
end
else
    indxx = find(lons>lonv(1)-360 & lons<lonv(2)-360);
    indxy = find(lats>latv(1) & lats<latv(2));

lon_map = ncread(filename,'lon',indxx(1),indxx(end)-indxx(1))+360;
lat_map = ncread(filename,'lat',indxy(1),indxy(end)-indxy(1));

sst_map = ncread(filename,'analysed_sst',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]); 
sst_map_ext = sst_map;
end

sst_map = sst_map_ext';
sst = sst_map;
lon = lon_map;
lat = lat_map;
end