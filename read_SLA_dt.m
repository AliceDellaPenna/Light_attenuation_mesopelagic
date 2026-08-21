function [lon,lat,sla,u,v] = read_SLA_dt(filename,lonv,latv)



lons = ncread(filename,'longitude');
lats = ncread(filename,'latitude');

if lonv(1)<=180 
indxx = find(lons>lonv(1) & lons<lonv(2));
indxy = find(lats>latv(1) & lats<latv(2));

lon_map = ncread(filename,'longitude',indxx(1),indxx(end)-indxx(1));
lat_map = ncread(filename,'latitude',indxy(1),indxy(end)-indxy(1));

sla_map = ncread(filename,'sla',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]);
ugos_map = ncread(filename,'ugos',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]);
vgos_map = ncread(filename,'vgos',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]);

if lonv(1)<180 & lonv(2)>180
    lonv2 = lonv(2)-360;
    indxx2 = find(lons<lonv2 & lons>-180);
    lon_map2 = ncread(filename,'longitude',indxx2(1),indxx2(end)-indxx2(1));
    sla_map2= ncread(filename,'sla',[indxx2(1) indxy(1) 1],[indxx2(end)-indxx2(1),indxy(end)-indxy(1) 1]);
    ugos_map2= ncread(filename,'ugos',[indxx2(1) indxy(1) 1],[indxx2(end)-indxx2(1),indxy(end)-indxy(1) 1]);
    vgos_map2= ncread(filename,'vgos',[indxx2(1) indxy(1) 1],[indxx2(end)-indxx2(1),indxy(end)-indxy(1) 1]);

    lon_map = [lon_map; lon_map2+360];
    sla_map_ext = [sla_map; sla_map2];
    ugos_map_ext = [ugos_map; ugos_map2];
    vgos_map_ext = [vgos_map; vgos_map2];
    
else
    sla_map_ext = sla_map;
    ugos_map_ext = ugos_map;
   vgos_map_ext = vgos_map;
end
else
    indxx = find(lons>lonv(1)-360 & lons<lonv(2)-360);
    indxy = find(lats>latv(1) & lats<latv(2));

lon_map = ncread(filename,'longitude',indxx(1),indxx(end)-indxx(1))+360;
lat_map = ncread(filename,'latitude',indxy(1),indxy(end)-indxy(1));

sla_map = ncread(filename,'sla',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]); 
ugos_map = ncread(filename,'ugos',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]); 
vgos_map = ncread(filename,'vgos',[indxx(1) indxy(1) 1],[indxx(end)-indxx(1),indxy(end)-indxy(1) 1]); 

sla_map_ext = sla_map;
    ugos_map_ext = ugos_map;
   vgos_map_ext = vgos_map;
end

sla_map = sla_map_ext';
u = ugos_map_ext';
v = vgos_map_ext';
sla = sla_map;
lon = lon_map;
lat = lat_map;

end