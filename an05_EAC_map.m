%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An05 - Maps locations of measurements from the Triaxus deployment in the
% Tasman Sea / East Australian Current
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set graphical preferences
set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')

% Select the paths to all needed scripts - this will likely have to change

addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Code\cmocean')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\BioSWOT\floats\Useful_packages/m_map')

% Select the geogrpahical area of interest with a buffer to avoid blank
% spaces around the map
lonv = [153,156];
latv = [-39,-36.5];
lonvp = [153,155.9];
latvp = [-39,-36.6];

% Reads acoustic data to extract locations
data120 = readtable('IN2023_V06_D20231020_T022530_Transducer_120000_sliced_transectb.csv');

% Reads satellite data - the folders will need to change depending on your
% path
folder_data = 'C:\Users\adel386\Downloads\';
SLA_map = sprintf('%s/dt_global_allsat_phy_l4_20231021_20241023.nc',folder_data);
SST_map = sprintf('%s/20231021120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc',folder_data);

% Reads satellite maps
[lon_sst,lat_sst,sst]= read_SST_dt(SST_map,lonv,latv);
[lon_sla,lat_sla,sla]= read_SLA_dt(SLA_map,lonv,latv);
[long,latg] = meshgrid(lon_sla,lat_sla);
m_proj('lambert','long',lonvp,'lat',latvp);

% Loads all Triaxus data - probably not all of this is needed but it's a
% good way to check that everything is there
flist = glob('in2023*.mat');
all_mtimes = [];
all_abundance = [];
all_biomass = [];
all_pressure = [];
all_SMEP = [];
all_NBSS = [];
all_lons = [];
all_lats = [];
for ff=2:length(flist)
    f = load(flist{ff});
    all_abundance = [all_abundance;f.s.Abundance];
    all_biomass = [all_biomass; f.s.Biomass];
    all_mtimes = [all_mtimes; f.s.datenum];
    all_pressure = [all_pressure; f.s.pressure];
    all_SMEP = [all_SMEP; f.s.SMEP];
    all_NBSS = [all_NBSS; f.s.NBSS_Slope];
    all_lons = [all_lons; f.s.longitude];
    all_lats = [all_lats; f.s.latitude];
end

%% Select the data we are actually using
dt0 = datenum([2023 10 20]);
dtf = datenum([2023 10 22]);
dayrefs_t = datevec(all_mtimes);
morning_t = find(dayrefs_t(:,4)>=20 & dayrefs_t(:,4)<23);
noonish_t = find(dayrefs_t(:,4)>=23 | dayrefs_t(:,4)<4);
afternoon_t = find(dayrefs_t(:,4)>=4 & dayrefs_t(:,4)<8);
idx_t = find(all_mtimes>=dt0 & all_mtimes<=dtf & dayrefs_t(:,4)>=20 | dayrefs_t(:,4)<8);
dayrefs_e = datevec(data120.Time_E);
morning_e = find(dayrefs_e(:,4)>=20 & dayrefs_e(:,4)<23);
noonish_e = find(dayrefs_e(:,4)>=23 | dayrefs_e(:,4)<4);
afternoon_e = find(dayrefs_e(:,4)>=4 & dayrefs_e(:,4)<8);
idx_e = find(datenum(data120.Time_E)>=dt0 & datenum(data120.Time_E)<=dtf & dayrefs_e(:,4)>=20 | dayrefs_e(:,4)<8)

%% Plot map
figure(1),clf,m_pcolor(lon_sst,lat_sst,sst-273.15)
clim([16 18])
colorbar
cmocean('dense','negative')
hold on,m_contour(lon_sla,lat_sla,sla,'w','linewidth',2)
hold on,m_contour(lon_sla,lat_sla,sla,'k')
hold on,m_plot(data120.Lon_S(idx_e),data120.Lat_S(idx_e),'w','linewidth',10)
hold on,m_plot(all_lons(idx_t),all_lats(idx_t),'r.','markersize',12)
m_grid('box', 'fancy', 'tickdir', 'in');
set(gcf, 'InvertHardCopy', 'off')
set(gcf, 'Color', 'w')
print(figure(1),'-dpng','EAC_map_eddy.png')


%% Tiny world inlet
figure(2),clf
%m_proj('azimuthal equal-area','radius',156,'lat',-46,'long',65,'rot',0);
ax4 = subplot(1,1,1);
m_proj('ortho','lat',-38','long',154');
colormap(ax4,[m_colmap('blues')]);  
      m_elev('image');
      m_coast('patch',[0.75 0.65 0.65]);
%print(figure(2),'-dpng','tiny_world_EAC.png')
%m_grid('xticklabel',[],'yticklabel',[],'linestyle','-','ytick',[-60:30:60]);