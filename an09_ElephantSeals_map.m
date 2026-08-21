%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An09 - Plot map showing the matches of ships and seals from the Kerguelen
% region.
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Graphical preferences and path

% Graphical preferences
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');

set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')

% Set paths - these will likely need to change
addpath('/Users/adel386/OneDrive - The University of Auckland/Desktop/Code')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\BioSWOT\floats\Useful_packages/m_map')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop/Code/cmocean')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\Light_Mesopelagic\GEBCO_Kerguelen\GEBCO_23_Apr_2025_e8d27c89343e')


%% Read data 

% Read data from the spreadsheet including the matches between ships and
% seals
matches = readtable('seals_shipacoustics_match_results_2025-03-13-FILTER.csv');

% Read the IMOS acoustic data from the ship with the highest number of
% matches - the full path will likely need to change depending on where
% this file is stored
magic_ship = 'IMOS_SOOP-BA_AE_20121215T021630Z_VHLU_FV02_Austral-Leader-II-ES60-38_END-20121222T190951Z_C-20160222T052709Z.nc';
magic_ship_full = 'D:\Acoustic\IMOS\38kHz\IMOS_SOOP-BA_AE_20121215T021630Z_VHLU_FV02_Austral-Leader-II-ES60-38_END-20121222T190951Z_C-20160222T052709Z.nc';
magic_ship_full(1) = 'E';
lon_ship = ncread(magic_ship,'LONGITUDE');
lat_ship = ncread(magic_ship,'LATITUDE');
mtime = ncread(magic_ship,'TIME')+datenum([1950 1 1]);

% Read SST data
sst_map = '20121222120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc';
sst = ncread(sst_map,'analysed_sst')-273.15;
lon_sst = ncread(sst_map,'lon');
lat_sst = ncread(sst_map,'lat');

%% Plot figure
m_proj('lambert','long',[67 84],'lat',[-54 -46]);
figure(1),clf,m_pcolor(lon_sst,lat_sst,sst')
cmocean('dense','negative')
caxis([1 8])
hold on,m_contour(lon_sst,lat_sst,sst',[4,8,10],'k')
cb = colorbar;
title(cb,'SST (^\circ C)')

%% Add bathymetry
bathy_file = 'gebco_2024_n-40.0_s-55.0_w60.0_e90.0.nc';
[lon_b,lat_b,elev]=read_bathy(bathy_file,[67 84],[-54 -46]);
hold on,m_contour(lon_b,lat_b,elev',[-1000 -200],'color',[0.9 0.9 0.9],'linewidth',1.5);

%% Ship information
idx_ship = find(mtime>=datenum([2012 12 20]) & mtime<= datenum([2012 12 22 12 0 0]));
idx_ship_match = find(strcmp(matches.ship_fname,magic_ship_full));
%hold on,m_plot(lon_ship,lat_ship,'k') 
hold on,m_plot(lon_ship(idx_ship),lat_ship(idx_ship),'ko','markerfacecolor','k')  % overall track
hold on,m_plot(matches.ship_lon(idx_ship_match),matches.ship_lat(idx_ship_match),'ro','markerfacecolor','r'); % Daily locations

%% Seal information - reading and plotting from the MEOP dataset using the spreadsheet as a map
useals = unique(matches.id(idx_ship_match));
for qq=1:length(useals)
seal_track = sprintf('seals/ncARGO_traj/%s_traj.nc',useals{qq});
lon_seal = ncread(seal_track,'LONGITUDE');
lat_seal = ncread(seal_track,'LATITUDE');
mtime_seal = ncread(seal_track,'TIME')+datenum([1950 1 1]);
idx_match = matches.date(find(strcmp(matches.id,useals{qq}) & strcmp(matches.ship_fname,magic_ship_full)));
idx_idx = find(strcmp(matches.id,useals{qq}) & strcmp(matches.ship_fname,magic_ship_full));
idx_time_seal = find(mtime_seal>=datenum(min(idx_match)-1.5) & mtime_seal<=datenum(max(idx_match)+1.5));
hold on,m_plot(lon_seal(idx_time_seal),lat_seal(idx_time_seal),'w.','markersize',20) % Plot track
hold on,m_plot(matches.mean_lon(idx_idx),matches.mean_lat(idx_idx),'bo','markerfacecolor','b') % Daily location
end
m_gshhs_f('patch',[1 1 1]) 
m_grid('box','fancy');
print(figure(1),'-dpng','match_map1.png')

%% Tiny world inlet
figure(2),clf
ax4 = subplot(1,1,1);
m_proj('ortho','lat',-46','long',65');
colormap(ax4,[m_colmap('blues')]);  
      m_elev('image');
      m_coast('patch',[0.75 0.65 0.65]);
print(figure(2),'-dpng','tiny_world.png')


