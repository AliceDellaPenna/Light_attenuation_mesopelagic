%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An03 - Plots glider deployment map for the data used in our study
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ------------------------------------------------------------------------
% Figure formatting and toolbox setup
%% ------------------------------------------------------------------------

% Set default figure appearance for publication-quality graphics
set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')

% Add boundedline package dependencies - edit and uncomment based on your
% paths
% addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
% addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
% addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
% addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')
% 
% % Add cmocean colormap package
% addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Code\cmocean')
% 
% % Add M_Map toolbox
% addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\BioSWOT\floats\Useful_packages/m_map')
% 
% % Additional utility functions
% addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Teaching\PhD_Students\Inka')

%% ------------------------------------------------------------------------
% Define map domain
%
% Geographic bounds surrounding the Samoa glider deployment region.
%% ------------------------------------------------------------------------

lonvp = [-173.5,-171.9];
latvp = [-14.1,-12.9];

%% ------------------------------------------------------------------------
% Read glider positions
%
% Surface GPS positions recorded each time the glider surfaced.
%% ------------------------------------------------------------------------

track = readtable( ...
    'glidersurfacings.csv', ...
    'ReadVariableNames',true, ...
    'Delimiter',',');

lon_glider = track.lon_dd;
lat_glider = track.lat_dd;

%% Alternative option:
% Extract positions directly from the NetCDF file if desired.
%
% data = 'selkie_20250824T0237_Samoa_2025.nc';
% lon_glider = ncread(data,'longitude');
% lat_glider = ncread(data,'latitude');

%% ------------------------------------------------------------------------
% Read GEBCO bathymetry
%
% Load regional bathymetric data used as map background.
%% ------------------------------------------------------------------------

gebco_file = 'C:\Users\adel386\Downloads\GEBCO/gebco_2026_samoa.nc';

lon_gebco = ncread(gebco_file,'lon');
lat_gebco = ncread(gebco_file,'lat');
elev      = ncread(gebco_file,'elevation');

%% ------------------------------------------------------------------------
% Set map projection
%
% Lambert conformal projection is appropriate for regional-scale mapping
% around Samoa.
%% ------------------------------------------------------------------------

m_proj('lambert', ...
       'long',lonvp, ...
       'lat',latvp);

%% ------------------------------------------------------------------------
% Plot bathymetry and glider track
%% ------------------------------------------------------------------------

figure(1),clf

% Plot bathymetry
m_pcolor(lon_gebco,lat_gebco,elev')

% Bathymetry colour range:
% 0 m     = sea level
% -5000 m = deep ocean
caxis([-5000 0])

colorbar

% Ocean-focused colormap
cmocean('-deep')

%% ------------------------------------------------------------------------
% Overlay glider trajectory
%% ------------------------------------------------------------------------

hold on

m_plot( ...
    lon_glider, ...
    lat_glider, ...
    'wo', ...
    'markerfacecolor','k', ...
    'markersize',5)

%% ------------------------------------------------------------------------
% Add coastline and land mask
%% ------------------------------------------------------------------------

% Fill land in black
m_gshhs_f('patch',[0 0 0]);

%% ------------------------------------------------------------------------
% Add map grid and labels
%% ------------------------------------------------------------------------

m_grid( ...
    'box','fancy', ...
    'tickdir','in')

%% ------------------------------------------------------------------------
% Improve figure appearance
%% ------------------------------------------------------------------------

% Keep white background when exporting
set(gcf,'InvertHardCopy','off')
set(gcf,'Color','w')

%% ------------------------------------------------------------------------
% Save figure
%
% Creates publication-ready PNG file.
%% ------------------------------------------------------------------------

print(figure(1),'-dpng','map_glider.png')