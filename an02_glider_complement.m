%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An02 - Plots all relevant non-light variables from the glider deployed by
% InkFish and equipped with the Wildlife computers tag
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

% The script:
% 1) Identifies bbp700 spikes
% 2) Bins observations into 5 m depth intervals
% 3) Groups profiles by time of day
% 4) Compares spike density among diel periods
% 5) Plots profiles of oxygen, bbp, and density gradients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load needed packages



addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop/Code/cmocean')

% This needs to be edited to include the boundedlinepackage and the
% subfolders github.com/kakearney/boundedline-pkg 
% Used for plotting confidence_interval


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ------------------------------------------------------------------------
% Load data and required colormaps
%% ------------------------------------------------------------------------

clear
close all

% NetCDF file
data = 'selkie_20250824T0237_Samoa_2025.nc';

%% Read variables from NetCDF file

temperature = ncread(data,'temperature');
chl         = ncread(data,'chlorophyll');
salinity    = ncread(data,'salinity');
density     = ncread(data,'density');
cdom        = ncread(data,'cdom');
oxygen      = ncread(data,'oxygen_concentration');
bs700       = ncread(data,'backscatter_700');

depth       = ncread(data,'depth');
dir         = ncread(data,'profile_direction');

% Convert Unix time (seconds since 1970) to MATLAB datenum
mtime = double(ncread(data,'time'))./60./60./24 + datenum([1970 1 1]);

profile_index = ncread(data,'profile_index');

%% Correct profile numbering
% Profile numbers restart after a mission interruption.
% Offset later profiles so all profile IDs remain unique.

profile_index(902979:end) = 252 + profile_index(902979:end);

%% ------------------------------------------------------------------------
% Visual inspection of backscatter through time and depth - uncomment if
% you want to check
%% ------------------------------------------------------------------------

% figure(3),clf
% 
% scatter(mtime,-depth,20,log10(bs700),'filled')
% 
% xlim([datenum([2025 9 8]) datenum([2025 10 6])])
% datetick('x','keeplimits')
% 
% ylim([-700 0])
% ylabel('depth (m)')
% colorbar
% 
% cmocean('turbid')

%% Define deployment period of interest
% Time interval corresponding to the light-sensor deployment
dt0 = datenum([2025 9 9]);
dtf = datenum([2025 10 4]);

%% ------------------------------------------------------------------------
% Select only observations collected during deployment period
%% ------------------------------------------------------------------------

idx = find(mtime >= dt0 & mtime <= dtf);

temperature   = temperature(idx);
chl           = chl(idx);
bs700         = bs700(idx);
depth         = depth(idx);
dir           = dir(idx);
mtime         = mtime(idx);
profile_index = profile_index(idx);
salinity      = salinity(idx);
density       = density(idx);
cdom          = cdom(idx);
oxygen        = oxygen(idx);

%% ------------------------------------------------------------------------
% Detect bbp700 spikes associated with large scatterers
%
% Idea:
% 1. Smooth each profile with a Hampel filter
% 2. Calculate profile-wide MAD (robust variability estimate)
% 3. Define spikes as observations exceeding:
%
% smoothed profile + 3 × MAD
%
% Positive residuals are interpreted as potential large scatterers.
%% ------------------------------------------------------------------------

zbins = 0:5:700; % 5 m depth bins
profy = unique(profile_index); % list of unique profiles

for prof = 1:length(profy)

idx_prof = find(profile_index == profy(prof));

% Only analyse one profile direction
% (likely downcasts)
if mean(dir(idx_prof)) > 0

% Mean timestamp of profile
dtms_prof(prof) = nanmean(mtime(idx_prof));

% Robust smoothing of bbp profile
bbp700_smooth(idx_prof) = hampel(bs700(idx_prof),15);

% Robust estimate of variability
bbp700_mad(prof) = mad(bs700(idx_prof),1);

% Spike magnitude above local background
peaks(idx_prof) = bs700(idx_prof) - bbp700_smooth(idx_prof)' - 3*bbp700_mad(prof);

else

% Ignore profiles travelling in opposite direction

dtms_prof(prof) = NaN;
bbp700_smooth(idx_prof) = NaN*bs700(idx_prof);
bbp700_mad(prof) = NaN;
peaks(idx_prof) = NaN*bs700(idx_prof);

end

end


%% ------------------------------------------------------------------------
% Bin biological and physical variables into 5 m depth intervals
%
% For each profile:
%    - calculate spike density
%    - average density
%    - average CDOM
%    - average bbp700
%    - average oxygen
%
% Result:
%    depth x profile matrices
%% ------------------------------------------------------------------------

zbins = 0:5:700;
profy = unique(profile_index);

for prof = 1:length(profy)

    idx_prof = find(profile_index == profy(prof));

    if mean(dir(idx_prof)) > 0

        dtms_prof(prof) = nanmean(mtime(idx_prof));

        for zz = 1:length(zbins)-1

            idx_db = find( ...
                depth(idx_prof) >= zbins(zz) & ...
                depth(idx_prof) <  zbins(zz+1));

            % Fraction of samples containing spikes
            peaks_binned(zz,prof) = ...
                sum(length(find(peaks(idx_prof(idx_db)))) > 0) ...
                ./ length(idx_prof(idx_db));

            % Mean physical and biogeochemical properties
            density_binned(zz,prof) = ...
                nanmean(density(idx_prof(idx_db)));

            cdom_binned(zz,prof) = ...
                nanmean(cdom(idx_prof(idx_db)));

            bbp_binned(zz,prof) = ...
                nanmean(bs700(idx_prof(idx_db)));

            oxygen_binned(zz,prof) = ...
                nanmean(oxygen(idx_prof(idx_db)));

        end

    else

        % Fill rejected profiles with NaNs

        for zz = 1:length(zbins)-1

            peaks_binned(zz,prof)   = NaN;
            density_binned(zz,prof) = NaN;
            cdom_binned(zz,prof)    = NaN;
            bbp_binned(zz,prof)     = NaN;
            oxygen_binned(zz,prof)  = NaN;

        end

    end
end

% Convert spike fraction to spike density (m^-1)
peaks_binned = peaks_binned ./ 5;


%% ------------------------------------------------------------------------
% Separate profiles into diel periods
%
% All times are UTC.
%
% morning   = 18:00-20:00
% night     = 20:00-03:00
% afternoon = 03:00-05:00
%
% (These names correspond to local sunrise/sunset conditions rather
% than clock time.)
%% ------------------------------------------------------------------------

id_good   = find(~isnan(dtms_prof));
dtms_prof = dtms_prof(id_good);

dayrefs = datevec(dtms_prof);

morning   = find(dayrefs(:,4) >= 18 & dayrefs(:,4) < 20);
noonish   = find(dayrefs(:,4) >= 20 | dayrefs(:,4) < 3);
afternoon = find(dayrefs(:,4) >= 3  & dayrefs(:,4) < 5);

%% ------------------------------------------------------------------------
% Figure 6
%
% Mean vertical distribution of bbp700 spikes
%
% Spikes are interpreted as potential large particles, zooplankton or
% micronekton targets. Here we compare the depth distribution of spike
% density among different times of day.
%
% Shading corresponds to the standard error of the mean.
%% ------------------------------------------------------------------------

figure(6),clf

boundedline(nanmean(peaks_binned(:,morning),2),-zbins(1:end-1), nanstd(peaks_binned(:,morning),[],2)./sqrt(length(morning)), 'orientation','horiz','b','alpha')

hold on

plot(nanmean(peaks_binned(:,morning),2), -zbins(1:end-1), 'b','LineWidth',2)

boundedline(nanmean(peaks_binned(:,noonish(1:end-1)),2),-zbins(1:end-1), nanstd(peaks_binned(:,noonish),[],2)./sqrt(length(noonish)), 'orientation','horiz','r','alpha')

plot(nanmean(peaks_binned(:,noonish(1:end-1)),2),-zbins(1:end-1), 'r','LineWidth',2)

boundedline(nanmean(peaks_binned(:,afternoon),2), -zbins(1:end-1),nanstd(peaks_binned(:,afternoon),[],2)./sqrt(length(afternoon)), 'orientation','horiz','k','alpha')

plot(nanmean(peaks_binned(:,afternoon),2), -zbins(1:end-1), 'k','LineWidth',2)

xlabel('density of bbp700 spikes (m^{-1})')
ylabel('depth (m)')



%% ------------------------------------------------------------------------
% Figure 8
%
% Vertical density gradient
%
% The density gradient (dρ/dz) is used as a measure of stratification.
%
%% ------------------------------------------------------------------------

figure(8),clf

boundedline(nanmean(diff(density_binned(:,morning))./5,2), -zbins(1:end-2), nanstd(diff(density_binned(:,morning))./5,[],2) ./sqrt(length(morning)), 'orientation','horiz','b','alpha')

hold on

plot(nanmean(diff(density_binned(:,morning))./5,2), -zbins(1:end-2), 'b','LineWidth',2)

boundedline(nanmean(diff(density_binned(:,noonish))./5,2), -zbins(1:end-2), nanstd(diff(density_binned(:,noonish))./5,[],2) ./sqrt(length(noonish)), 'orientation','horiz','r','alpha')

plot(nanmean(diff(density_binned(:,noonish))./5,2), -zbins(1:end-2), 'r','LineWidth',2)

boundedline(nanmean(diff(density_binned(:,afternoon))./5,2), -zbins(1:end-2), nanstd(diff(density_binned(:,afternoon))./5,[],2) ./sqrt(length(afternoon)), 'orientation','horiz','k','alpha')

plot(nanmean(diff(density_binned(:,afternoon))./5,2), -zbins(1:end-2), 'k','LineWidth',2)

xlabel('vertical density gradient (kg m^{-4})')
ylabel('depth (m)')

%% ------------------------------------------------------------------------
% Figure 11
%
% Mean bbp700 profiles
%% ------------------------------------------------------------------------

figure(11),clf


boundedline(nanmean(bbp_binned(:,morning),2), -zbins(1:end-1), nanstd(bbp_binned(:,morning),[],2)./sqrt(length(morning)), 'orientation','horiz','b','alpha')
hold on
plot(nanmean(bbp_binned(:,morning),2), -zbins(1:end-1),'b','LineWidth',2)


boundedline(nanmean(bbp_binned(:,noonish),2), -zbins(1:end-1), nanstd(bbp_binned(:,noonish),[],2)./sqrt(length(noonish)), 'orientation','horiz','r','alpha')

plot(nanmean(bbp_binned(:,noonish),2), -zbins(1:end-1),'r','LineWidth',2)


boundedline(nanmean(bbp_binned(:,afternoon),2), -zbins(1:end-1), nanstd(bbp_binned(:,afternoon),[],2)./sqrt(length(afternoon)), 'orientation','horiz','k','alpha')

plot(nanmean(bbp_binned(:,afternoon),2), -zbins(1:end-1),'k','LineWidth',2)
xlim([0.5 2]*10^(-4))
xlabel('bbp700 (m^{-1})')
ylabel('depth (m)')

%% ------------------------------------------------------------------------
% Figure 12
%
% Mean oxygen concentration profiles
%% ------------------------------------------------------------------------

figure(12),clf

boundedline(nanmean(oxygen_binned(:,morning),2), -zbins(1:end-1), nanstd(oxygen_binned(:,morning),[],2)./sqrt(length(morning)), 'orientation','horiz','b','alpha')

hold on

plot(nanmean(oxygen_binned(:,morning),2),-zbins(1:end-1), 'b','LineWidth',2)

boundedline(nanmean(oxygen_binned(:,noonish),2), -zbins(1:end-1), nanstd(oxygen_binned(:,noonish),[],2)./sqrt(length(noonish)), 'orientation','horiz','r','alpha')

plot(nanmean(oxygen_binned(:,noonish),2), -zbins(1:end-1), 'r','LineWidth',2)

boundedline(nanmean(oxygen_binned(:,afternoon),2),-zbins(1:end-1), nanstd(oxygen_binned(:,afternoon),[],2)./sqrt(length(afternoon)), 'orientation','horiz','k','alpha')

plot(nanmean(oxygen_binned(:,afternoon),2), -zbins(1:end-1), 'k','LineWidth',2)

xlabel('oxygen concentration (\mu mol kg^{-1})')
ylabel('depth (m)')
