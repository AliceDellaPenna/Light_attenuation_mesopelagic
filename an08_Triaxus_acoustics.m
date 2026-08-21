%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An08 - Plot sonar acoustics matching Triaxus deployment in the
% Tasman Sea / East Australian Current
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Graphical preferences and path

% Graphical preferences
set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')

% Paths - these will need to change
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')

%% Load data

% Multi frequency sonar acoustics
data120 = readtable('IN2023_V06_D20231020_T022530_Transducer_120000_sliced_transect.csv');
data70 = readtable('IN2023_V06_D20231020_T022530_Transducer_70000_sliced_transect.csv');
data38 = readtable('IN2023_V06_D20231020_T022530_Transducer_38000_sliced_transect.csv');
data18 = readtable('IN2023_V06_D20231020_T022530_Transducer_18000_sliced_transect.csv');

% Triaxus data: load and concatenate from four deployments

flist = glob('in2023*.mat');


all_mtimes = [];
all_abundance = [];
all_biomass = [];
all_pressure = [];
all_SMEP = [];
all_NBSS = [];
all_lons = [];
all_lats = [];

for ff=1:length(flist)

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


%% Selects the time of interest and separate by the time of the day
% Remember that the times are in UTC and that this needs to be done for
% both datasets

dt0 = datenum([2023 10 20 14 50 49]);
dtf = datenum([2023 10 21 11 32 53]);

dayrefs_t = datevec(all_mtimes);

morning_t = find(dayrefs_t(:,4)>=20 & dayrefs_t(:,4)<23);
noonish_t = find(dayrefs_t(:,4)>=23 | dayrefs_t(:,4)<4);
afternoon_t = find(dayrefs_t(:,4)>=4 & dayrefs_t(:,4)<8);

idx_t = find(all_mtimes>=dt0 & all_mtimes<=dtf & dayrefs_t(:,4)>=20 | dayrefs_t(:,4)<8);

dayrefs_e = datevec(data120.Time_E);

morning_e = find(dayrefs_e(:,4)>=20 & dayrefs_e(:,4)<23);
noonish_e = find(dayrefs_e(:,4)>=23 | dayrefs_e(:,4)<4);
afternoon_e = find(dayrefs_e(:,4)>=4 & dayrefs_e(:,4)<8);

idx_e = find(datenum(data120.Time_E)>=dt0 & datenum(data120.Time_E)<=dtf & (dayrefs_e(:,4)>=20 | dayrefs_e(:,4)<8));

dayrefs_e18 = datevec(data18.Time_E);
idx_e18 = find(datenum(data18.Time_E)>=dt0 & datenum(data18.Time_E)<=dtf & (dayrefs_e18(:,4)>=20 | dayrefs_e18(:,4)<8));

%%  Filter the acoustic datasets based on the time selection
mtime120 = data120.Time_S(idx_e);
range120 = data120.Depth_mean(idx_e);
Sv120 = data120.sv(idx_e);
NASC120 = data120.NASC(idx_e);

mtime70 = data70.Time_S(idx_e);
range70 = data70.Depth_mean(idx_e);
Sv70 = data70.sv(idx_e);
NASC70 = data70.NASC(idx_e);

mtime38 = data38.Time_S(idx_e);
range38 = data38.Depth_mean(idx_e);
Sv38 = data38.sv(idx_e);
NASC38 = data38.NASC(idx_e);

mtime18 = data18.Time_S(idx_e18);
range18 = data18.Depth_mean(idx_e18);
Sv18 = data18.sv(idx_e18);
NASC18 = data18.NASC(idx_e18);


%% And remove the data from the time the Triaxus was out of the water
idx_gap = find(datenum(mtime120)>=datenum([2023 10 20 23 0 0]) & datenum(mtime120)<=datenum([2023 10 21 3 0 0]));
idx_gap18 = find(datenum(mtime18)>=datenum([2023 10 20 23 0 0]) & datenum(mtime18)<=datenum([2023 10 21 3 0 0]));

Sv38(idx_gap)= NaN;
Sv120(idx_gap) = NaN;
Sv70(idx_gap) = NaN;
Sv18(idx_gap18) = NaN;


%% Convert sv to Sv and filter acoustics by time of the day

Sv120 = 10*log10(Sv120); Sv120(Sv120>-60)=NaN;
Sv70 = 10*log10(Sv70); Sv70(Sv70>-60)=NaN;
Sv38 = 10*log10(Sv38); Sv38(Sv38>-60)=NaN;
Sv18 = 10*log10(Sv18); Sv18(Sv18>-60)=NaN;


dayrefs = datevec(mtime70);
dayrefs18 = datevec(mtime18);

morning = find(dayrefs(:,4)>=20 & dayrefs(:,4)<23);
noonish = find(dayrefs(:,4)>=23 | dayrefs(:,4)<4);
afternoon = find(dayrefs(:,4)>=4 & dayrefs(:,4)<8);

morning18 = find(dayrefs18(:,4)>=20 & dayrefs18(:,4)<23);
noonish18 = find(dayrefs18(:,4)>=23 | dayrefs18(:,4)<4);
afternoon18 = find(dayrefs18(:,4)>=4 & dayrefs18(:,4)<8);


%% Bin the data

zbins = [0:5:300];

for zz=1:length(zbins)-1
idx_morning = find(range70(morning)>=zbins(zz) & range70(morning)<zbins(zz+1));
idx_noonish = find(range70(noonish)>=zbins(zz) & range70(noonish)<zbins(zz+1));
idx_afternoon = find(range70(afternoon)>=zbins(zz) & range70(afternoon)<zbins(zz+1));


idx_morning18 = find(range18(morning18)>=zbins(zz) & range18(morning18)<zbins(zz+1));
idx_noonish18 = find(range18(noonish18)>=zbins(zz) & range18(noonish18)<zbins(zz+1));
idx_afternoon18 = find(range18(afternoon18)>=zbins(zz) & range18(afternoon18)<zbins(zz+1));

prof_morning120(zz) = nanmean(Sv120(morning(idx_morning)));
prof_noonish120(zz) = nanmean(Sv120(noonish(idx_noonish)));
prof_afternoon120(zz) = nanmean(Sv120(afternoon(idx_afternoon)));

prof_morning70(zz) = nanmean(Sv70(morning(idx_morning)));
prof_noonish70(zz) = nanmean(Sv70(noonish(idx_noonish)));
prof_afternoon70(zz) = nanmean(Sv70(afternoon(idx_afternoon)));

prof_morning38(zz) = nanmean(Sv38(morning(idx_morning)));
prof_noonish38(zz) = nanmean(Sv38(noonish(idx_noonish)));
prof_afternoon38(zz) = nanmean(Sv38(afternoon(idx_afternoon)));

prof_morning18(zz) = nanmean(Sv18(morning18(idx_morning18)));
prof_noonish18(zz) = nanmean(Sv18(noonish18(idx_noonish18)));
prof_afternoon18(zz) = nanmean(Sv18(afternoon18(idx_afternoon18)));

prof_morning120_std(zz) = nanstd(Sv120(morning(idx_morning)));
prof_noonish120_std(zz) = nanstd(Sv120(noonish(idx_noonish)));
prof_afternoon120_std(zz) = nanstd(Sv120(afternoon(idx_afternoon)));

prof_morning70_std(zz) = nanstd(Sv70(morning(idx_morning)));
prof_noonish70_std(zz) = nanstd(Sv70(noonish(idx_noonish)));
prof_afternoon70_std(zz) = nanstd(Sv70(afternoon(idx_afternoon)));

prof_morning38_std(zz) = nanstd(Sv38(morning(idx_morning)));
prof_noonish38_std(zz) = nanstd(Sv38(noonish(idx_noonish)));
prof_afternoon38_std(zz) = nanstd(Sv38(afternoon(idx_afternoon)));

prof_morning18_std(zz) = nanstd(Sv18(morning18(idx_morning18)));
prof_noonish18_std(zz) = nanstd(Sv18(noonish18(idx_noonish18)));
prof_afternoon18_std(zz) = nanstd(Sv18(afternoon18(idx_afternoon18)));

prof_morning120_n(zz) = length(Sv120(morning(idx_morning)));
prof_noonish120_n(zz) = length(Sv120(noonish(idx_noonish)));
prof_afternoon120_n(zz) = length(Sv120(afternoon(idx_afternoon)));

prof_morning70_n(zz) = length(Sv70(morning(idx_morning)));
prof_noonish70_n(zz) = length(Sv70(noonish(idx_noonish)));
prof_afternoon70_n(zz) = length(Sv70(afternoon(idx_afternoon)));

prof_morning38_n(zz) = length(Sv38(morning(idx_morning)));
prof_noonish38_n(zz) = length(Sv38(noonish(idx_noonish)));
prof_afternoon38_n(zz) = length(Sv38(afternoon(idx_afternoon)));

prof_morning18_n(zz) = length(Sv18(morning18(idx_morning18)));
prof_noonish18_n(zz) = length(Sv18(noonish18(idx_noonish18)));
prof_afternoon18_n(zz) = length(Sv18(afternoon18(idx_afternoon18)));


end


%% Supplementary materials figure with the binned acoustic response at different frequencies


figure(6),clf,subplot(1,4,1),boundedline(prof_morning120,-zbins(1:end-1),prof_morning120_std./sqrt(prof_morning120_n), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning120,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish120,-zbins(1:end-1),prof_noonish120_std./sqrt(prof_noonish120_n), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish120,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon120,-zbins(1:end-1),prof_afternoon120_std./sqrt(prof_afternoon120_n), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon120,-zbins(1:end-1),'k','LineWidth',2)
title('120 kHz')
ylim([-265 0])

subplot(1,4,2),boundedline(prof_morning70,-zbins(1:end-1),prof_morning70_std./sqrt(prof_morning70_n), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning70,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish70,-zbins(1:end-1),prof_noonish70_std./sqrt(prof_noonish70_n), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish70,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon70,-zbins(1:end-1),prof_afternoon70_std./sqrt(prof_afternoon70_n), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon70,-zbins(1:end-1),'k','LineWidth',2)
title('70 kHz')
ylim([-265 0])

subplot(1,4,3),boundedline(prof_morning38,-zbins(1:end-1),prof_morning38_std./sqrt(prof_morning38_n), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning38,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish38,-zbins(1:end-1),prof_noonish38_std./sqrt(prof_noonish38_n), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish38,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon38,-zbins(1:end-1),prof_afternoon38_std./sqrt(prof_afternoon38_n), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon38,-zbins(1:end-1),'k','LineWidth',2)
title('38 kHz')
ylim([-265 0])

subplot(1,4,4),boundedline(prof_morning18,-zbins(1:end-1),prof_morning18_std./sqrt(prof_morning18_n), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning18,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish18,-zbins(1:end-1),prof_noonish18_std./sqrt(prof_noonish18_n), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish18,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon18,-zbins(1:end-1),prof_afternoon18_std./sqrt(prof_afternoon18_n), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon18,-zbins(1:end-1),'k','LineWidth',2)
title('18 kHz')
ylim([-265 0])

print(figure(6),'-dpng','multi_frequency_EAC.png')


%% Makes figure for the background of the contours with the attenuation

mtime18 = data18.Time_S(idx_e18);
range18 = data18.Depth_mean(idx_e18);
Sv18 = data18.sv(idx_e18);
NASC18 = data18.NASC(idx_e18);
Sv18 = 10*log10(Sv18); Sv18(Sv18>-60)=NaN;

% Grid the acoustic data

[mtime18g,range18g] = meshgrid(datenum(unique(mtime18(isfinite(mtime18)))),unique(range18(isfinite(range18))));

good_data = isfinite(range18) & isfinite(mtime18);

Sv18_gridded = griddata(datenum(mtime18(good_data)),range18((good_data)),Sv18((good_data)),datenum(mtime18g),range18g);

% Load and grid the light attenuation contours
light = load('EAC_PAR_profiles.mat'); % these are saved as the output of an06

section_clc = light.section.*(1-light.cast_mask);
section_clc(section_clc==0) = NaN;
dayrefs = datevec(light.dtms_prof);

idx_s = find(light.dtms_prof>7000);

[light_mtime,light_depth] = meshgrid(datenum(light.dtms_prof(idx_s)),-light.profile(light.qq).depth(1:end-1));

% Plot echogram with overlaid contours

figure(12),clf,h=pcolor(mtime18g,-range18g,Sv18_gridded);
set(h, 'EdgeColor', 'none')
hold on,contour(light_mtime,light_depth,real(light.section(:,idx_s)),[0:0.05:0.2],'w','LineWidth',2)
hold on,scatter(datenum(mtime18(idx_gap)),-range18(idx_gap),20,Sv18(idx_gap),'filled')
set(h, 'EdgeColor', 'none')
ylim([-265 0])
caxis([-85 -65])
cmocean('thermal')
datetick('x','keeplimits')
ylabel('depth (m)')
colorbar

