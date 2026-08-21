%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An07 - Plot complementary measurements from Triaxus deployment in the
% Tasman Sea / East Australian Current
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths and graphical preferences
% Select the paths to all needed scripts - this will likely have to change

addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop/Code/')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\BioSWOT\floats\Useful_packages/gsw_matlab_v3_06_12/')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\BioSWOT\floats\Useful_packages/gsw_matlab_v3_06_12/library/')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')

% Graphical preferences
set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')


%% Reads and compile the data

flist = glob('in2023*.mat'); % This may change depending on your file structure


all_mtimes = [];
all_abundance = [];
all_biomass = [];
all_pressure = [];
all_SMEP = [];
all_NBSS = [];
all_temperature = [];
all_salinity = [];
all_bs = [];
all_cast = [];

for ff=1:length(flist)

    f = load(flist{ff});

    all_abundance = [all_abundance;f.s.Abundance];
    all_biomass = [all_biomass; f.s.Biomass];
    all_mtimes = [all_mtimes; f.s.datenum];
    all_pressure = [all_pressure; f.s.pressure];
    all_SMEP = [all_SMEP; f.s.SMEP];
    all_NBSS = [all_NBSS; f.s.NBSS_Slope];
    all_temperature = [all_temperature; f.s.temperature];
    all_salinity = [all_salinity; f.s.salinity];
    bs_interp = interp1(f.EcoTrip.datenum, f.EcoTrip.Backscatter, f.s.datenum);
    all_bs = [all_bs; bs_interp];
    all_cast = [all_cast; f.s.cast_no];

end
all_density = gsw_rho(all_salinity,all_temperature,all_pressure);
dayrefs = datevec(all_mtimes);

all_biomass = log10(all_biomass);


%% Bin the data and separate upcast and downcast

smooth_depth = movmean(all_pressure,10);

[pks,locs] = findpeaks(smooth_depth);
chng = find(diff(locs));

idx_neg = find(diff(smooth_depth)>0.2);
idx_pos = find(diff(smooth_depth)<-0.1);

morning_up = find(dayrefs(idx_pos,4)>=20 & dayrefs(idx_pos,4)<23);
noonish_up = find(dayrefs(idx_pos,4)>=23 | dayrefs(idx_pos,4)<4);
afternoon_up = find(dayrefs(idx_pos,4)>=4 & dayrefs(idx_pos,4)<8);

morning_down = find(dayrefs(idx_neg,4)>=20 & dayrefs(idx_neg,4)<23);
noonish_down = find(dayrefs(idx_neg,4)>=23 | dayrefs(idx_neg,4)<4);
afternoon_down = find(dayrefs(idx_neg,4)>=4 & dayrefs(idx_neg,4)<8);

zbins = [0:20:300];

for zz=1:length(zbins)-1

    % LOPC
    idx_mrn_up = find(all_pressure(idx_pos(morning_up))>=zbins(zz) & all_pressure(idx_pos(morning_up))<zbins(zz+1));
    idx_non_up = find(all_pressure(idx_pos(noonish_up))>=zbins(zz) & all_pressure(idx_pos(noonish_up))<zbins(zz+1));
    idx_aft_up = find(all_pressure(idx_pos(afternoon_up))>=zbins(zz) & all_pressure(idx_pos(afternoon_up))<zbins(zz+1));

    idx_mrn_down = find(all_pressure(idx_neg(morning_down))>=zbins(zz) & all_pressure(idx_neg(morning_down))<zbins(zz+1));
    idx_non_down = find(all_pressure(idx_neg(noonish_down))>=zbins(zz) & all_pressure(idx_neg(noonish_down))<zbins(zz+1));
    idx_aft_down = find(all_pressure(idx_neg(afternoon_down))>=zbins(zz) & all_pressure(idx_neg(afternoon_down))<zbins(zz+1));
    
    prof_morning_upcast(zz) = nanmean(all_biomass(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast(zz) = nanmean(all_biomass(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast(zz) = nanmean(all_biomass(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast(zz) = nanmean(all_biomass(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast(zz) = nanmean(all_biomass(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast(zz) = nanmean(all_biomass(idx_neg(afternoon_down(idx_aft_down))));

    prof_morning_upcast_std(zz) = nanstd(all_biomass(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast_std(zz) = nanstd(all_biomass(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast_std(zz) = nanstd(all_biomass(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast_std(zz) = nanstd(all_biomass(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast_std(zz) = nanstd(all_biomass(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast_std(zz) = nanstd(all_biomass(idx_neg(afternoon_down(idx_aft_down))));

    prof_morning_upcast_n(zz) = length(all_biomass(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast_n(zz) = length(all_biomass(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast_n(zz) = length(all_biomass(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast_n(zz) = length(all_biomass(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast_n(zz) = length(all_biomass(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast_n(zz) = length(all_biomass(idx_neg(afternoon_down(idx_aft_down))));

    % Density 
    prof_morning_upcast_dens(zz) = nanmean(all_density(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast_dens(zz) = nanmean(all_density(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast_dens(zz) = nanmean(all_density(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast_dens(zz) = nanmean(all_density(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast_dens(zz) = nanmean(all_density(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast_dens(zz) = nanmean(all_density(idx_neg(afternoon_down(idx_aft_down))));

    prof_morning_upcast_std_dens(zz) = nanstd(all_density(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast_std_dens(zz) = nanstd(all_density(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast_std_dens(zz) = nanstd(all_density(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast_std_dens(zz) = nanstd(all_density(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast_std_dens(zz) = nanstd(all_density(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast_std_dens(zz) = nanstd(all_density(idx_neg(afternoon_down(idx_aft_down))));

    prof_morning_upcast_n_dens(zz) = length(all_density(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast_n_dens(zz) = length(all_density(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast_n_dens(zz) = length(all_density(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast_n_dens(zz) = length(all_density(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast_n_dens(zz) = length(all_density(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast_n_dens(zz) = length(all_density(idx_neg(afternoon_down(idx_aft_down))));

        % Density gradient
    prof_morning_upcast_densg(zz) = nanmean(abs(diff(all_density(idx_pos(morning_up(idx_mrn_up))))));
    prof_noonish_upcast_densg(zz) = nanmean(abs(diff(all_density(idx_pos(noonish_up(idx_non_up))))));
    prof_afternoon_upcast_densg(zz) = nanmean(abs(diff(all_density(idx_pos(afternoon_up(idx_aft_up))))));

    prof_morning_downcast_densg(zz) = nanmean(abs(diff(all_density(idx_neg(morning_down(idx_mrn_down))))));
    prof_noonish_downcast_densg(zz) = nanmean(abs(diff(all_density(idx_neg(noonish_down(idx_non_down))))));
    prof_afternoon_downcast_densg(zz) = nanmean(abs(diff(all_density(idx_neg(afternoon_down(idx_aft_down))))));

    prof_morning_upcast_std_densg(zz) = nanstd(abs(diff(all_density(idx_pos(morning_up(idx_mrn_up))))));
    prof_noonish_upcast_std_densg(zz) = nanstd(abs(diff(all_density(idx_pos(noonish_up(idx_non_up))))));
    prof_afternoon_upcast_std_densg(zz) = nanstd(abs(diff(all_density(idx_pos(afternoon_up(idx_aft_up))))));

    prof_morning_downcast_std_densg(zz) = nanstd(abs(diff(all_density(idx_neg(morning_down(idx_mrn_down))))));
    prof_noonish_downcast_std_densg(zz) = nanstd(abs(diff(all_density(idx_neg(noonish_down(idx_non_down))))));
    prof_afternoon_downcast_std_densg(zz) = nanstd(abs(diff(all_density(idx_neg(afternoon_down(idx_aft_down))))));

    prof_morning_upcast_n_densg(zz) = length(all_density(idx_pos(morning_up(idx_mrn_up))));
    prof_noonish_upcast_n_densg(zz) = length(all_density(idx_pos(noonish_up(idx_non_up))));
    prof_afternoon_upcast_n_densg(zz) = length(all_density(idx_pos(afternoon_up(idx_aft_up))));

    prof_morning_downcast_n_densg(zz) = length(all_density(idx_neg(morning_down(idx_mrn_down))));
    prof_noonish_downcast_n_densg(zz) = length(all_density(idx_neg(noonish_down(idx_non_down))));
    prof_afternoon_downcast_n_densg(zz) = length(all_density(idx_neg(afternoon_down(idx_aft_down))));


end



%% Plotting

figure(1),clf,boundedline(prof_morning_downcast,-zbins(1:end-1),prof_morning_downcast_std./sqrt(length(prof_morning_downcast_n)), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning_downcast,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish_downcast,-zbins(1:end-1),prof_noonish_downcast_std./sqrt(length(prof_noonish_downcast_n)), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish_downcast,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon_downcast,-zbins(1:end-1),prof_afternoon_downcast_std./sqrt(length(prof_afternoon_downcast_n)), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon_downcast,-zbins(1:end-1),'k','LineWidth',2)
ylim([-265 0])
%xlim([-0.01 0.25])
xlabel('log(biomass) (mg m^{-3})')
ylabel('depth (m)')


figure(2),clf,boundedline(prof_morning_upcast,-zbins(1:end-1),prof_morning_upcast_std./sqrt(length(prof_morning_upcast_n)), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning_upcast,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish_upcast,-zbins(1:end-1),prof_noonish_upcast_std./sqrt(length(prof_noonish_upcast_n)), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish_upcast,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon_upcast,-zbins(1:end-1),prof_afternoon_upcast_std./sqrt(length(prof_afternoon_upcast_n)), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon_upcast,-zbins(1:end-1),'k','LineWidth',2)
ylim([-265 0])
%xlim([-0.01 0.25])
xlabel('log(LOPC biomass) (mg.m^{-3})')
ylabel('depth (m)')
%print(figure(2),'-dpng','upcast_LOPC.png')

%% Density layers

figure(3),clf,boundedline(prof_morning_downcast_densg,-zbins(1:end-1),prof_morning_downcast_std_densg./sqrt(length(prof_morning_downcast_n_densg)), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning_downcast_densg,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish_downcast_densg,-zbins(1:end-1),prof_noonish_downcast_std_densg./sqrt(length(prof_noonish_downcast_n_dens)), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish_downcast_densg,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon_downcast_densg,-zbins(1:end-1),prof_afternoon_downcast_std_densg./sqrt(length(prof_afternoon_downcast_n_densg)), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon_downcast_densg,-zbins(1:end-1),'k','LineWidth',2)
ylim([-265 0])
%xlim([-0.01 0.25])
%xlabel('log(biomass) (mg.m^{-3})')
ylabel('depth (m)')
xlabel('density gradient (kg m^{-4})')
%print(figure(1),'-dpng','downcast_LOPC.png')


figure(4),clf,boundedline(prof_morning_upcast_densg,-zbins(1:end-1),prof_morning_upcast_std_densg./sqrt(length(prof_morning_upcast_n_densg)), 'orientation', 'horiz','b','alpha')
hold on, plot(prof_morning_upcast_densg,-zbins(1:end-1),'b','LineWidth',2)
boundedline(prof_noonish_upcast_densg,-zbins(1:end-1),prof_noonish_upcast_std_densg./sqrt(length(prof_noonish_upcast_n_densg)), 'orientation', 'horiz','r','alpha')
hold on, plot(prof_noonish_upcast_densg,-zbins(1:end-1),'r','LineWidth',2)
boundedline(prof_afternoon_upcast_densg,-zbins(1:end-1),prof_afternoon_upcast_std_densg./sqrt(length(prof_afternoon_upcast_n_densg)), 'orientation', 'horiz','k','alpha')
hold on, plot(prof_afternoon_upcast_densg,-zbins(1:end-1),'k','LineWidth',2)
ylim([-265 0])
%xlim([-0.01 0.25])
%xlabel('log(LOPC biomass) (mg.m^{-3})')
ylabel('depth (m)')
