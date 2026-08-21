%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An12 - Plot light attenuation profiles and sonar acoustic data from the Kerguelen
% region from equipped southern elelphant seals and vessels of opportunity
% respectively.
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
addpath('/Users/adel386/OneDrive - The University of Auckland/Desktop/Teaching/Courses/MAR399/Chatham_Rise/Matlab_code/libraries/m_map/')
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop/Code/cmocean')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')

%% Read sonar acoustic data

matches = readtable('seals_shipacoustics_match_results_2025-03-13-FILTER.csv');
%magic_ship = 'E:\Acoustic\IMOS\38kHz\IMOS_SOOP-BA_AE_20121215T021630Z_VHLU_FV02_Austral-Leader-II-ES60-38_END-20121222T190951Z_C-20160222T052709Z.nc';
magic_ship = '\IMOS_SOOP-BA_AE_20121215T021630Z_VHLU_FV02_Austral-Leader-II-ES60-38_END-20121222T190951Z_C-20160222T052709Z.nc';
magic_ship_D = 'E:\Acoustic\IMOS\38kHz\IMOS_SOOP-BA_AE_20121215T021630Z_VHLU_FV02_Austral-Leader-II-ES60-38_END-20121222T190951Z_C-20160222T052709Z.nc';
magic_ship_D(1) = 'E';

ship_name_edit = magic_ship;

lon_ship = ncread(ship_name_edit,'LONGITUDE');
lat_ship = ncread(ship_name_edit,'LATITUDE');


Sv = ncread(ship_name_edit,'Sv');
mtime = ncread(ship_name_edit,'TIME')+datenum([1950 1 1]);
depth = ncread(ship_name_edit,'DEPTH');
Sv_2 = 10*log10(Sv);

%% Plot echogram
figure(1),clf,h=pcolor(mtime,-depth,Sv_2);
set(h, 'EdgeColor', 'none')
cb = colorbar;
title(cb,'S_v (dB re m^{-1})')
caxis([-90 -60]) % pick another example
xlim([datenum([2012 12 21 1 0 0]) datenum([2012 12 22 12 0 0])])
cmocean('thermal')
datetick('x','keeplimits')
ylabel('depth (m)')
ylim([-800 0])


%% Load seal data

idx_seals = find(strcmp(matches.ship_fname,magic_ship_D));
seal_indiv = unique(matches.id(idx_seals));
f=load('all_kd_micronekton.mat');


%% We add the metadata information and extract matches
check = 0;
for pp=1:size(f.profiles,2)
    idx_tag = find(strcmp(f.profiles(pp).tag, matches.id));
    if ~isempty(idx_tag)
        
        for rr=1:length(idx_tag)
            if floor(f.profiles(pp).mtime)-datenum(matches.date(idx_tag(rr)))==0
                check=check+1;
                disp(datevec(f.profiles(pp).mtime))
                f.profiles(pp).matching_ship = matches.ship_fname(idx_tag(rr));
                f.profiles(pp).matching_ship_time = matches.ship_date(idx_tag(rr));
                disp('match found!')
            end
        end
    end
end

%% Extract matching SST to classify profiles


sst_map = '20121222120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc';
sst = ncread(sst_map,'analysed_sst')-273.15;
lon_sst = ncread(sst_map,'lon');
lat_sst = ncread(sst_map,'lat');
interp_sst = interp2(lon_sst,lat_sst,sst',lon_ship,lat_ship);

%% We compile the matches into averages by zone

kk=0;
avg_profile = zeros(2,length(f.profiles(1).depth));
sum_profile = zeros(2,length(f.profiles(1).depth));
profiles_az = [];
profiles_pfz = [];
for ll = 1:size(f.profiles,2)
    if ~isempty(f.profiles(ll).matching_ship)
        disp('match is there')
        %disp(f.profiles(ll).matching_ship)
        if(strcmp(f.profiles(ll).matching_ship,magic_ship_D))
            disp('yay')
            kk=kk+1;
            if ~isempty(f.profiles(ll).kd_micronekton)
        
                    if (strcmp(f.profiles(ll).sst_region,'PFZ'))
                        profiles_pfz = [profiles_pfz;f.profiles(ll).kd_micronekton];
                        avg_profile(1,~isnan(f.profiles(ll).kd_micronekton)) = avg_profile(1,~isnan(f.profiles(ll).kd_micronekton))+f.profiles(ll).kd_micronekton(~isnan(f.profiles(ll).kd_micronekton));
                        sum_profile(1,:) = sum_profile(1,:)+~isnan(f.profiles(ll).kd_micronekton);
                    elseif (strcmp(f.profiles(ll).sst_region,'AZ'))
                        avg_profile(2,~isnan(f.profiles(ll).kd_micronekton)) = avg_profile(2,~isnan(f.profiles(ll).kd_micronekton))+f.profiles(ll).kd_micronekton(~isnan(f.profiles(ll).kd_micronekton));
                        sum_profile(2,:) = sum_profile(2,:)+~isnan(f.profiles(ll).kd_micronekton);
                        profiles_az = [profiles_az;f.profiles(ll).kd_micronekton];
                end
            else('bad profile')
            end
        end
    end
end




%% Plot the average profiles in acoustics and light attenuation for the SAZ and the PFZ

figure(10),clf,
subplot(1,2,2)
boundedline(nanmean(profiles_az,1), -double(f.profiles(ll).depth),nanstd(profiles_az,0,1)./sqrt(size(profiles_az,1)), 'orientation', 'horiz','r','alpha')
hold on
plot(nanmean(profiles_az,1),-f.profiles(ll).depth,'r','linewidth',1.5)
boundedline(nanmean(profiles_pfz,1), -double(f.profiles(ll).depth),nanstd(profiles_pfz,0,1)./sqrt(size(profiles_pfz,1)), 'orientation', 'horiz','alpha')
hold on
plot(nanmean(profiles_pfz,1),-f.profiles(ll).depth,'b','linewidth',1.5)
ylim([-600 0])
ylabel('depth (m)')
xlabel ('Kd_{micronekton} (m^{-1})')
xlim([0 0.03])
dn = ncread(ship_name_edit,'day')';

surf_temp = interp_sst;
idx_pfz = find(surf_temp<=8 & surf_temp>=4 & dn'==1);
idx_az = find(surf_temp<4 & dn'==1);
subplot(1,2,1)
standard_error_pfz = std(Sv_2(:,idx_pfz),1,2,'omitnan')./sqrt(length(idx_pfz));
nanni = ~isnan(standard_error_pfz);
boundedline(nanmean(Sv_2(nanni,idx_pfz),2), -double(depth(nanni)),standard_error_pfz(nanni), 'orientation', 'horiz')
hold on,
plot(nanmean(Sv_2(:,idx_pfz),2),-depth,'b','linewidth',1.5)
standard_error_az = std(Sv_2(:,idx_az),1,2,'omitnan')./sqrt(length(idx_az));
nanni = ~isnan(standard_error_az);
boundedline(nanmean(Sv_2(nanni,idx_az),2), -double(depth(nanni)),standard_error_az(nanni), 'orientation', 'horiz','r')
plot(nanmean(Sv_2(:,idx_az),2),-depth,'r','LineWidth',1.5)
ylabel('depth (m)')
xlabel ('S_v (dB re m^{-1})')
ylim([-600 0])
print(figure(10),'-dpng','profiles_together_case_study1.png')
