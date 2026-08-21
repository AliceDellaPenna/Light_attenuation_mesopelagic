%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An06 - Generate light profiles from the Triaxus deployment in the
% Tasman Sea / East Australian Current
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import packages and select graphical preferences

set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')

% These will need to change
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')


%% Reads the data
% The path to this will need to change
data = 'C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\Light_Mesopelagic\DAB\Data\L0_triaxus\Triaxus_deployments8_9_10_11_EAC.nc';

pres = ncread(data,'pres_all');
mtime = ncread(data,'mtime_all');
PAR = ncread(data,'PAR_all');

lon = ncread(data,'lon_all');
lat = ncread(data,'lat_all');


%% Identify dives

dt0 = datenum([2023 10 20]);
dtf = datenum([2023 10 22]);
numdays = dtf-dt0 +1;

idx = find(mtime>=dt0 & mtime<=dtf);
dtms = mtime(idx);
depth = pres(idx);
light = PAR(idx);

smooth_depth = movmean(depth,50);

[pks,locs] = findpeaks(smooth_depth);
chng = find(diff(locs)>100);

%% Separates eddy core from eddy peripheries

ec_lon = 154.2;
ec_lat = -37.8;

for qq= 1:length(lon)
    dist_eddy_center(qq) = m_lldist([ec_lon,lon(qq)],[ec_lat,lat(qq)]);
end


dist_focus= dist_eddy_center(idx);


%% Generate light attenuation profiles


for dts = 1:length(chng)-1
    % Let's put a filter to ignore the night time profiles

    % We extract the profile
    idx_t = locs(chng(dts)):locs(chng(dts+1));

    % Now we want to separate the upcasts and downcasts
    prof_depth = depth(idx_t);
    prof_light = light(idx_t);
  
   
    profile(dts).dtm = nanmean(datenum(dtms(idx_t)));
    profile(dts).dist = nanmean(dist_focus(idx_t));

    %if max(prof_light>100)

    smooth_depth = movmean(prof_depth,50);

    if nanmean(diff(smooth_depth))>0 % we identify if it's an upcast or downcast

        profile(dts).direction = 1; % downcast
    else
        profile(dts).direction = 0; % upcast
    end
    %figure(7),hold on,plot(nanmean(dtms(idx_t)),profile(dts).direction,'ko')

    % Remove all the light measurements where the value is higher deeper
    % than up the same profile

    % I need to decide if I want to do this

    % Smooths light profiles
    prof_light_fit = movmean(prof_light,100);
   

    % Create average binned profiles
    z_layer = 5;
    dz = 0:z_layer:700;


    for zz = 1:length(dz)-1
        idx_b = find(prof_depth>=dz(zz) & prof_depth<dz(zz+1));
        avg_prof_light(zz) = geomean(prof_light(idx_b));
    end

    prof_bin_smooth = movmean(avg_prof_light,5);


    % Sanity check
    %figure(5),clf,plot(light(idx_t),-depth(idx_t),'mo')
    %figure(5),hold on,plot(avg_prof_up,-dz(1:end-1),'k','linewidth',2)
    %figure(5),hold on,plot(avg_prof_down,-dz(1:end-1),'b','linewidth',2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate kd
    %dL = diff(prof_bin_smooth);
    %kd = (1-10.^(dL./20))./z_layer; 

    % Using another equation given that the PAR measurements are in
    % absolute units

    kd = -diff(log(prof_bin_smooth))./(diff(dz(1:end-1)));

    % Sanity check
    % figure(6),clf,plot(kd_up,-dz(1:end-2),'r')
    % figure(6),hold on,plot(kd_down,-dz(1:end-2),'b')
    % xlim([0 0.1])

    % Smooth the kd profiles
    smooth_kd = movmean(kd,5);
   

    % Sanity check
    % figure(6),hold on,plot(smooth_kd_up,-dz(1:end-2),'r','LineWidth',2)
    % figure(6),hold on,plot(smooth_kd_down,-dz(1:end-2),'b','LineWidth',2)
    % xlim([0 0.1])



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(~isnan(smooth_kd))
    % Stores profiles as data structures
    profile(dts).kd = smooth_kd;


    profile(dts).depth = dz(1:end-1);
    profile(dts).light = prof_bin_smooth;
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Separates particles and micronekton components

    id_clear = find(~isnan(profile(dts).light));
    if ~isempty(id_clear)
        dz = profile(dts).depth(id_clear(1));

        %surface_value = profile(dts).light(id_clear(1))-20*log(1-dz*profile(dts).kd(id_clear(1)));
        surface_value = profile(dts).light(id_clear(1));



        % Because these are shallower data instead of selecting the
        % euphotic zone we pick an isolume - this value does not impact the
        % result much but does change the range of values we get

        id_dark = find(profile(dts).light<surface_value-40);

        % if we can find the limit of the euphotic zone we pick
        % the last observations from there
        if ~isempty(id_dark)
            idx_zeu = id_dark(1)-1;
        else % otherwise - we just assume it's the surface
            idx_zeu = 1;
        end

        % adds a check that is not a NaN - if that's the case we
        % replace it with the first non-NaN value
        if isnan(profile(dts).kd(idx_zeu))
            nnan = find(~isnan(profile(dts).kd));
            idx_zeu = nnan(1);
        end

        %hold on,subplot(1,2,1),plot(profile(qq).kd(idx_zeu),-profile(qq).depth(idx_zeu),'ro','markerfacecolor','r')

        % then we try to identify some local minima in the kd curve
        % below the euphotic zone

        profile(dts).kd_filt = profile(dts).kd;
        profile(dts).kd_filt(1:idx_zeu)=NaN;

        [temp,tmpidx] = findpeaks(-profile(dts).kd_filt);



        % to find the second point to fit the curve
        if ~isempty(tmpidx) % if there is at least a minimum
            last = tmpidx(1);  % we pick the first one below the zeuphotic

        else % otherwise we pick the absolute mimum of kd
            nanni = find(~isnan(profile(dts).kd));
            [~,last] = min(profile(dts).kd);

        end

        % we add a check to make sure the two points are not too close - we
        % want them to be at least 50 or 100 m away from teach other

        first = idx_zeu;
        profile(dts).zeu_idx = idx_zeu;

        if (abs(profile(dts).depth(first)-profile(dts).depth(last))<51 & length(tmpidx)>1)
            last = tmpidx(2);
        end


        z1 = profile(dts).depth(first);
        kdz1 = profile(dts).kd(first);
        z2 = profile(dts).depth(last);
        kdz2 = profile(dts).kd(last);

        %figure(1), hold on, subplot(1,2,1),plot(profile(qq).kd(last),-profile(qq).depth(last),'mo','MarkerFaceColor','m');

        lambda = 1./(z2-z1).*log(kdz1./kdz2);
        profile(dts).retrieved = kdz1.*exp(-(profile(dts).depth-z1).*lambda);

        if profile(dts).retrieved(end)-profile(dts).retrieved(1)>0
            profile(dts).retrieved = NaN*profile(dts).retrieved;
        end
    
    
    profile(dts).kd_mknekton = profile(dts).kd-profile(dts).retrieved(1:end-1);
    profile(dts).kd_mknekton(1:last) = NaN;
    end

    
    dayref = datevec(profile(dts).dtm);
    

    if ~isempty(profile(dts).kd_mknekton)
        % Sanity check plots -uncomment to see all profiles
    % figure(8),clf,subplot(1,1,1),plot(profile(dts).kd,-profile(dts).depth(1:end-1),'b')
    % hold on,plot(profile(dts).kd_mknekton,-profile(dts).depth(1:end-1),'k')
    % hold on,plot(profile(dts).retrieved,-profile(dts).depth,'r')
    % xlim([-0.05 0.2])
    % title(sprintf('%02d,%02d:%02d',profile(dts).direction,dayref(4),dayref(5)))
    % 
    % 
    % xlim([-0.05 0.2])
    % print(figure(8),'-dpng',sprintf('plts/EAC_%02d.png',dts))
    end
    end
    %end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We make some plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for qq=1:size(profile,2)
    if ~isempty(profile(qq).kd_mknekton) 
light_section(:,qq) = profile(qq).light;
kd_section(:,qq) = profile(qq).kd;
section(:,qq) = profile(qq).kd_mknekton;
dtms_prof(qq) = profile(qq).dtm;
cast_mask(qq) = profile(qq).direction;
dist_mask(qq) = profile(qq).dist;
    end
end

night = find(dtms_prof<10);
section(:,night)=NaN;

section_clc = section.*(1-cast_mask);
section_clc(section_clc==0) = NaN;
dayrefs = datevec(dtms_prof);

%% Times are in UTC

morning = find(dayrefs(:,4)>=20 & dayrefs(:,4)<23);
noonish = find(dayrefs(:,4)>=23 | dayrefs(:,4)<4);
afternoon = find(dayrefs(:,4)>=4 & dayrefs(:,4)<8);


section_clc(section_clc==0) = NaN;
dayrefs = datevec(dtms_prof);


%% Plotting for paper figure

section_clc(section_clc<0) = 0;

figure(1),clf,boundedline(nanmean(section_clc(:,morning),2),-double(profile(dts).depth(1:end-1)),nanstd(section_clc(:,morning),[],2)./sqrt(length(morning)), 'orientation', 'horiz','b','alpha')
hold on, plot(nanmean(section_clc(:,morning),2),-profile(dts).depth(1:end-1),'b','LineWidth',2)
boundedline(nanmean(section_clc(:,noonish),2),-double(profile(dts).depth(1:end-1)),nanstd(section_clc(:,noonish),[],2)./sqrt(length(noonish)), 'orientation', 'horiz','r','alpha')
hold on,plot(nanmean(real(section_clc(:,noonish)),2),-profile(dts).depth(1:end-1),'r','LineWidth',2)
boundedline(nanmean(section_clc(:,afternoon),2),-double(profile(dts).depth(1:end-1)),nanstd(section_clc(:,afternoon),[],2)./sqrt(length(afternoon)), 'orientation', 'horiz','k','alpha')
hold on,plot(nanmean(section_clc(:,afternoon),2),-profile(dts).depth(1:end-1),'k','LineWidth',2)
ylim([-265 0])
xlim([-0.01 0.25])
xlabel('Kd_{micronekton} (m^{-1})')
ylabel('depth (m)')





