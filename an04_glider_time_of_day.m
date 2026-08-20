%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An04 - Plots light profiles from the glider equipped with the Wildlife
% Computers tag binned by times of the day to show the movements associated
% with layers of strong light attenuation
% Code by Alice Della Penna (alice.penna@auckland.ac.nz)
% Code associated to submitted manuscript: "What we do in the shadows:
% using light attenuation to observe pelagic ecosystems"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
% Load required packages and custom functions
%% ------------------------------------------------------------------------
% The lines below will need to be changed based on your path
addpath('C:\Users\adel386\OneDrive - The University of Auckland\Desktop\Research\Light_Mesopelagic\')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/Inpaint_nans')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/boundedline')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/catuneven')
addpath('C:/Users/adel386/OneDrive - The University of Auckland/Desktop/Code/boundedlinepackage/singlepatch')

% Sets graphical preferences
set(0,'defaultaxesfontsize',20)
set(0,'defaulttextfontsize',20)
set(0,'defaulttextfontname','helvetica')

%% ------------------------------------------------------------------------
% Read archival light-tag file
%
% Variables:
%   Time       MATLAB timestamps
%   Depth      depth (m)
%   LightLevel raw light measurements
%% ------------------------------------------------------------------------

tb = readtable('2290026_a3-Archive.csv');

dtms  = datenum(tb.Time);
depth = tb.Depth;
light = tb.LightLevel;


%% ------------------------------------------------------------------------
% Read glider surfacing positions
%
% The tag does not contain geolocation continuously (but only when the glider surfaces, 
% so latitude and longitude are interpolated from the glider GPS fixes obtained during
% surfacings.
%% ------------------------------------------------------------------------

track = readtable('glidersurfacings.csv','Delimiter',',');

lon_interp = interp1( ...
    datenum(track.dt,'yyyy-mm-dd HH:MM:SS'), ...
    track.lon_dd, ...
    dtms);

lat_interp = interp1( ...
    datenum(track.dt,'yyyy-mm-dd HH:MM:SS'), ...
    track.lat_dd, ...
    dtms);

%% ------------------------------------------------------------------------
% Restrict analysis to the deployment period
%
% This removes pre-deployment and post-recovery records that are not
% relevant to the biological analysis.
%% ------------------------------------------------------------------------

dt0 = datenum([2025 9 9]);
dtf = datenum([2025 10 4]);

numdays = dtf - dt0 + 1;

idx = find(dtms >= dt0 & dtms <= dtf);

dtms  = dtms(idx);
depth = depth(idx);
light = light(idx);

%% ------------------------------------------------------------------------
% Locate turning points between adjacent profiles
%
% Depth is smoothed to reduce noise before searching for maxima
% corresponding to profile endpoints.
%% ------------------------------------------------------------------------

smooth_depth = movmean(depth,50);

[pks,locs] = findpeaks(smooth_depth);

chng = find(diff(locs) > 100);

%% ------------------------------------------------------------------------
% Loop through all identified profiles
%
% For each profile:
%
%   1. Determine cast direction
%   2. Bin light data vertically
%   3. Estimate Kd
%   4. Smooth Kd profile
%   5. Separate background attenuation and organism signal
%% ------------------------------------------------------------------------
% Runs a for loop over the profiles when data were collected

for dts = 1:length(chng)-1

    % We extract the profile
    idx_t = locs(chng(dts)):locs(chng(dts+1));

    % Now we want to separate the upcasts and downcasts
    prof_depth = depth(idx_t);
    prof_light = light(idx_t);
    
    if max(prof_depth)>=100
    smooth_depth = movmean(prof_depth,50);
    

    if nanmean(diff(smooth_depth))>0.05 % we identify if it's an upcast or downcast

        profile(dts).direction = 1; % downcast
    else
        profile(dts).direction = 0; % upcast
    end

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
    dL = diff(prof_bin_smooth);
    kd = (1-10.^(dL./20))./z_layer; 

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

    % Stores profiles as data structures
    profile(dts).kd = smooth_kd;


    profile(dts).depth = dz(1:end-1);
    profile(dts).light = prof_bin_smooth;
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Separates particles and micronekton components

    id_clear = find(~isnan(profile(dts).light));
    if ~isempty(id_clear)
        dz = profile(dts).depth(id_clear(1));

        surface_value = profile(dts).light(id_clear(1))-20*log(1-dz*profile(dts).kd(id_clear(1)));



        % Definition of euphotic zone using Martini's funny maths;
        % -40 is two orders of magnitude
        % -20 will be for the gliders and floats

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

    profile(dts).dtm = nanmean(datenum(dtms(idx_t)));
    dayref = datevec(profile(dts).dtm);
    

    if ~isempty(profile(dts).kd_mknekton)
    % Uncomment portion below if you want to visualise each profile
    % figure(8),clf,subplot(1,1,1),plot(profile(dts).kd,-profile(dts).depth(1:end-1),'b')
    % hold on,plot(profile(dts).kd_mknekton,-profile(dts).depth(1:end-1),'k')
    % hold on,plot(profile(dts).retrieved,-profile(dts).depth,'r')
    % xlim([-0.05 0.1])
    % title(sprintf('%02d,%02d:%02d',profile(dts).direction,dayref(4),dayref(5)))
    % 
    % 
    % xlim([-0.05 0.1])
    % print(figure(8),'-dpng',sprintf('plts/glider_%02d.png',dts))
    end
    end
end


%% We produce plots classified by time of the day
%
section = NaN*ones(length(profile(dts).depth)-1,size(profile,2));

for qq=1:dts
    if ~isempty(profile(qq).kd_mknekton)
        light_section(:,qq) = profile(qq).light;
        kd_section(:,qq) = profile(qq).kd;
        section(:,qq) = profile(qq).kd_mknekton;
        dtms_prof(qq) = profile(qq).dtm;
        cast_mask(qq) = profile(qq).direction;
    end
end

dayrefs = datevec(dtms_prof); % convert in an easier format to extract hours

night = find(light_section(1,:)<100); % define night as profiles where surface light is below a threshold
section(:,night)=NaN;  % We remove any spurious night profile


section_clc = section.*(1-cast_mask); % we select the casts in the right direction


section_clc(1:25,:) = NaN; % remove 'near surface' data - hard to compare to anything and not what we want to focus on anyways
section_clc(section_clc<0) = 0;
%section_clc(section_clc>1) = NaN;

% We bin by hours to get hourly average profiles

ttz=[0:23];

n = zeros(139,24);
prof_by_time = zeros(139,24);
for qq=1:size(section,2)
    profilez = section_clc(:,qq);
    hh = dayrefs(qq,4)+1;
    
    
    idx = find(~isnan(profilez));
    if ~isempty(idx)
    prof_by_time(idx,hh) = prof_by_time(idx,hh)+profilez(idx);
    n(idx,hh) = n(idx,hh)+1;
    end
    
end

for tt=1:length(ttz)
avg_prof_time(:,tt) = prof_by_time(:,tt)./n(:,tt);
end

% And produce the plots

figure(15),clf,subplot(1,3,1),plot(avg_prof_time(:,19:21),-profile(dts).depth(1:end-1),'linewidth',3)
xlim([-0.001 0.04])
ylim([-700 0])
ylabel('depth(m)')
subplot(1,3,2),plot(avg_prof_time(:,22:24),-profile(dts).depth(1:end-1),'linewidth',3)
hold on,plot(avg_prof_time(:,1:3),-profile(dts).depth(1:end-1),'linewidth',2)
xlim([-0.001 0.04])
ylim([-700 0])
xlabel('Kd_{micronekton} (m^{-1})')
subplot(1,3,3),plot(avg_prof_time(:,4:6),-profile(dts).depth(1:end-1),'linewidth',3)
xlim([-0.001 0.04])
ylim([-700 0])
print(figure(15),'-dpng','hourly_kd.png')





