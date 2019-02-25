%% Loic VIENS
% 21/10/2018
% To compute the long-period ground motions of the 2004 off the Kii peninsula earthquake using impulse response functions (IRFs) extracted from the ambient seismic field.
% This example computes the vertical to vertical long-period ground motions of the Mw 7.2 earthquake at two Hi-net stations (KNHH and NAGH).
% A simple (slightly modified) elliptical source model  based on the source inversion studies by Yagi (2004) and Okuwaki and Yagi (2018) is used.
% For each subfault, a triangle function is used for the moment rate function.
% The rupture velocity is 3.1 km/s (can be changed).
% All the waveforms have been filtered between 4 and 10 s using a 2-pass Butterworth filter.
% The results are the same as the Figure 7 of Viens and Denolle (Submitted to BSSA)

% Requirements:
%     - Amplitude calibrated IRFs (in the "GFs" folder)
%     - Earthquake waveforms (in the "Earthquake" folder)
%     - Metadata for the stations (location)
%Functions: func_readsac.m (to read sac files)
%           RS_spectra.m (to compute response spectra)

% Output: simulated waveforms
%%
clear all
close all
clc

%% Folders
dir_ini = '/Users/loic/Documents/DONET_simuls/Example_Github';
dir_eq = [dir_ini  '/Earthquake/']; % Earthquake data
dir_GF = [dir_ini '/GFs/' ]; % IRFs

%% Load station locations
load([dir_ini '/Meta_data.mat'])
stations_names = hinet_data.station ;
Latitude = hinet_data.Latitude ;
Longitude =  hinet_data.Longitude ;

%% Filter parameters
perio2 = 4; % filtered data period 1 (in second)
perio1 = 10; % filtered data period 2 (in second)

%% GFs
station_source = {'M.KMD14'} ; % Virtual source
station = { 'N.KNHH', 'N.NAGH'} ; % Receiver
compo1 = 'U'; % Vertical component
compo2 = 'U'; % Vertical component
limstack = {'06_09'}; % GFs are computed from the ambient seismic field recorded between June 1st and September 30, 2015
V_surf = 3.2; % Surface wave velocity used to time shift the IRFs (in km/s)

%% Length of the simulation, sampling rate,...
delta = 10 ; % Sampling rate in Hz
dt = 1/delta; % dt in s
time1 = 250; % duration of the simulation in s
tps_200s=time1*delta;  % Number of data point considered
t2=0:1/delta:time1-1/delta; % time vector

%% Source model
Eq1.Mw = '7.2'; % Mw
Eq1.M0 =  7.54e+26 ; % M0 in Nm
Eq1.date = '2004_249'; % Date of the earthquake
Eq1.eq_lat= 33.033 ; % Earthquake latitude
Eq1.eq_lon=	136.789 ; % Earthquaake longitude
Eq1.eq_depthsurf = 7; % Depth of the top of the fault plane
dip_angle = 40 ; % Dip angle (in Degree)
strike = 280 ; % Strike angle (in Degree)

Rupt_velocity = 3.1; % Rupture velocity
para_source =[.5 , .82, .18]  ; % Location of the epicenter
halfdurtriangle = 1.6 * delta/2; % Half duration of the source time function (1.6/ 2 = .8 s)
L_l = 54; % length of the fault plane
W_l = 38; % Width of the fault plane

L_s= 2 ; % Length of a subfault (in km)
W_s= 2 ; % Width of a subfault (in km)

N_L = L_l/L_s ; % Number of subfaults in the strike direction
N_W = W_l/W_s ; % Number of subfaults in the dip direction
disp([' - Mw '  Eq1.Mw ': Length: '  num2str(L_l) ', Width: ' num2str(W_l)])
disp([' - Number of subfaults: ' num2str(N_L) ' x ' num2str(N_W)])

%% Get the location of the fault plane limits
[latest,lonestn] = reckon( Eq1.eq_lat,Eq1.eq_lon, km2deg( cos(deg2rad(dip_angle))*W_l*para_source(1)),strike-90); % to get the point at the surface
[lat_north_est,lon_north_est] = reckon(latest,lonestn, km2deg(L_l*para_source(2)),strike-180);
[lat_south_est,lon_south_est] = reckon(latest,lonestn, km2deg(L_l*para_source(3)),strike);
[lat_north_west,lon_north_west] = reckon(lat_north_est,lon_north_est, km2deg( cos(deg2rad(dip_angle))*W_l),strike+90);
[lat_south_west,lon_south_west] = reckon(lat_south_est,lon_south_est, km2deg( cos(deg2rad(dip_angle))*W_l),strike+90);

[dist,~] = distance(lat_north_est,lon_north_est,lat_north_west,lon_north_west);
dist2=deg2km(dist);

%% Fault plan (subfault positions)

[dist,AZ] = distance(lat_south_est,lon_south_est,lat_north_est,lon_north_est);
dist2=deg2km(dist);

size4_7 = 0;
for i = 1:1:N_L
    
    clear latout; clear lonout ; clear latout2 ; clear lonout2
    commpt3 = 1;
    size4_7b =  W_s;
    [latout,lonout] = reckon(lat_south_est, lon_south_est, km2deg(size4_7), AZ);
    
    for j = 1:N_W-1
        [latout2(j),lonout2(j)] = reckon(latout, lonout, km2deg(size4_7b*cosd(dip_angle)), AZ+270);
        size4_7b = W_s+ size4_7b;
    end
    finlat(1,i) = latout;
    finlat(2:N_W,i) = latout2';
    finlon(1,i) =   lonout;
    finlon(2:N_W,i) = lonout2';
    size4_7 = L_s +size4_7 ;
end

%% Depth of each subfault by taking the dip angle
for i = 1:length(finlon(:,1))
    [dist,~] = distance(finlat(1,1),finlon(1,1),finlat(i,1),finlon(i,1));
    dist2=deg2km(dist) * tand(dip_angle);
    Prof_s(i,1:length(finlon(1,:))) =  -(dist2 + Eq1.eq_depthsurf);
end



%% Subfault-epicenter distance
for j = 1:N_W %ligne
    for i = 1:N_L %col
        dd_true_epi1(j,i) = deg2km(distance( finlat(j,i),finlon(j,i) ,  Eq1.eq_lat , Eq1.eq_lon));
    end
end
minMatrix = min(dd_true_epi1(:));
[Pos_w,Pos_l] = find(dd_true_epi1==minMatrix);
Epi(1) = finlon(Pos_w,Pos_l);
Epi(2) = finlat(Pos_w,Pos_l);
disp(['Epicenter location (W, L): ' num2str(Pos_w) ', ' num2str(Pos_l) '; distance to true epicenter: ' num2str(round(minMatrix,2)) ' km' ])

%% Compute distance epicenter other point source
for j = 1:N_W %ligne
    for i = 1:N_L %col
        dist_epi(j,i) = deg2km(distance(Epi(1,2),Epi(1,1),finlat(j,i),finlon(j,i)) )   ;
    end
end


%% Set the source model
target = Eq1.M0/10^7; % M0
val_stf = 2*10^17: .1*10^16 : 5*10^17;

[x, y, z] = ellipsoid(30,0,0,N_L*2/2,N_W/2,1); % Create the Ellipse
xf = x(end/2+1:end,:);
yf = y(end/2+1:end,:);
zf = z(end/2+1:end,:);

xne = xf.^1.9/70 -14; % Modify the ellipse here

x = -round(N_L/2):round(N_L/2);
y = -round(N_W/2):round(N_W/2);
[X,Y] = meshgrid(x,y);

%% The zd matrix contains the normalized amplitude of the source time function for each subfault.
zd = griddata(xne,yf,zf, X , Y,'linear');
zd(:,1) =[];
zd(1,:) =[];
zd(:,end) =[];
zd(end,:) =[];

%% Create the source time function for each subfault and find the true amplitude of the source time (triangle) functions for every subfault using the normalized amplitude matrix "zd"
% so the total seismic moment is  equal to the "target" variable

for tto = 1:length(val_stf)
    stf1test = triangularPulse(0, halfdurtriangle, halfdurtriangle*2, 0:1: 50)* val_stf(tto);
    diraftot = zeros(1,1000);
    for j = 1:N_W %ligne
        for i = 1:N_L %col
            if ~isnan(zd(j,i))
                if dist_epi(j,i)> 0.1 % To take the points which are not the hypocenter (Because the hypocenter is not exactly 0 due to some approximations (but very close to zero))
                    dira = zeros(1,1000);
                    beg(j,i) = round(dist_epi(j,i)/Rupt_velocity(round(length(Rupt_velocity)/2))*delta); % Beginin of the dirac function by taking into account the rupture velocity
                    dira(beg(j,i):beg(j,i)+length(stf1test)-1)=stf1test*zd(j,i);
                    diraftot = dira +diraftot;
                else             % To start with the hypocenter
                    dira = zeros(1,1000);
                    beg(j,i) = 1;
                    dira(beg(j,i):beg(j,i)+length(stf1test)-1)=stf1test*zd(j,i);
                    diraftot = dira +diraftot;
                end
            end
        end
    end
    valf(tto) = sum(diraftot)*dt ;
end
[a,val_b] = min(abs(valf - target));
valf(val_b)
stf1 = triangularPulse(0, halfdurtriangle, halfdurtriangle*2, 0:1: 50)* val_stf(val_b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         Compute waveforms                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computation of the large event
for iii = 1: length(station) % Loop on the number of receiver stations
    %% Load IRF data
    file = [ dir_GF   deblank(station_source{1}) '_' deblank(station{iii}) '_' compo1 '_' ...
        compo2  '_' num2str(perio2) '_' num2str(perio1) 's_deconv_06_09_10Hz.sac' ] ;
    datt=func_readsac(file);
    dat=datt.trace;
    
    %% Get the receiver station coordinates
    rec_oo=find(strcmp(deblank(station{iii}),stations_names)==1);
    
    %% Get the virtual source coordinates
    sour_oo=find(strcmp(station_source{1},stations_names)==1);
    %% Virtual source receiver distance
    [dist_sta_sta(iii), ~] = distance(Latitude(sour_oo),Longitude(sour_oo),Latitude(rec_oo),Longitude(rec_oo));
    rk = deg2km(dist_sta_sta(iii));
    %% Epicenter-receiver distance
    [Epicentral_dist(iii)] = deg2km(distance(Epi(1,2),Epi(1,1),Latitude(rec_oo),Longitude(rec_oo))) ;
    %% Compute distances: hypocenter- other subfault centers and difference dist hypocenter- Receiver and epicenter-receiver
    for j = 1:N_W
        for i = 1:N_L
            %% Epicenter to subfault distances
            dist_epi(j,i) = deg2km(distance(Epi(1,2),Epi(1,1),finlat(j,i),finlon(j,i)))    ;
            %% distance station and each point source + difference dist hypocenter- Receiver and epicenter-receiver
            [dist_4(j,i),~] = distance(finlat(j,i),finlon(j,i),Latitude(rec_oo),Longitude(rec_oo)) ;
            r_over_rij(j,i) =sqrt(deg2km(dist_sta_sta(iii)))/ sqrt(deg2km(dist_4(j,i))) ;
        end
    end
    
    %% Taper the GFs with a Hanning function
    t_min = round(deg2km(dist_sta_sta(iii))/6)*delta; %-> 6 km/s, to smooth the data that travel faster than 6 km/s
    hann1=hann(t_min*2);
    dat(1:length(hann1)/2)=hann1(1:end/2).*dat(1:length(hann1)/2);
    
    %% Time shift and amplitude correction of the GFs
    for j = 1:N_W %ligne
        for i = 1:N_L %col
            rij(j,i) = deg2km(distance(Latitude(rec_oo),Longitude(rec_oo),finlat(j,i),finlon(j,i)) );
            shifg(j,i) = round((rij(j,i) -rk)/V_surf*delta); % Time shift GFs with a constant surface wave velocity
            dattmp = zeros(1,600*delta);
            dattmp(shifg(j,i):shifg(j,i)+tps_200s-1) = dat(1:tps_200s) ; % Shift the waveforms
            finshift= dattmp(1:tps_200s)*r_over_rij(j,i); % Geometrical spreading correction
            for pp1 = 1:length(finshift(1,:))
                ZI2(j,i,pp1)  = finshift(pp1);
            end
        end
    end
    
    %% Loop over the rupture velocity
    for yyy = 1 :length(Rupt_velocity)
        final = 0;
        diraftot = zeros(1,50*delta);
        for j = 1:N_W %ligne
            for i = 1:N_L %col
                clear q2
                ZI3 = [];
                for k =1:length(ZI2(1,1,:))
                    ZI3(k) =  ZI2(j,i,k) ;
                end
                
                if ~isnan(zd(j,i))
                    if dist_epi(j,i)> 0.005
                        dira = zeros(1,length(diraftot));
                        beg(j,i) = round(dist_epi(j,i)/Rupt_velocity(yyy)*delta);
                        dira(beg(j,i):beg(j,i)+length(stf1)-1)=stf1*zd(j,i);
                        q2 = conv(ZI3,dira,'full'); % convolution of the STF with the GF
                        final = q2 + final;  % sumation
                        diraftot = dira +diraftot;
                    else
                        dira = zeros(1,length(diraftot));
                        beg(j,i) = 1;
                        dira(beg(j,i):beg(j,i)+length(stf1)-1)=stf1*zd(j,i);
                        q2 = conv(ZI3,dira,'full');
                        final = q2 + final;
                        diraftot = dira +diraftot;
                    end
                end
            end
        end
        
        final_s(iii,:) = final(1:tps_200s);
        
        %% Plot simulated and observed waveforms
        file_eq = [ dir_eq    deblank(station{iii}) '_' compo2  '_' Eq1.date '_Mw_' Eq1.Mw '_10Hz.sac' ] ;
        dat_eq0=func_readsac(file_eq);
        dat_eq(iii,:)=[dat_eq0.trace; 0];
        dat_eqf = dat_eq;
        subplot(2,1,iii)
        plot(t2, final_s(iii,:) , 'r', 'linewidth', 1.5)
        hold on
        grid on
        plot(t2, dat_eq(iii,:), 'k', 'linewidth', 1.5)
        ylabel('Velocity (cm/s)')
        xlabel('Time (s)')
        legend('Simulation', 'Observation')
        title({[station{iii}(3:end) ' station, Epicentral distance: ' num2str(round(Epicentral_dist(iii))) ' km']; ['Filter ' num2str(perio2) '-' num2str(perio1) ' s']})
        
        
    end
    clear ZI2
end

%% Plot  subfault moment rate function maximum amplitude and total STF

t_sou = 0: 1/delta : length(diraftot)/delta - 1/delta;
fig = figure('position' , [440   403   960   395]);
subplot(1,3,1:2)
imagesc(0:L_l, 0:W_l, zd*val_stf(val_b))
xlabel('Distance along strike (km)')
ylabel('Distance along dip (km)')
colormap(flipud((hot)))
c = colorbar('southoutside');
c.Label.String = 'Moment rate function maximum amplitude (Nm/s)';
c.Label.FontSize = 12;
set(gca, 'fontsize', 12)

subplot(1,3,3);
plot(t_sou,diraftot, 'k' ,'linewidth', 2)
grid on
xlim([0 22])
xlabel('Time (s)')
ylabel('Moment rate (Nm/s)')
title({'Total source time function'; 'for the Mw 7.2 event'})
set(gca, 'fontsize', 12)

%% Compute and plot acceleration response spectra for the earthquake and the simulation.

figure;
T = 4:.1:10; % Period vector (in s)
xi = .05 ;% Damping: 0.05 -> 5%
for i =1 :2
    acc_sim =  diff(final_s(i,:)).*1/dt; % Differenciate once to retrieve the acceleration waveforms of the simulations
    acc_eq = diff(dat_eq(i,:)).*1/dt;% Differenciate once to retrieve the acceleration waveforms of the recorded earthquake
    
    % Compute the response spectra using the Duhamel's integral technique.
    [S(i,:)] = RS_spectra(acc_sim, delta, T, xi, 'SA');
    [E(i,:)] = RS_spectra(acc_eq, delta, T, xi, 'SA');
    
    subplot(2,1,i) % Plot
    semilogy(T, S(i,:), 'r', 'linewidth', 1.5)
    hold on
    plot(T, E(i,:), 'k', 'linewidth', 1.5)
    
    ylabel('SA (cm/s^2)')
    grid on
    ylim([.3 5])
    title({[station{i}(3:end) ' station ']})
    set(gca, 'fontsize', 12)

end
legend('Simulation', 'Earthquake')
xlabel('Period (s)')
set(gca, 'fontsize', 12)



