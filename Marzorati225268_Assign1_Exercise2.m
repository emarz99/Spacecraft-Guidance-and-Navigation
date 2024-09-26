clc, clearvars, close all
format long g

% SGN Assignment 1
% Emanuele Marzorati, 225268, 10724126
% Exercise 2

%% INITIALIZE

% Load kernels
cspice_furnsh('ex02.tm');

% Selecting frame and center 
center = 'SSB';
frame = 'ECLIPJ2000'; 
% Define target
object = '20099942';

% Define list of celestial bodies:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initialize propagation data by writing this function
bodies = nbody_init(labels);

%% EXERCISE 2.1

% Define input UTC (coordinated universal time = civil time at Greenwich):
ref_epoch_str = '2029-Jan-01 00:00:00.0000 TDB';
final_epoch_str = '2029-Jul-31 00:00:00.0000 TDB';

% Time conversion 
et0 = cspice_str2et(ref_epoch_str);
etf = cspice_str2et(final_epoch_str);

% Initial state
x0 = cspice_spkezr('20099942',et0,frame,'NONE',center);

% Integration in time
options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
[tt, xx] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), [et0 etf], x0, options);

% Retriving Apophis state
rv_apophis = cspice_spkezr('20099942',tt',frame,'NONE',center)';

%a=rv(:,1:3)
% Computing vector of Apophis wrt the bodies in every instant of time 
rr_sun_obj   = rv_apophis(:,1:3)' - cspice_spkpos('10'  ,tt',frame,'NONE',center) ;
rr_moon_obj  = rv_apophis(:,1:3)' - cspice_spkpos('Moon' ,tt',frame,'NONE',center);
rr_earth_obj = rv_apophis(:,1:3)' - cspice_spkpos('Earth',tt',frame,'NONE',center) ;
%rr_earth_obj = cspice_spkpos(object,tt,frame,'NONE','Earth') ;

% Computing distance
dist_sun_km = sqrt(sum(rr_sun_obj.^2,1));
dist_moon_km = sqrt(sum(rr_moon_obj.^2,1));
dist_earth_km = sqrt(sum(rr_earth_obj.^2,1));

% Converting to AU
dist_sun_AU = cspice_convrt(dist_sun_km,'km','au');
dist_moon_AU = cspice_convrt(dist_moon_km,'km','au');
dist_earth_AU = cspice_convrt(dist_earth_km,'km','au');

% Computing new needed vectors
rr_obj_sun = - rr_sun_obj;
rr_obj_earth = - rr_earth_obj;

%Computing the angle in two different ways
costheta = max(min(dot(rr_obj_sun,rr_obj_earth)/(norm(rr_obj_sun)*norm(rr_obj_earth)),1),-1);
alpha = real(acosd(costheta));
%alpha = rad2deg(atan2(norm(cross(rr_obj_sun,rr_obj_earth)),dot(rr_obj_sun,rr_obj_earth)));

%Plotting
subplot(2,2,1)
plot(tt, dist_sun_AU,'r')
title('Relative distance Apophis-Sun')
xlabel('Epoch [MJD2000]')
ylabel('Distance [AU]')

subplot(2,2,2)
plot(tt, dist_moon_AU,'k')
title('Relative distance Apophis-Moon')
xlabel('Epoch [MJD2000]')
ylabel('Distance [AU]')

subplot(2,2,3)
plot(tt, dist_earth_AU,'b')
title('Relative distance Apophis-Earth')
xlabel('Epoch [MJD2000]')
ylabel('Distance [AU]')

figure
plot(tt,alpha)
title('Angle Earth-Apophis-Sun')
xlabel('Epoch [MJD2000]')
ylabel('Angle [deg]')

% We take the point of smaller distance
[CA_dist,idx] = min(dist_earth_AU);
et_CA = tt(idx);
CA_epoch = cspice_et2utc(et_CA,'C', 1e-3);
fprintf('Close approach epoch [UTC]: %s', CA_epoch);

% Creating new time span
f_initial_window = et_CA - 6*3600;
t_final_window   = et_CA + 6*3600;

tt_window = linspace(f_initial_window,t_final_window,100);

t_initial_window_utc = cspice_et2utc(f_initial_window, 'C', 1e-3);
t_final_window_utc = cspice_et2utc(t_final_window, 'C', 1e-3);

% Position vector of Apophis in every instant of time window wrt Earth
% Baricentre
xx_earth_obj = cspice_spkpos(object,tt_window,frame,'NONE','Earth');
% [radius,lon,lat] = cspice_reclat(xx_earth_obj);
% groundtrack_lon = lon*(180/pi);
% groundtrack_lat = lat*(180/pi);

% Groundtrack plotting
groundtrack_lat = zeros(size(tt_window));
groundtrack_lon = zeros(size(tt_window));

for i=1:length(tt_window)
      et = tt_window(i);

      rot_J2_ECEF = cspice_pxform(frame,'IAU_EARTH',et);
      xx_ecef_obj =  rot_J2_ECEF*xx_earth_obj(:,i);
 
      coord = ecef2lla((xx_ecef_obj').*10^3);
 
      groundtrack_lat(i) = coord(1);
      groundtrack_lon(i) = coord(2);
end

% plotting
figure()
hold on
worldmap('World')
load coastlines
plotm(coastlat,coastlon,'k','Linewidth',1.2)
geoshow('landareas.shp','Facecolor',[0.8 0.7 0.8])
geoshow(groundtrack_lat,groundtrack_lon,'Linewidth',1.3,'Color','r')
xlabel('longitude [deg]')
ylabel('latitude [deg]')
title('Ground Tracks')
hold off

%% EXERCISE 2.2 - 2.3



cspice_kclear();

%% FUNCTIONS 

% OBJECTIVE FUNCTION



% FUNCTION N-BODIES
function [bodies] = nbody_init(labels)

    % Initialize cell array bodies
    bodies = cell(size(labels));
    
    % Loop over labels
    for i=1:size(labels)
        % Store body label
        bodies{i}.name = labels{i};
        % Store body gravitational constant
        bodies{i}.GM = cspice_bodvrd(labels{i},'GM',1);
    end
end

% RHS FOR N-BODY PROBLEM
function [dxdt] = nbody_rhs(t, x, bodies, frame)

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the object position from state x
rr_ssb_obj = x(1:3);

% Loop over all bodies
for i = 1:length(bodies)

    % Retrieve position and velocity of i-th celestial body wrt Solar
    % System Barycentre in inertial frame
    rv_ssb_body = cspice_spkezr(bodies{i}.name, t, frame,'NONE', 'SSB') ;

    % Extract object position wrt. i-th celestial body
    rr_body_obj = rr_ssb_obj - rv_ssb_body(1:3) ;

    % Compute square distance and distance
    dist2 = dot(rr_body_obj, rr_body_obj); 
    dist = sqrt(dist2) ;

    % Compute the gravitational acceleration using Newton's law
    aa_grav = - bodies{i}.GM * rr_body_obj /(dist*dist2);

    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + aa_grav; 
end

end
