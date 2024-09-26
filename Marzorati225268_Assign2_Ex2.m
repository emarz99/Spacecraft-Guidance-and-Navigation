clc; clearvars; close all
format long g

% SGN Assignment 2
% Emanuele Marzorati, 225268, 10724126
% Exercise 2
% Load kernels
cspice_furnsh('assignment02.tm');
addpath('sgp4');
%%
% Initial Data

% Time Epoch
et0 = cspice_str2et('2010-08-12T05:30:00.000'); % [s]
etf = cspice_str2et('2010-08-12T11:00:00.000'); % [s]

% Initial State of Satellite 1
%rr0 = [4622.232026629, 5399.3369588058, -0.0212138165769957]; % [km]
%rv0 = [0.812221125483763, -0.721512914578826, 7.42665302729053]; % [km/s] 

% Initial State of Satellite 2 TUNABLE FOR POINT 2.4
% EXERCISE 2.2.A,  EXERCISE 2.3, TLE TO BE SWITCHED TO SECOND SATELLITE
%%
rr0 = [4621.69343340281, 5399.26386352847, -3.09039248714313]; % [km]
rv0 = [0.813960847513811,  -0.719449862738607, 7.42706066911294]; % [km/s] 

x0 = [rr0, rv0]';

% Time of separation
et_sep = cspice_str2et('2010-08-12 05:27:39.114 UTC'); % [s]

% Earth Gravitational Constant
GM = cspice_bodvrd('Earth', 'GM', 1); % [km^3/s^2]

% Considered Stations
stationName1 = 'KOUROU';
stationName2 = 'SVALBARD';

% Minimum Elevations
min_el_1 = 10*cspice_rpd(); % [rad]
min_el_2 = 5*cspice_rpd();  % [rad]

% Not considering J2 effect
tag_j2 = false;

% Noise matrix
noise_m_1 = diag([(0.01)^2, (0.1*cspice_rpd())^2, (0.1*cspice_rpd())^2]); % [km,rad,rad]^2
noise_m_2 = diag([(0.01)^2, (0.125*cspice_rpd())^2, (0.125*cspice_rpd())^2]); % [km,rad,rad]^2
%%
% EXERCISE 2.1

% Creating Time Grid
npoints = round((etf-et0)/60.0)+1; 
et_vec = linspace(et0, etf, npoints);
dt = et_vec(2)-et_vec(1);

% Propagating states on the grid
rv = zeros(6, npoints); 

% Propagation
for i=1:npoints
    % Time constraint
    tf = et_vec(i);
    % Propagation
    [xf,~,~,~] = keplerian_propagator(et_sep, x0, tf, GM, tag_j2);
    % Storing
    rv(:,i) = xf; 
end

% Computing the angles of Satellite 1 wrt the stations
[range1, azimuth1, elevation1] = antenna_pointing(stationName1, et_vec, rv);
[range2, azimuth2, elevation2] = antenna_pointing(stationName2, et_vec, rv);

% Creating logical array of the visibility windows
i_visibility1 = elevation1 >= min_el_1; 
i_visibility2 = elevation2 >= min_el_2;
 
% Compute indices of visible measurements 
sat1_s1_vw = find(i_visibility1);
sat1_s2_vw = find(i_visibility2);
% Compute indices of not visible measurements 
vw_station1 = find(~i_visibility1);
vw_station2 = find(~i_visibility2);

% Get the predicted azimuth and elevation during the visibility windows
% Kourou
azimuth1_copy = azimuth1; % creating a copy of the array
azimuth1_copy(vw_station1) = 0; % turn into 0 the values of non-visibility
elevation1_copy = elevation1; % duplicate the array
elevation1_copy(vw_station1) = 0; % turn into 0 the values of non-visibility
% Svalbard
azimuth2_copy = azimuth2;
azimuth2_copy(vw_station2) = 0;
elevation2_copy = elevation2;
elevation2_copy(vw_station2) = 0; 

% Elevation in Svalbard [deg]
figure()
hold on
plot(et_vec(i_visibility2)/cspice_spd(),elevation2(i_visibility2)*cspice_dpr(), 'c*', LineWidth=2)
plot(et_vec/cspice_spd(),elevation2*cspice_dpr(), 'm', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('Elevation [deg]')
legend(["Visibility","Elevation"])
title({'Svalbard, Elevation in time epoch'})
grid on

% Elevation in Kourou [deg]
figure()
hold on
plot(et_vec(i_visibility1)/cspice_spd(),elevation1(i_visibility1)*cspice_dpr(), 'c*', LineWidth=2)
plot(et_vec/cspice_spd(),elevation1*cspice_dpr(), 'm', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('Elevation [deg]')
legend(["Visibility","Elevation"])
title({'Kourou, Elevation in time epoch'})
grid on

% Azimuth in Svalbard [deg]
figure()
hold on
plot(et_vec(i_visibility2)/cspice_spd(),azimuth2(i_visibility2)*cspice_dpr(), 'c*', LineWidth=2)
plot(et_vec/cspice_spd(),azimuth2*cspice_dpr(), 'm', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('Azimuth [deg]')
legend(["Visibility","Azimuth"])
title({'Svalbard, Azimuth in time epoch'})
grid on

% Azimuth in Kourou [deg] 
figure()
hold on
plot(et_vec(i_visibility1)/cspice_spd(),azimuth1(i_visibility1)*cspice_dpr(), 'c*', LineWidth=2)
plot(et_vec/cspice_spd(),azimuth1*cspice_dpr(), 'm', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('Azimuth [deg]')
legend(["Visibility","Azimuth"])
title({'Azimuth, Visibility in time epoch'})
grid on
%%
% EXERCISE 2.2.A

typerun    = 'u';  
opsmode    = 'a';  
whichconst =  72; 

% TLE 1 
%longstr1 = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
%longstr2 = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';

% TLE 2 %% TO BE TUNED FOR EXERCISE 2.4 WITH ALSO 
longstr1 = '1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996';
longstr2 = '2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55';

satrec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

% Get TLE epoch 
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s', sat_epoch_str);

% SET TIME SPAN
et0 = cspice_str2et('2010-08-12T05:30:00.000'); 
etf = cspice_str2et('2010-08-12T11:00:00.000');
npoints = round((etf-et0)/60.0)+1;
et_vec = linspace(et0, etf, npoints);

% SET NUTATION PARAMETERS
arcsec2rad = pi/(180*3600);
ddpsi = -0.073296*arcsec2rad; %  [rad] 
ddeps = -0.009373*arcsec2rad; %  [rad]

% initialize ECI state vector
rv_eci = zeros(6, npoints);

for i = 1:npoints
    % SGP4 propagation
    tsince_old = (et_vec(i) - sat_epoch_et)/60.0;
    [~,rteme,vteme] = sgp4(satrec,  tsince_old);
    % Time to convert properly
    ttt = cspice_unitim(et_vec(i), 'ET', 'TDT')/cspice_jyear()/100;
    % Conversion
    [rv_eci(1:3,i), rv_eci(4:6,i)] = teme2eci(rteme, vteme, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
end

% EXERCISE 2.2.B

% Using previous visibility windows from 2.1
[range1_spg4_n, azimuth1_spg4_n, elevation1_spg4_n, sat_vis_et1] = measurement_sim(satrec, sat_epoch_et, et_vec, stationName1, min_el_1, noise_m_1, i_visibility1);
[range2_spg4_n, azimuth2_spg4_n, elevation2_spg4_n, sat_vis_et2] = measurement_sim(satrec, sat_epoch_et, et_vec, stationName2, min_el_2, noise_m_2, i_visibility2);
%%
% EXERCISE 2.3

% Setting reference epoch to initial epoch
et_ref = et0;

% Real reference state for Satellite 1, Mango,
%spg4_state = sgp4_propagation(et_ref, 1);
% Real reference state for Satellite 2, Tango TO BE TUNED FOR EXERCISE 2.4 
spg4_state = sgp4_propagation(et_ref, 2);

% Reference state first guess (at TLE epoch)
state_first_guess = spg4_state;

% Real measurements from Station 1
meas_real1 = [range1_spg4_n, azimuth1_spg4_n, elevation1_spg4_n];
% Real measurements from Station 2
meas_real2 = [range2_spg4_n, azimuth2_spg4_n, elevation2_spg4_n];

% Weight for measurements Station 1
W_m_1 = inv(noise_m_1);
% Weight for measurements Station 2
W_m_2 = inv(noise_m_2);
%%
% EXERCISE 2.3.A
% Considering Kuorou ground station, not considering J2 acceleration
tag_j2 = false; % No J2

% Create the cost function for the first ground station
residual1 = @(x) costfunction_ONE_GS(x, et0, sat_vis_et1, meas_real1, W_m_1, GM, 'KOUROU', tag_j2);

% Solve the non-linear least squares
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[state_sat1_OGS, resnorm_sat1_OGS, residual_sat1_OGS, exitflag_sat1_OGS, ~, ~, jac_sat1_OGS] = lsqnonlin(residual1, state_first_guess, [], [], opt);

% Error
(spg4_state - state_sat1_OGS).'

% Computing covariance
Jac_OGS = full(jac_sat1_OGS);
P_sat1_OGS = resnorm_sat1_OGS/(length(residual_sat1_OGS)-length(state_sat1_OGS)).* inv(Jac_OGS'*Jac_OGS)

%%
% EXERCISE 2.3.B 
% Considering Kuorou and Svalbard ground stations
tag_j2 = false; % No J2

% Create the cost function
residual2 = @(x)costfunction_TWO_GS(x, et0, sat_vis_et1, sat_vis_et2, meas_real1, meas_real2, 'KOUROU', 'SVALBARD', W_m_1, W_m_2, GM, tag_j2)

% Solve the non-linear least squares
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[state_sat1_TGS, resnorm_sat1_TGS, residual_sat1_TGS, ~, ~, ~, jac_sat1_TGS] = lsqnonlin(residual2, state_first_guess, [], [], opt);

% Error
(spg4_state - state_sat1_TGS).'

% Computing covariance
Jac_TGS = full(jac_sat1_TGS);
P_sat1_TGS = resnorm_sat1_TGS / (length(residual_sat1_TGS)-length(state_sat1_TGS)).*inv(Jac_TGS'*Jac_TGS)
%%
% EXERCISE 2.3.C 
% Considering Kuorou and Svalbard ground stations, Accounting for J2

% Considering J2 effect
tag_j2 = true;

% Create the cost function
residual3 = @(x)costfunction_TWO_GS(x, et0, sat_vis_et1, sat_vis_et2, meas_real1, meas_real2, 'KOUROU', 'SVALBARD', W_m_1, W_m_2, GM, tag_j2);

% Solve the non-linear least squares
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[state_sat1_j2, resnorm_sat1_j2, residual_sat1_j2, exitflag_sat1_j2, ~, ~, jac_j2] = lsqnonlin(residual3, state_first_guess, [], [], opt);

% Error
(spg4_state - state_sat1_j2).'

% Compute resulting covariance
Jac_sat1_j2 = full(jac_j2);
P_sat1_j2 = resnorm_sat1_j2 / (length(residual_sat1_j2)-length(state_sat1_j2)) .* inv(Jac_sat1_j2'*Jac_sat1_j2);
%%
% clean kernel pool
cspice_kclear();
%% FUNCTIONS

function residual = costfunction_ONE_GS(x, et0, t_r_meas, r_meas, W_m, mu, stationName, j2)

% Initializing residuals
residual = zeros(size(r_meas));

% Propagating guess state in the time grid
for i = 1:length(t_r_meas) 
    [rv(:,i), ~, ~, ~] = keplerian_propagator(et0, x, t_r_meas(i), mu, j2);
end 

% Measuring the guess state in the instants of time of the true measurements
[range, azimuth, elevation] = antenna_pointing(stationName, t_r_meas, rv); 
meas_pred = [range', azimuth', elevation']; % combine into single matrix

% Computation of the residuals 
for i = 1:length(t_r_meas)
    diff_meas_weighted = W_m * (meas_pred(i,:) - r_meas(i,:))';
    residual(i,:) = diff_meas_weighted';
end

end

function residual = costfunction_TWO_GS(x, et0, t_r_meas1, t_r_meas2, r_meas1, r_meas2, stationName1, stationName2, W_m_1, W_m_2, mu, j2)

% Initializing residuals
residual = zeros(size([r_meas1; r_meas2]));

% Propagating guess state in the time grid
% First ground station
for i = 1:length(t_r_meas1) 
    [rv1(:,i), ~, ~, ~] = keplerian_propagator(et0, x, t_r_meas1(i), mu, j2);
end 
% Second ground station
for i = 1:length(t_r_meas2) 
    [rv2(:,i), ~, ~, ~] = keplerian_propagator(et0, x, t_r_meas2(i), mu, j2);
end 

% Measuring the guess state in the instants of time of the true measurements
[range1, azimuth1, elevation1] = antenna_pointing(stationName1, t_r_meas1, rv1); 
[range2, azimuth2, elevation2] = antenna_pointing(stationName2, t_r_meas2, rv2); 

% Combine into single matrix
meas_pred1 = [range1', azimuth1', elevation1']; 
meas_pred2 = [range2', azimuth2', elevation2']; 

% Computation of the residuals 
% First ground station
for i = 1:length(t_r_meas1) 
    diff_meas_weighted = W_m_1 * [meas_pred1(i,1)-r_meas1(i,1), angdiff(meas_pred1(i,2),r_meas1(i,2)), angdiff(meas_pred1(i,3),r_meas1(i,3))]';
    residual(i,:) = diff_meas_weighted';
end
% Second ground station
for i = 1:length(t_r_meas2) 
    diff_meas_weighted = W_m_2 * [meas_pred2(i,1)-r_meas2(i,1), angdiff(meas_pred2(i,2),r_meas2(i,2)), angdiff(meas_pred2(i,3),r_meas2(i,3))]';
    residual(length(t_r_meas1)+i,:) = diff_meas_weighted';
end

end


function [range, azimuth, elevation] = antenna_pointing(stationName, et_vec, xx)

% Define station name
    topoFrame = [stationName, '_TOPO'];

for i=1:length(et_vec)
    
    % Transformation from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et_vec(i));

    % Compute spacecraft position in ECI
    rv_station_eci = cspice_spkezr(stationName, et_vec(i), 'J2000', 'NONE', 'EARTH');

    % Compute station-satellite vector in ECI
    rv_station_sat_eci = xx(:,i) - rv_station_eci;

    % Convert state into topocentric frame
    rv_station_sat_topo = ROT_ECI2TOPO*rv_station_sat_eci;

    [range(i),azimuth(i),elevation(i)] = cspice_reclat(rv_station_sat_topo(1:3));

end

end

function [xf, tf, xx, tt]  = keplerian_propagator(t0, x0, tf, mu, J2)
    
% Perform integration
options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,mu,J2), [t0 tf], x0, options);

% Extract state vector and State Transition Matrix
xf = xx(end,1:6)';
tf = tt(end);

end

% Simulation of the measurements from the ground stations
function [sim_meas_r, sim_meas_az, sim_meas_el, t_vec] = measurement_sim(sat_rec, sat_t, window, station, el_min, cov_station, sat_vis)

% Constant for arcseconds to radians conversions
arcsec2rad = pi/(180*3600);

% Array of a priori visibility ephemeris times
t_vec = window(sat_vis);

% Visibility states in TEME
for i=1:length(t_vec)
    t_since = (t_vec(i)-sat_t)/60; % [min]
    [~, rrteme_vw(:,i), vvteme_vw(:,i)] = sgp4(sat_rec,t_since);
end

% Reference frame rotation, parameters for nutation correction
ddpsi = -0.073296*arcsec2rad; % [rad] 
ddeps = -0.009373*arcsec2rad; % [rad] 
for i=1:length(t_vec)
    ateme = [0;0;0];
    ttt = cspice_unitim(t_vec(i), 'ET', 'TDT')/cspice_jyear()/100; 
    [r_eci(:,i), v_eci(:,i), ~] = teme2eci(rrteme_vw(:,i), vvteme_vw(:,i), ateme, ttt, ddpsi, ddeps);
end

% Derive the associated measurements (without noise) from the ground station
[range, azimuth, elevation] = antenna_pointing(station, t_vec, [r_eci;v_eci]); % [km,rad,rad]

% Update the measurements adding noise (gaussian random distribution)
n = mvnrnd([range',azimuth',elevation'], cov_station);
r_noise = n(:,1);
az_noise = n(:,2);
el_noise = n(:,3);

% Filter out the measurements that do not fulfill the visibility condition
k = el_noise >= el_min; % [logical array]
sim_meas_r = r_noise(k);
sim_meas_az = az_noise(k);
sim_meas_el = el_noise(k);
t_vec = t_vec(k);

end

function [dxdt] = keplerian_rhs(t, x, mu, J2)

% Initialize right-hand-side
dxdt = zeros(6,1);
    
% Extract position
rr = x(1:3);
    
% Compute square distance and distance
dist = norm(rr);

% Compute the gravitational acceleration using Newton's law
a_grav = - mu*rr/dist^3;
        
% Extract the resultant dxdt[6,1]
dxdt(1:3) = x(4:6);
dxdt(4:6) = a_grav;

% Computing J2 accelleration
if J2 == true
    % Value and conversion parameters 
    j2 = 0.0010826269;
    r_earth = cspice_bodvrd('EARTH','RADII',3); % [km] Earth's radius
    et = t; % [s] COnversion ephemeris time

    % Rotation matrix ECI-ECEF
    ROT_ECI2ECEF = cspice_pxform('J2000','ITRF93',et);
        
    % Position in ECEF
    rr_ECEF = ROT_ECI2ECEF*rr;
    dist_ECEF = norm(rr_ECEF);

    % J2 acceleration in ECEF
    j2_x = 3/2*mu*j2*rr_ECEF(1)/dist_ECEF^3*(r_earth(1)/dist_ECEF)^2*(5*(rr_ECEF(3)/dist_ECEF)^2-1);
    j2_y = 3/2*mu*j2*rr_ECEF(2)/dist_ECEF^3*(r_earth(1)/dist_ECEF)^2*(5*(rr_ECEF(3)/dist_ECEF)^2-1);
    j2_z = 3/2*mu*j2*rr_ECEF(3)/dist_ECEF^3*(r_earth(1)/dist_ECEF)^2*(5*(rr_ECEF(3)/dist_ECEF)^2-3);
    a_j2 = [j2_x;j2_y;j2_z];

    % J2 acceleration in ECI
    ROT_ECEF2ECI = cspice_pxform('ITRF93','J2000',et);
    J2_ECI_acceleration = ROT_ECEF2ECI*a_j2;

    % Add to the keplerian acceleration
    dxdt(4:6) = dxdt(4:6) + J2_ECI_acceleration;

end


end

function [state_sgp4] = sgp4_propagation(et, sat)

whichconst =  72;

% Constant for arcseconds to radians conversions
arcsec2rad = pi/(180*3600);

% Get TLE data
if sat == 1
        [sat_rec, ~, ~] = read_3LE(36599, 'tle\36599.3le', whichconst);
elseif sat == 2
        [sat_rec, ~, ~] = read_3LE(36827, 'tle\36827.3le', whichconst);
end
    
% Get TLE epoch 
[year,mon,day,hr,min,sec] = invjday(sat_rec.jdsatepoch, sat_rec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% Propagate
t_since = (et-sat_epoch_et)/60; 
[~, rteme, vteme] = sgp4(sat_rec, t_since); % TEME Reference Frame
% Nutation Correction Parameters [rad] 
ddpsi = -0.073296*arcsec2rad;  
ddeps = -0.009373*arcsec2rad; 
ateme = [0;0;0]; % Acceleration for conversion
% Precession Correction
ttt = cspice_unitim(et, 'ET', 'TDT')/cspice_jyear()/100; 
[reci, veci, ~] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps); % ECI Reference Frame
state_sgp4 = [reci; veci];

end