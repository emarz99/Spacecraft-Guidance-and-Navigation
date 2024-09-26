clc; clearvars; cspice_kclear; close all
format long g
%%

% SGN Assignment 2
% Emanuele Marzorati, 225268, 10724126
% Exercise 3
cspice_furnsh('assignment02.tm');
addpath('sgp4');
%%
% Initial Data

% Ground station
station = 'SVALBARD';

% Minimum elevation of the station
el_min = 5*cspice_rpd(); % [rad]

% Noise matrix of the station
sigma = diag([(0.01)^2, (0.1*cspice_rpd())^2, (0.1*cspice_rpd())^2]); % [km,rad,rad]^2

% Reference epoch
et_sep = cspice_str2et('2010-08-12 05:27:39.114 UTC'); % [s]

% Mean state of Satellite 1 
r0 = [4622.232026629; 5399.3369588058; -0.0212138165769957]; % [km]
v0 = [0.812221125483763; -0.721512914578826; 7.42665302729053]; % [km/s]
rv0 = [r0; v0]; % [km; km/s]

% Covariance at separation
P0 = [5.6e-7, 3.5e-7, -7.1e-8, 0, 0, 0;
      3.5e-7, 9.7e-7, 7.6e-8, 0, 0, 0;
      -7.1e-8, 7.6e-8, 8.1e-8, 0, 0, 0;
      0, 0, 0, 2.8e-11, 0, 0;
      0, 0, 0, 0, 2.7e-11, 0;
      0, 0, 0, 0, 0, 9.6e-12];

% Earth gravity constant
GM = cspice_bodvrd('Earth','GM',1); % [km^3/s^2]
%%
% EXERCISE 3.1.A

% Observation window start and ending epochs
t0_window = cspice_str2et('2010-08-12 05:30:00.000 UTC');
tf_window = cspice_str2et('2010-08-12 06:30:00.000 UTC');

% Observation window time vector
npoints = round((tf_window-t0_window)/5)+1;
etvec_window = linspace(t0_window, tf_window, npoints);

% Propagate
for i=1:npoints
    
    % Define the time boundaries
    et_END = etvec_window(i);

    % Propagation with keplerian motion
    [xf, ~, ~, ~] = KeplerianPropagator(et_sep, rv0, et_END, GM, true);
    rv_sat_sat1_window(:,i) = xf;
end

% Measuring state, retreiving azimuth and elevation measurenments
[~, az_sat_sat1, el_sat_sat1] = station_meas(station, etvec_window, rv_sat_sat1_window); % [rad] Svalbard

% Visibility windows logical array  
visible_meas = el_sat_sat1 >= el_min; 

% Indices of visible measurenment
visible_meas_id = find(visible_meas);
% Indices of not visible measurenment
nvisible_meas_id = find(~visible_meas);

% EXERCISE 3.1.B

% Reading Two-Line Element and retreiving epoch from TLE for Satellite 1
whichconst =  72; % Parameter for SGP4
[sat_sat1_spg4, ~, ~] = read_3LE(36599, 'tle\36599.3le', whichconst); 
[year,mon,day,hr,min,sec] = invjday(sat_sat1_spg4.jdsatepoch, sat_sat1_spg4.jdsatepochf);
sat_sat1_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_sat1_epoch_et = cspice_str2et(sat_sat1_epoch_str);

% Creating measurements as captured from Svalbard Station using SPG4
% propagation for Satellite 1
[sat_sat1_r_visible, sat_sat1_az_visible, sat_sat1_el_visible, sat_sat1_etvec_visible] = measure(sat_sat1_spg4, sat_sat1_epoch_et, etvec_window, station, el_min, sigma, visible_meas);
sat_sat1_meas = [sat_sat1_r_visible';sat_sat1_az_visible';sat_sat1_el_visible'];

% EXERCISE 3.1.C 

% Mean and covariance of separation epoch
etvec_sat_sat1_UKF = [sat_sat1_epoch_et, sat_sat1_etvec_visible]; % extending the time array to incorporate the separation epoch
state_sat_sat1 = rv0;
k = 1e4;
sigma1_t0 = k*P0;

% UKF absolute state of Mango
alpha = 1e-3; 
[state_sat_sat1_filter, cov_sat_sat1_filter] = UnscentedKalmanFilter(state_sat_sat1, sigma1_t0, etvec_sat_sat1_UKF, station, sigma, sat_sat1_meas, alpha, GM);

% Get real state along the visibility window
for i = 1:length(sat_sat1_etvec_visible)
    xx_sat_sat1_real(:,i) = SPG4propagation(sat_sat1_etvec_visible(i), 1);
end

% Error computation
for i = 1:length(sat_sat1_etvec_visible)

    % For the state
    error_sat_sat1(:,i) = state_sat_sat1_filter(:,i) - xx_sat_sat1_real(:,i);
    error_r(i) = norm(error_sat_sat1(1:3,i));
    error_v(i) = norm(error_sat_sat1(4:6,i));

    % For the covariance
      % For the position
      sigma_r(i) = 3*sqrt(trace(cov_sat_sat1_filter(1:3,1:3,i)));
      % For the velocity
      sigma_v(i) = 3*sqrt(trace(cov_sat_sat1_filter(4:6,4:6,i)));
end

% Plotting the results
% Position
figure()
semilogy(sat_sat1_etvec_visible/cspice_spd(), error_r, 'r', LineWidth=3.5)
hold on
semilogy(sat_sat1_etvec_visible/cspice_spd(), sigma_r, 'b',LineWidth=3.5)
xlabel('Epoch [MJD2000]', FontSize=15)
ylabel('Position: Error and 3 Sigma [km]', FontSize=15)
legend(["Position error","3 Sigma Position"],'Location','best',FontSize=15)

% Velocity
figure()
semilogy(sat_sat1_etvec_visible/cspice_spd(), error_v, 'r', LineWidth=3.5)
hold on
semilogy(sat_sat1_etvec_visible/cspice_spd(), sigma_v, 'b', LineWidth=3.5)
xlabel('Time [MJD2000]',  FontSize=15)
ylabel('Velocity: Error and 3 Sigma [km/s]', FontSize=15)
legend(["Velocity error","3 Sigma Velocity"],'Location','best',FontSize=15)
%%
% EXERCISE 3.2.A

% Propagate TLEs sat_sat1 and sat_sat2 to t0 with SGP4

% Reading Two-Line Element and retreiving epoch from TLE for Satellite 1
whichconst =  72; % Parameter for SGP4
[sat_sat1_spg4, ~, ~] = read_3LE(36599, 'tle\36599.3le', whichconst);
[year,mon,day,hr,min,sec] = invjday(sat_sat1_spg4.jdsatepoch, sat_sat1_spg4.jdsatepochf);
sat_sat1_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_sat1_epoch_et = cspice_str2et(sat_sat1_epoch_str); 
rv_sat_sat1_et0 = SPG4propagation(t0_window, 1); % SPG4 Propagation Satellite 1

% Reading Two-Line Element and retreiving epoch from TLE for Satellite 2
whichconst =  72; % Parameter for SGP4
[sat_sat2_spg4, ~, ~] = read_3LE(36827, 'tle\36827.3le', whichconst);
[year,mon,day,hr,min,sec] = invjday(sat_sat2_spg4.jdsatepoch, sat_sat2_spg4.jdsatepochf);
sat_sat2_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_sat2_epoch_et = cspice_str2et(sat_sat2_epoch_str); 
rv_sat_sat2_et0 = SPG4propagation(t0_window, 2); % SPG4 Propagation Satellite 2

% Computing relative state at t0 
rv_rel_et0 = rv_sat_sat2_et0 - rv_sat_sat1_et0; 
rv_rel_et0_ECI = rv_rel_et0; % The used reference frame is ECI

% State is rotated from ECI to LVLH at t0
r_sat_sat1_et0 = rv_sat_sat1_et0(1:3); 
v_sat_sat1_et0 = rv_sat_sat1_et0(4:6); 
% Get rotation matrix
ECI2LVLH_ROTM = rotMatrix(r_sat_sat1_et0, v_sat_sat1_et0, t0_window, GM);

% Relative state at t0 in LVLH (initial condition for integration)
rv_rel_et0_LVLH = ECI2LVLH_ROTM * rv_rel_et0_ECI;

% Mango mean motion at t0
MM_sat_sat1 = sqrt(GM/norm(r_sat_sat1_et0)^3);

% EXERCISE 3.2.B
% Simulating FFRF measurements integrating Chohessy-Whiltshire dynamics 

% Noise matrix
sigma_FFRF = diag([(1e-5)^2, (1*cspice_rpd())^2, (1*cspice_rpd())^2]); % [km,rad,rad]^2

% Relative state at t0 in LVLH is used as initial condition for integration
state_rel_LVLH(:,1) = rv_rel_et0_LVLH; 

% Propagation in time vector
for i = 2:length(etvec_window)

    % Define the time limits
    et_END = etvec_window(i);

    % Propagation with Chohessy-Whiltshire dynamics
    [xf, ~, ~, ~] = CWpropagation(t0_window, rv_rel_et0_LVLH, et_END, MM_sat_sat1);
    state_rel_LVLH(:,i) = xf; 

end

% Measurements of FFRF are simulated for Satellite 1 and Satellite 2 
[r_rel, az_rel, el_rel] = sim_meas_FFRF(etvec_window, state_rel_LVLH, sigma_FFRF);
rel_measurenment = [r_rel';az_rel';el_rel'];

% EXERCISE 3.2.C
% Estimation of relative state mean and covariance in LVLH with UKF 

% Computation of considered interval of estimation 
visible_meas_id_LAST = visible_meas_id(end); % index of the end of the considered interval of time
window_time = etvec_window(visible_meas_id_LAST+1 : visible_meas_id_LAST+1+20*60/5); % considered interval epochs
sat_sat1_meas_interval = rel_measurenment(:,(visible_meas_id_LAST+1 : visible_meas_id_LAST+1+20*60/5)); % considered interval measurements 

% The relative state at the end of the visibility window is used as initial mean for UKF estimation 
rel_t0_LVLH = state_rel_LVLH(:,visible_meas_id_LAST); 

% Initial covariance for UKF estimation equal to the time step just before to the interval

% Unscented Transform propagation on the visibility window in ECI
% Satellite 1 
n_step = 1;
t_sat1 = etvec_window(visible_meas_id_LAST) - sat_sat1_epoch_et; % from separation time to end of visibility window
x_sat1_sep_ECI = SPG4propagation(sat_sat1_epoch_et, 1); % sgp4 state at separation time
alpha = 1e-3;
[xx_sat1_end_ECI, sigma_sat1_end_UT] = UTtransform(sat_sat1_epoch_et, x_sat1_sep_ECI, P0, n_step, t_sat1, GM, alpha);
% Satellite 2 
n_step = 1;
t_sat2 = etvec_window(visible_meas_id_LAST) - sat_sat2_epoch_et; % from separation time to end of visibility window
x_sat2_sep_ECI = SPG4propagation(sat_sat2_epoch_et, 2); % sgp4 state at separation time
alpha = 1e-3;
[xx_sat2_end_ECI, sigma_sat2_end_UT] = UTtransform(sat_sat2_epoch_et, x_sat2_sep_ECI, P0, n_step, t_sat2, GM, alpha);

% Relative state covariance at end of visibility window in ECI
sigma_rel_end_UT = sigma_sat2_end_UT - sigma_sat1_end_UT;

% SGP4 final state of Satellite 1 in ECI 
x_sat1_LAST_SPG4 = SPG4propagation(etvec_window(visible_meas_id_LAST), 1);

%  SGP4 final state of Satellite 1 is rotated between ECI and LVLH 
r_sat1_LAST_SPG4 = x_sat1_LAST_SPG4(1:3); % ECI Position
v_sat1_LAST_SPG4 = x_sat1_LAST_SPG4(4:6); % ECI Velocity
ECI2LVLH_ROTM = rotMatrix(r_sat1_LAST_SPG4, v_sat1_LAST_SPG4, etvec_window(visible_meas_id_LAST), GM);

% Relative state covariance at end of visibility window in LVLH
Sigma_rel_LAST = ECI2LVLH_ROTM * sigma_rel_end_UT * ECI2LVLH_ROTM';

% The time step immediately before the interval of time vector for
% estimation with UKF is used for the initial covariance
%P_init_UKF = Sigma_rel_LAST
% Used covariance
P_init_UKF = diag([1e-6, 1e-6, 1e-6, 1e-11, 1e-11, 1e-12]);

% Interval of time vector for estimation with UKF 
vec_t_UKF = [etvec_window(visible_meas_id_LAST), window_time];

% Relative State of Satellite 1 and Satellite 2 using UKF 
alpha = 1e-3; 
[state_rel_filter_LVLH, sigma_rel_filter_LVLH] = UnscentedKalmanFilterCW(rel_t0_LVLH, P_init_UKF, vec_t_UKF, sigma_FFRF, sat_sat1_meas_interval, alpha, GM, MM_sat_sat1);

% For each time step the real relative state is defined 
x_rel_real_LVLH = state_rel_LVLH(:,(visible_meas_id_LAST+1 : visible_meas_id_LAST+1+20*60/5));

% Computing the evolution of the error
for i = 1:length(window_time)

    % For the relative state 
    error_rel(:,i) = state_rel_filter_LVLH(:,i) - x_rel_real_LVLH(:,i);
    error_r(i) = norm(error_rel(1:3,i)); % Position
    error_v(i) = norm(error_rel(4:6,i)); % Velocity
     
    % For the covariance
    sigma_r(i) = 3*sqrt(trace(sigma_rel_filter_LVLH(1:3,1:3,i))); % Position
    sigma_v(i) = 3*sqrt(trace(sigma_rel_filter_LVLH(4:6,4:6,i))); % Velocity

end

% Plotting of position
figure()
semilogy(window_time/cspice_spd(), error_r,  'r', LineWidth=3.5)
hold on
semilogy(window_time/cspice_spd(), sigma_r, 'b', LineWidth=3.5)
xlabel('Time [MJD2000]',FontSize=15)
ylabel('Position: Error and 3 Sigma  [km]', FontSize=15)
legend(["Position error","3 Sigma Position"],'Location','best',FontSize=15)

% Plotting of velocity
figure()
semilogy(window_time/cspice_spd(), error_v, 'r', LineWidth=3.5)
hold on
semilogy(window_time/cspice_spd(), sigma_v, 'b', LineWidth=3.5)
xlabel('Time [MJD2000]',FontSize=15)
ylabel('Velocity: Error and 3 Sigma [km/s]', FontSize=15)
legend(["Velocity error","3 Sigma Velocity"],'Location','best',FontSize=15)
%%
% EXERCISE 3.3

% Linear approach is exploited for absolute covariance of Mango propagation
n_step = 1+20*60/5; 
t_sat1 = 5; 
alpha = 1e-3;
[rv_sat1_twindow, sigma_sat1_twindow] = UTtransform(etvec_window(visible_meas_id_LAST), state_sat_sat1_filter(:,end), cov_sat_sat1_filter(:,:,end), n_step, t_sat1, GM, alpha);

% Rotate relative covariance from LVLH to ECI over the interval
for i=1:length(window_time)

    % State of satellite 1 ECI
    x_sat1_ECI = SPG4propagation(window_time(i), 1);
    
    % State is rotated from ECI to LVLH
    r_sat1_LAST_SPG4 = x_sat1_ECI(1:3); % Position ECI
    v_sat1_LAST_SPG4 = x_sat1_ECI(4:6); % Velocity ECI
    ECI2LVLH_ROTM = rotMatrix(r_sat1_LAST_SPG4, v_sat1_LAST_SPG4, window_time(i), GM);

    % Covariance ECI frame
    sigma_filter_rel_ECI(:,:,i) = ECI2LVLH_ROTM\sigma_rel_filter_LVLH(:,:,i)/ECI2LVLH_ROTM';

end

% Absolute covariance of Tango in ECI-2000
sigma_sat2_filter = sigma_sat1_twindow + sigma_filter_rel_ECI;

% Plot time evolution of 3 Sigma Position and Velocity
for i = 1:length(window_time)
    sigma_r(i) = 3*sqrt(trace(sigma_sat2_filter(1:3,1:3,i)));
    sigma_v(i) = 3*sqrt(trace(sigma_sat2_filter(4:6,4:6,i)));
end

% Plot Tango Absolute Covariance Position evolution in time
figure()
semilogy(window_time/cspice_spd(),sigma_r,LineWidth=3.5)
xlabel('Time [MJD2000]', 'Interpreter','latex',FontSize=15)
ylabel('3 Sigma Position [km]', FontSize=15)

% Plot Tango Absolute Covariance Velocity evolution in time
figure()
semilogy(window_time/cspice_spd(),sigma_v,LineWidth=3.5)
xlabel('Time [MJD2000]', 'Interpreter','latex',FontSize=15)
ylabel('3 Sigma Velocity [km/s]', FontSize=15)

cspice_kclear()

%% FUNCTIONS

% Keplerian Right End Side
function [dxdt] = keplerianRHS(t, x, GM, tag_j2)

% Initialize RHS
dxdt = zeros(6,1);
    
% position
rr = x(1:3);
    
% Computing distance
dist = norm(rr);
    
% Keplerian Acceleration 
a_grav = - GM*rr/dist^3;
    
% Assembling RHS
dxdt(1:3) = x(4:6);
dxdt(4:6) = a_grav;

% Accounting for J2 effect
if tag_j2 == true

    J2 = 0.0010826269; % J2 Value
    Re = cspice_bodvrd('EARTH','RADII',3); % [km] Retreiving Earth's radius with CSpice
    conversion_time = t; % [s] Ephemeris time

    % Rotation matrix from ECI to ECEF
    rotm_ECI_ECEF = cspice_pxform('J2000','ITRF93',conversion_time);
        
    % Position in ECEF
    rr_ECEF = rotm_ECI_ECEF*rr;
    dist_ECEF = norm(rr_ECEF);

    % J2 acceleration in ECEF
    ax_j2 = 3/2*GM*J2*rr_ECEF(1)/dist_ECEF^3*(Re(1)/dist_ECEF)^2*(5*(rr_ECEF(3)/dist_ECEF)^2-1);
    ay_j2 = 3/2*GM*J2*rr_ECEF(2)/dist_ECEF^3*(Re(1)/dist_ECEF)^2*(5*(rr_ECEF(3)/dist_ECEF)^2-1);
    az_j2 = 3/2*GM*J2*rr_ECEF(3)/dist_ECEF^3*(Re(1)/dist_ECEF)^2*(5*(rr_ECEF(3)/dist_ECEF)^2-3);
    a_j2 = [ax_j2;ay_j2;az_j2];

    % J2 acceleration in ECI
    rotm_ECEF_ECI = cspice_pxform('ITRF93','J2000',conversion_time);
    a_ECI = rotm_ECEF_ECI*a_j2;

    % Add to the keplerian acceleration
    dxdt(4:6) = dxdt(4:6) + a_ECI;

end

end

% Keplerian Propagator, State Only
function [xf, tf, xx, tt]  = KeplerianPropagator(t0, x0, tf, GM, tag_j2)

% Propagating RHS
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode78(@(t,x) keplerianRHS(t,x,GM,tag_j2), [t0 tf], x0, options);

% Extract state vector 
xf = xx(end,1:6)';
tf = tt(end);

end

% Unscent Transform Propagatation, State and Covariance

% Right End Side for Cohessy Wiltshire 
function [dxdt] = cwRHS(t, x, MM_sat_sat1)

% Initialize right-hand-side
dxdt = zeros(6,1);
    
% Computing the acceleration
a_x = 3*MM_sat_sat1^2*x(1) + 2*MM_sat_sat1*x(5);
a_y = -2*MM_sat_sat1*x(4);
a_z = -MM_sat_sat1^2*x(3);
aa = [a_x;a_y;a_z];
    
% Assembling RHS
dxdt(1:3) = x(4:6);
dxdt(4:6) = aa;

end

function [state_UT, sigma_UT] = UTtransform(et0, rv0, P0, N, T, GM, alpha)

% Parameters
k = 0;
n = size(P0,1);
lambda = alpha^2*(n+k)-n; 
beta = 2; 

% Initialization
x_sat0 = rv0;
P_sat0 = P0;

% Loop over time
for i=1:N

    % Integration Time Grid
    t0 = et0 + (i-1)*T;
    tf = et0 + i*T;

    sqrt_sat = sqrtm((n+lambda)*P_sat0);

    % Selecting the sigma points
    sigmaPoints = zeros(n,2*n+1);
    sigmaPoints(:,1) = x_sat0;
    for x=1:n
        sigmaPoints(:,x+1) = x_sat0 + sqrt_sat(:,x);
    end
    for x=n+1:2*n
        sigmaPoints(:,x+1) = x_sat0 - sqrt_sat(:,x-n);
    end

    % Compute weights
    W_m = zeros(2*n+1,1);
    W_m(1) = lambda/(n+lambda);
    W_m(2:end) = 1/(2*(n+lambda));
    W_c = zeros(2*n+1,1);
    W_c(1) = lambda/(n+lambda) + (1-alpha^2+beta);
    W_c(2:end) = 1/(2*(n+lambda));

    % Propagate sigma points
    sigma_points = zeros(n,2*n+1);
    for x=1:2*n+1
        [yf, ~, ~, ~] = KeplerianPropagator(t0, sigmaPoints(:,x), tf, GM, true);
        sigma_points(:,x) = yf;
    end

    % Weighted Sampled Mean  
    x_mean = 0;
    for x=1:2*n+1
        x_mean = x_mean + W_m(x)*sigma_points(:,x);
    end
    state_UT(:,i) = x_mean;

    % Covariance
    cov_sat = zeros(n);
    for x=1:2*n+1
        cov_sat = cov_sat + W_c(x)*(sigma_points(:,x)-x_mean)*(sigma_points(:,x)-x_mean)';
    end
    sigma_UT(:,:,i) = cov_sat;

    % For loop updates: Initial Mean, Initial Covariance
    x_sat0 = x_mean;
    P_sat0 = cov_sat;

end

end

% Rotation Matrix from ECI to LVLV
function [Rot_Matrix] = rotMatrix(r_sat,v_sat,et,GM)

% Dividing state of the spacecraft into position and velocity
rr = r_sat;
vv = v_sat;

% Computing Norm of the position
r = norm(rr);

% Computing the angular momentum
HH = cross(rr,vv);
H = norm(HH);
dxdt = keplerianRHS(et, [rr;vv], GM, true);
aa = dxdt(4:6);
cross_rr_aa = cross(rr,aa);

% r
i_rot = rr/r; 
k_rot = HH/H; 
j_rot = cross(k_rot,i_rot); 
r_rot = [i_rot';j_rot';k_rot'];
    
% r prime
i_rot_prime = 1/r*(vv-dot(rr/r,vv)*rr/r);
k_rot_prime = 1/H*(cross_rr_aa-dot(HH/H,cross_rr_aa)*HH/H);
j_rot_prime = cross(k_rot_prime,i_rot) + cross(k_rot,i_rot_prime);
r_rot_prime = [i_rot_prime';j_rot_prime';k_rot_prime'];

% Assembling the rotation matrix
Rot_Matrix = [r_rot, zeros(3); r_rot_prime, r_rot];

end

function [xf] = SPG4propagation(et, sat) 

% Needed constants and parameters
arcsec2rad = pi/(180*3600);
whichconst =  72; 

% Get TLE data
if sat== 1
    [sat_rec, ~, ~] = read_3LE(36599, 'tle\36599.3le', whichconst);
elseif sat == 2
    [sat_rec, ~, ~] = read_3LE(36827, 'tle\36827.3le', whichconst);
end
    
% Retreiving TLE epoch 
[year,mon,day,hr,min,sec] = invjday(sat_rec.jdsatepoch, sat_rec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% Propagation 
t_since = (et-sat_epoch_et)/60;
[~, rteme, vteme] = sgp4(sat_rec, t_since); % State in TEME
% Parameters to correct nutation, to be tuned according to time frame, 
ddpsi = -0.073296*arcsec2rad; % [rad] Delta Psi 
ddeps = -0.009373*arcsec2rad; % [rad] Delta Epsilon 
ttt = cspice_unitim(et, 'ET', 'TDT')/cspice_jyear()/100; % Accounting for Precession
[reci, veci, ~] = teme2eci(rteme, vteme, [0;0;0], ttt, ddpsi, ddeps); % State in ECI
xf = [reci; veci];

end

% Functions Simulating Measurenments of FFRF 
function [r_simulated, az_simulated, el_simulated] = sim_meas_FFRF(t_vec, rel_state, sigma_FFRF)

% Computing range, aziuth and elevation 
for i = 1:length(t_vec)
    [r_measurenment(i),az_measurenment(i),el_measurenment(i)] = cspice_reclat(rel_state(1:3,i));
end

% Measurenment are simulated adding noise 
R = mvnrnd([r_measurenment',az_measurenment',el_measurenment'], sigma_FFRF);
r_simulated  = R(:,1);
az_simulated = R(:,2);
el_simulated = R(:,3);

end

% Relative State Propagator using Cohessy Wiltshire Dynamics
function [xf, tf, xx, tt]  = CWpropagation(t0,x0,tf,MM_sat_sat1)

tof = tf - t0;
    
% Perform integration
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode78(@(t,x) cwRHS(t,x,MM_sat_sat1), [0 tof], x0, options);

% Extract state vector and State Transition Matrix
xf = xx(end,1:6)';
tf = tt(end);

end




% Uscented Kalman Filter for exploiting the state of the satellites  
function [state_sat_filter, sigma_sat_est] = UnscentedKalmanFilter(xx_sat_t0, sigma_sat_t0, etvec_UKF, gs, cov, sat_sat1_meas, alpha, GM)
    
% Parameters for Unscent Transform
k = 0;
n = size(sigma_sat_t0,1);
lambda = alpha^2*(n+k)-n; 
beta = 2; 
m = size(sat_sat1_meas,1); % Measurements Vector

% Weights Computation
W_m = zeros(2*n+1,1);
W_m(1) = lambda/(n+lambda);
W_m(2:end) = 1/(2*(n+lambda));
W_c = zeros(2*n+1,1);
W_c(1) = lambda/(n+lambda) + (1-alpha^2+beta);
W_c(2:end) = 1/(2*(n+lambda));

% Previous step, a posteriori mean and covariance
rv_post_plus = xx_sat_t0;
sigma_post_plus = sigma_sat_t0;

% Loop over measurements
for i = 1:length(etvec_UKF)-1
    
    % Times of previous and next measurements
    et_prev = etvec_UKF(i);
    et_next = etvec_UKF(i+1);

    % Computing the sigma points and associated weights
    matrix_decomposed = sqrtm((n+lambda)*sigma_post_plus);
    sigma_x_post = zeros(n,2*n+1);
    sigma_x_post(:,1) = rv_post_plus;
    for j=1:n
        sigma_x_post(:,j+1) = rv_post_plus + matrix_decomposed(:,j);
    end
    for j=n+1:2*n
        sigma_x_post(:,j+1) = rv_post_plus - matrix_decomposed(:,j-n);
    end

    % Propagating the sigma points
    sigma_x_prop = zeros(n,2*n+1);
    for j=1:2*n+1
        [xf, ~, ~, ~] = KeplerianPropagator(et_prev, sigma_x_post(:,j), et_next, GM, true);
        sigma_x_prop(:,j) = xf;
    end
        
    % Computing the measurements associated to the sigma points
    meas_sigma_x = zeros(m,2*n+1);
    for j=1:2*n+1
        [range,azimuth,elevation] = station_meas(gs, et_next, sigma_x_prop(:,j));
        meas_sigma_x(:,j) = [range;azimuth;elevation];
    end
        
    % Computing the mean a priori 
    sigma_x_prop_minus = zeros(n,1);
    meas_sigma_x_minus = zeros(m,1);
    for j=1:2*n+1
        sigma_x_prop_minus = sigma_x_prop_minus + W_m(j)*sigma_x_prop(:,j);
        meas_sigma_x_minus = meas_sigma_x_minus + W_m(j)*meas_sigma_x(:,j);
    end
        
    % Compute a priori covariances
    sigma_next_minus = zeros(n);
    sigma_est_plus = zeros(m) + cov;
    meas_sigma_est_plus = zeros(n,m);
    for j=1:2*n+1
        sigma_next_minus = sigma_next_minus + W_c(j)*(sigma_x_prop(:,j)-sigma_x_prop_minus)*(sigma_x_prop(:,j)-sigma_x_prop_minus)';
        sigma_est_plus = sigma_est_plus + W_c(j)*(meas_sigma_x(:,j)-meas_sigma_x_minus)*(meas_sigma_x(:,j)-meas_sigma_x_minus)';
        meas_sigma_est_plus = meas_sigma_est_plus + W_c(j)*(sigma_x_prop(:,j)-sigma_x_prop_minus)*(meas_sigma_x(:,j)-meas_sigma_x_minus)';
    end

    % Computing the estimation a posteriori 
    K_next = meas_sigma_est_plus/sigma_est_plus;
    diff_meas = [sat_sat1_meas(1,i)-meas_sigma_x_minus(1), angdiff(sat_sat1_meas(2,i),meas_sigma_x_minus(2)), angdiff(sat_sat1_meas(3,i),meas_sigma_x_minus(3))];
    sigma_x_prop_plus = sigma_x_prop_minus + K_next*diff_meas';
    P_next_plus = sigma_next_minus - K_next*sigma_est_plus*K_next';

    % Assembling data structure
    state_sat_filter(:,i) = sigma_x_prop_plus;
    sigma_sat_est(:,:,i) = P_next_plus;

    % Updates
    rv_post_plus = sigma_x_prop_plus;
    sigma_post_plus = P_next_plus;
end

end

% Function creating the simulated measurement from ground station 
function [r_sim_vw, az_sim_vw, el_sim_vw, t_vect_window] = measure(sat_rec, sat_epoch_et, etvec_window, gs, el_station, sigma, sat_station_vw)

% Constant for arcseconds to radians conversions
arcsec2rad = pi/(180*3600);

% Array of a priori visibility ephemeris times
t_vect_window = etvec_window(sat_station_vw);

% Visibility states in TEME
for x=1:length(t_vect_window)
    t_since = (t_vect_window(x)-sat_epoch_et)/60; % [min]
    [~, r_window_TEME(:,x), v_window_TEME(:,x)] = sgp4(sat_rec,t_since);
end

% Conversion from TEME to ECI
% Parameters needed for nutation correction 
ddpsi = -0.073296*arcsec2rad; % [rad] Delta-Psi 
ddeps = -0.009373*arcsec2rad; % [rad] Delta-Epsilon 
for x=1:length(t_vect_window) 
    ttt = cspice_unitim(t_vect_window(x), 'ET', 'TDT')/cspice_jyear()/100; % correction for precession
    [r_eci_vw(:,x), v_eci_vw(:,x), ~] = teme2eci(r_window_TEME(:,x), v_window_TEME(:,x), [0;0;0], ttt, ddpsi, ddeps);
end

    % Creating the measurements from the ground station
    [r_sim, az_sim, el_sim] = station_meas(gs, t_vect_window, [r_eci_vw;v_eci_vw]); % [km,rad,rad]

    % Simulating true measurements adding gaussian noise to the created
    % measurements
    R = mvnrnd([r_sim',az_sim',el_sim'], sigma);
    r_sim_noised = R(:,1);
    az_sim_noised = R(:,2);
    el_sim_noised = R(:,3);

    % Filtering not visible measurement
    i = el_sim_noised >= el_station; 
    r_sim_vw = r_sim_noised(i);
    az_sim_vw = az_sim_noised(i);
    el_sim_vw = el_sim_noised(i);
    t_vect_window = t_vect_window(i);

end

% Measuring: Range, Azimuth, Elevation, from s/c state and
% station name
function [r,az,el] = station_meas(station_name, t_vector, state)

% Station reference frame 
topoFrame = [station_name, '_TOPO'];

% Loop over time
for i=1:length(t_vector)
    
    % Transformation from ECI to topocentric frame
    rot_m = cspice_sxform('J2000', topoFrame, t_vector(i));

    % Compute station position in ECI
    r_s_ECI = cspice_spkezr(station_name, t_vector(i), 'J2000', 'NONE', 'EARTH');

    % Compute station-satellite vector in ECI
    r_ss_ECI = state(:,i) - r_s_ECI;

    % Convert state into topocentric frame
    r_s_TOPO = rot_m*r_ss_ECI;

    % Compute range, azimuth and elevation using cspice_reclat
    [r(i),az(i),el(i)] = cspice_reclat(r_s_TOPO(1:3));

end

end

% Uscented Kalman Filter for exploiting the relative state of the satellites
function [state_rel_filter_LVLH, sigma_rel_filter_LVLH] = UnscentedKalmanFilterCW(x_rel_init, P_rel_init, vec_t_UKF, sigma_FFRF, sat_sat1_meas, alpha, GM, MM_sat_sat1)

% Parameters for the transform
k = 0;
n = size(P_rel_init,1);
lambda = alpha^2*(n+k)-n; % scaling factor
beta = 2; % for the weights
m = size(sat_sat1_meas,1); % size of the measurement vector

% Compute weights
W_m = zeros(2*n+1,1);
W_m(1) = lambda/(n+lambda);
W_m(2:end) = 1/(2*(n+lambda));
W_c = zeros(2*n+1,1);
W_c(1) = lambda/(n+lambda) + (1-alpha^2+beta);
W_c(2:end) = 1/(2*(n+lambda));

% Taking the a posteriori mean and covariance from the previous step 
rv_post_plus = x_rel_init;
sigma_post_plus = P_rel_init;

% For all the measurements
for i = 1:length(vec_t_UKF)-1

    % Time of previous and next measurements
    et_prev = vec_t_UKF(i);
    et_next = vec_t_UKF(i+1);

    % Creating the sigma points
    matrix_decomposed = sqrtm((n+lambda)*sigma_post_plus);
    sigma_x_post = zeros(n,2*n+1);
    sigma_x_post(:,1) = rv_post_plus;
    for j=1:n
        sigma_x_post(:,j+1) = rv_post_plus + matrix_decomposed(:,j);
    end
    for j=n+1:2*n
        sigma_x_post(:,j+1) = rv_post_plus - matrix_decomposed(:,j-n);
    end

    % Propagating the sigma points
    sigma_x_prop = zeros(n,2*n+1);
    for j=1:2*n+1
        [xf, ~, ~, ~] = CWpropagation(et_prev, sigma_x_post(:,j), et_next, MM_sat_sat1);
        sigma_x_prop(:,j) = xf;
    end
        
    % Computing the measurements associated to the sigma points
    meas_sigma_x = zeros(m,2*n+1);
    for j=1:2*n+1
        [r,az,el] = cspice_reclat(sigma_x_prop(1:3,j));
        meas_sigma_x(:,j) = [r;az;el];
    end
        
    % Computing the mean a priori 
    sigma_x_prop_minus = zeros(n,1);
    meas_sigma_x_minus = zeros(m,1);
    for j=1:2*n+1
        sigma_x_prop_minus = sigma_x_prop_minus + W_m(j)*sigma_x_prop(:,j);
        meas_sigma_x_minus = meas_sigma_x_minus + W_m(j)*meas_sigma_x(:,j);
    end
         
    % Compute a priori covariances
    sigma_next_minus = zeros(n);
    sigma_est_plus = zeros(m) + sigma_FFRF;
    meas_sigma_est_plus = zeros(n,m);
    for j=1:2*n+1
        sigma_next_minus = sigma_next_minus + W_c(j)*(sigma_x_prop(:,j)-sigma_x_prop_minus)*(sigma_x_prop(:,j)-sigma_x_prop_minus)';
        sigma_est_plus = sigma_est_plus + W_c(j)*(meas_sigma_x(:,j)-meas_sigma_x_minus)*(meas_sigma_x(:,j)-meas_sigma_x_minus)';
        meas_sigma_est_plus = meas_sigma_est_plus + W_c(j)*(sigma_x_prop(:,j)-sigma_x_prop_minus)*(meas_sigma_x(:,j)-meas_sigma_x_minus)';
    end

    % Computing the estimation a posteriori 
    K_next = meas_sigma_est_plus/sigma_est_plus;
    diff_meas = [sat_sat1_meas(1,i)-meas_sigma_x_minus(1), angdiff(sat_sat1_meas(2,i),meas_sigma_x_minus(2)), angdiff(sat_sat1_meas(3,i),meas_sigma_x_minus(3))];
    sigma_x_prop_plus = sigma_x_prop_minus + K_next*diff_meas';
    P_next_plus = sigma_next_minus - K_next*sigma_est_plus*K_next';

    % Assembling data structure
    state_rel_filter_LVLH(:,i) = sigma_x_prop_plus;
    sigma_rel_filter_LVLH(:,:,i) = P_next_plus;

    % Updates
    rv_post_plus = sigma_x_prop_plus;
    sigma_post_plus = P_next_plus;

end

end

function [state_UT, sigma_UT] = UTtransform(et0, rv0, P0, N, T, GM, alpha)

% Parameters
k = 0;
n = size(P0,1);
lambda = alpha^2*(n+k)-n; 
beta = 2; 

% Initialization
x_sat0 = rv0;
P_sat0 = P0;

% Loop over time
for i=1:N

    % Integration Time Grid
    t0 = et0 + (i-1)*T;
    tf = et0 + i*T;

    sqrt_sat = sqrtm((n+lambda)*P_sat0);

    % Selecting the sigma points
    sigmaPoints = zeros(n,2*n+1);
    sigmaPoints(:,1) = x_sat0;
    for x=1:n
        sigmaPoints(:,x+1) = x_sat0 + sqrt_sat(:,x);
    end
    for x=n+1:2*n
        sigmaPoints(:,x+1) = x_sat0 - sqrt_sat(:,x-n);
    end

    % Compute weights
    W_m = zeros(2*n+1,1);
    W_m(1) = lambda/(n+lambda);
    W_m(2:end) = 1/(2*(n+lambda));
    W_c = zeros(2*n+1,1);
    W_c(1) = lambda/(n+lambda) + (1-alpha^2+beta);
    W_c(2:end) = 1/(2*(n+lambda));

    % Propagate sigma points
    sigma_points = zeros(n,2*n+1);
    for x=1:2*n+1
        [yf, ~, ~, ~] = KeplerianPropagator(t0, sigmaPoints(:,x), tf, GM, true);
        sigma_points(:,x) = yf;
    end

    % Weighted Sampled Mean  
    x_mean = 0;
    for x=1:2*n+1
        x_mean = x_mean + W_m(x)*sigma_points(:,x);
    end
    state_UT(:,i) = x_mean;

    % Covariance
    cov_sat = zeros(n);
    for x=1:2*n+1
        cov_sat = cov_sat + W_c(x)*(sigma_points(:,x)-x_mean)*(sigma_points(:,x)-x_mean)';
    end
    sigma_UT(:,:,i) = cov_sat;

    % For loop updates: Initial Mean, Initial Covariance
    x_sat0 = x_mean;
    P_sat0 = cov_sat;

end

end