clc, clear all, close all
format long g

% SGN Assignment 2
% Emanuele Marzorati, 225268, 10724126
% Exercise 1

%% INITIALIZE

% Load kernels
cspice_furnsh('assignment02.tm');
%%
% Initial Data

% Define input UTC (coordinated universal time = civil time at Greenwich):
ref_epoch_str = '2010-08-12 05:27:39.1140 TDB';

% Time conversion 
et0 = cspice_str2et(ref_epoch_str);

% Initial States
% Sat 1
r0_sat1 = [4622.232026629, 5399.3369588058, -0.0212138165769957];
v0_sat1 = [0.812221125483763,  -0.721512914578826, 7.42665302729053];
% Sat 2
r0_sat2 = [4621.69343340281, 5399.26386352847, -3.09039248714313];
v0_sat2 = [0.813960847513811,  -0.719449862738607, 7.42706066911294];

% Initial Covariance
P0 = [5.6e-7, 3.5e-7, -7.1e-8, 0, 0, 0;
      3.5e-7, 9.7e-7, 7.6e-8, 0, 0, 0;
     -7.1e-8, 7.6e-8, 8.1e-8, 0, 0, 0;
         0, 0, 0, 2.8e-11, 0, 0;
         0, 0, 0, 0, +2.7e-11, 0;
         0, 0, 0, 0, 0, +9.6e-12];

N_rev = 10; % Number of Revolutions  % Time of orbit is only the one of satellite one

GM = cspice_bodvrd('Earth', 'GM', 1); % Earth Gravitational Constant 

% Compute orbital period of Sat 1
T = Period(r0_sat1, v0_sat1, GM);

tof = N_rev*T; % Time of Flight Sat 1
etf = et0 + tof; % Final Time Sat 1

x0_1 = [r0_sat1, v0_sat1]'; % Initial State of Sat 1
x0_2 = [r0_sat2, v0_sat2]'; % Initial State of Sat 2
%%
% EXERCISE 1.1

% Time grid
tt = et0:T:(et0 + T*N_rev);

% For Satellite 1
[xx_sat1_lin, P_sat1_lin] = LinPropagation(x0_1,P0,et0,T,GM,N_rev); % Linear Evolution
[xx_sat1_UT, P_sat1_UT] = UnscentedPropagation(x0_1,P0,et0,T,GM,N_rev); % Uscented Transpose Evolution

% For the Position
for i=1:(N_rev+1)
    trace_lin_sat1_pos(i) = sqrt(trace(P_sat1_lin(1:3,1:3,i)));
    trace_unscented_sat1_pos(i) = sqrt(trace(P_sat1_UT(1:3,1:3,i)));
end
MSE_sat1_pos = 1/length(trace_lin_sat1_pos)*sum(abs(trace_lin_sat1_pos-trace_unscented_sat1_pos))

figure()
hold on 
plot(tt/cspice_spd(), trace_lin_sat1_pos,'Color','r', LineWidth=2)
plot(tt/cspice_spd(), trace_unscented_sat1_pos,'--','Color','b', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('$\sqrt{trace(P_r)}$ [km]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT"], 'Interpreter','latex','Location','best',FontSize=16)
title({'Mango Position Covariance '})
grid on

% For the Velocity
for i=1:(N_rev+1)
    trace_lin_sat1_vel(i) = sqrt(trace(P_sat1_lin(4:6,4:6,i)));
    trace_unscented_sat1_vel(i) = sqrt(trace(P_sat1_UT(4:6,4:6,i)));
end
MSE_sat1_vel = 1/length(trace_lin_sat1_vel)*sum(abs(trace_lin_sat1_vel-trace_unscented_sat1_vel))

figure()
hold on 
plot(tt/cspice_spd(), trace_lin_sat1_vel,'Color','r', LineWidth=2)
plot(tt/cspice_spd(), trace_unscented_sat1_vel,'--','Color','b', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('$\sqrt{trace(P_v)}$ [km/s]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT"],'Location','best',FontSize=16)
title({'Mango Velocity Covariance'})
grid on

% For Satellite 2
[xx_sat2_lin, P_sat2_lin] = LinPropagation(x0_2,P0,et0,T,GM,N_rev); % Linear Evolution
[xx_sat2_UT, P_sat2_UT] = UnscentedPropagation(x0_2,P0,et0,T,GM,N_rev); % Uscented Transpose Evolution

% For the Position
for i=1:(N_rev+1)
    trace_lin_sat2_pos(i) = sqrt(trace(P_sat2_lin(1:3,1:3,i)));
    trace_unscented_sat2_pos(i) = sqrt(trace(P_sat2_UT(1:3,1:3,i)));
end
MSE_sat2_pos = 1/length(trace_lin_sat2_pos)*sum(abs(trace_lin_sat2_pos-trace_unscented_sat2_pos))
figure()
hold on 
plot(tt/cspice_spd(), trace_lin_sat2_pos,'Color','r', LineWidth=2)
plot(tt/cspice_spd(), trace_unscented_sat2_pos,'--','Color','b', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('$\sqrt{trace(P_r)}$ [km]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT"],'Location','best',FontSize=16)
title({'Tango Position Covariance '})
grid on

% For the Velocity
for i=1:(N_rev+1)
    trace_lin_sat2_vel(i) = sqrt(trace(P_sat2_lin(4:6,4:6,i)));
    trace_unscented_sat2_vel(i) = sqrt(trace(P_sat2_UT(4:6,4:6,i)));
end
MSE_sat2_vel = 1/length(trace_lin_sat2_vel)*sum(abs(trace_lin_sat2_vel-trace_unscented_sat2_vel))
figure()
hold on 
plot(tt/cspice_spd(), trace_lin_sat2_vel,'Color','r', LineWidth=2)
plot(tt/cspice_spd(), trace_unscented_sat2_vel,'--','Color','b', LineWidth=2)
xlabel('Time [MJD2000]')
ylabel('$\sqrt{trace(P_r)}$ [km/s]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT"], 'Location','best',FontSize=16)
title({'Tango Velocity Covariance'})
grid on

%%
%% EXERCISE 1.2

% Linear Case
exitflag = 0;
i=2;
while  (i<11) && (exitflag==0)

    % Sum of covariance submatrices associated to position for linear propagation
    Psum_lin = P_sat1_lin(1:3,1:3,i) + P_sat2_lin(1:3,1:3,i); 

    % Norm of relative position for linear propagation
    distance_lin = norm(xx_sat1_lin(1:3,i) - xx_sat2_lin(1:3,i)); 

    delta_r_lin = 3*sqrt(max(eig(Psum_lin)));

    if distance_lin < delta_r_lin 
        exitflag = 1;
        fprintf('\nLinear Case, Number of revolutions before critical event: %+.3e \nMinimum Distance: %+.3e \nActual Distance: %+.3e',i-1,delta_r_lin,distance_lin);
    end
    i = i + 1;
end

% Unscented Trasform Case
exitflag = 0;
i=2;
while  (i<11) && (exitflag==0)

    % Sum of covariance submatrices associated to position for UT propagation
    Psum_UT = P_sat1_UT(1:3,1:3,i) + P_sat2_UT(1:3,1:3,i); 

    % Norm of relative position for UT propagation
    distance_UT = norm(xx_sat1_UT(1:3,i) - xx_sat2_UT(1:3,i));

    delta_r_UT = 3*sqrt(max(eig(Psum_UT)));

    if distance_UT < delta_r_UT
        exitflag = 1;
        fprintf('\nUT Case, Number of revolutions before critical event: %+.3e \nMinimum Distance: %+.3e \nActual Distance: %+.3e',i-1,delta_r_UT,distance_UT);
    end
    i = i + 1;
end

%%
% 3 - Monte Carlo Simulation
Samples = 100; % First Trial
% Samples = 1000; % Final Trial

% Satellite 1
[state_sat1_MC, P_sat1_MC, sigma1] = MCprop(et0, x0_1, P0, N_rev, T, GM, Samples);

% Satellite 2
[state_sat2_MC, P_sat2_MC, sigma2] = MCprop(et0, x0_2, P0, N_rev, T, GM, Samples);
%%
% 3.A - Plot and compare the expressions of the submatrices

% Satellite 1 Position
for i=1:N_rev+1
    Pr_sat1_lin(i) = 3*sqrt(max(eig(P_sat1_lin(1:3,1:3,i))));
    Pr_sat1_UT(i) = 3*sqrt(max(eig(P_sat1_UT(1:3,1:3,i))));
    Pr_sat1_MC(i) = 3*sqrt(max(eig(P_sat1_MC(1:3,1:3,i))));
end
% Mean absolute error
mae_r_sat1 = 1/length(Pr_sat1_lin)*sum(abs(Pr_sat1_lin-Pr_sat1_MC));

figure()
hold on
plot(tt/cspice_spd(),Pr_sat1_lin,'Color', 'r', LineWidth=1.5)
plot(tt/cspice_spd(),Pr_sat1_UT,'--','Color', 'b', LineWidth=1.5)
plot(tt/cspice_spd(),Pr_sat1_MC,'Color', 'k', LineWidth=1.5)
xlabel('Time [MJD2000]',FontSize=15)
ylabel('$3\sigma_r$ [km]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT","MC"])
title({'Mango Position Covariance'})
grid on


% Satellite 1 Velocity
for i=1:N_rev+1
    Pv_sat1_lin(i) = 3*sqrt(max(eig(P_sat1_lin(4:6,4:6,i))));
    Pv_sat1_UT(i) = 3*sqrt(max(eig(P_sat1_UT(4:6,4:6,i))));
    Pv_sat1_MC(i) = 3*sqrt(max(eig(P_sat1_MC(4:6,4:6,i))));
end
% Mean absolute error
mae_v_sat1 = 1/length(Pv_sat1_lin)*sum(abs(Pv_sat1_lin-Pv_sat1_MC));

figure()
hold on
plot(tt/cspice_spd(),Pv_sat1_lin,'Color', 'r', LineWidth=1.5)
plot(tt/cspice_spd(),Pv_sat1_UT,'--','Color', 'b', LineWidth=1.5)
plot(tt/cspice_spd(),Pv_sat1_MC,'Color', 'k', LineWidth=1.5)
xlabel('Time [MJD2000]',FontSize=15)
ylabel('$3\sigma_v$ [km/s]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT","MC"])
title({'Mango Velocity Covariance'})
grid on


% Satellite 2 Position
for i=1:N_rev+1
    Pr_sat2_lin(i) = 3*sqrt(max(eig(P_sat2_lin(1:3,1:3,i))));
    Pr_sat2_UT(i) = 3*sqrt(max(eig(P_sat2_UT(1:3,1:3,i))));
    Pr_sat2_MC(i) = 3*sqrt(max(eig(P_sat2_MC(1:3,1:3,i))));
end
% Mean absolute error
mae_r_sat2 = 1/length(Pr_sat2_lin)*sum(abs(Pr_sat2_lin-Pr_sat2_MC));

figure()
hold on
plot(tt/cspice_spd(),Pr_sat2_lin,'Color', 'r', LineWidth=1.5)
plot(tt/cspice_spd(),Pr_sat2_UT,'--','Color', 'b', LineWidth=1.5)
plot(tt/cspice_spd(),Pr_sat2_MC,'Color', 'k', LineWidth=1.5)
xlabel('Time [MJD2000]',FontSize=15)
ylabel('$3\sigma_r$ [km]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT","MC"])
title({'Tango Position Covariance'})
grid on

% Satellite 2 Velocity
for i=1:N_rev+1
    Pv_sat2_lin(i) = 3*sqrt(max(eig(P_sat2_lin(4:6,4:6,i))));
    Pv_sat2_UT(i) = 3*sqrt(max(eig(P_sat2_UT(4:6,4:6,i))));
    Pv_sat2_MC(i) = 3*sqrt(max(eig(P_sat2_MC(4:6,4:6,i))));
end
% Mean absolute error
mae_v_sat2 = 1/length(Pv_sat2_lin)*sum(abs(Pv_sat2_lin-Pv_sat2_MC));

figure()
hold on
plot(tt/cspice_spd(),Pv_sat2_lin,'Color', 'r', LineWidth=1.5)
plot(tt/cspice_spd(),Pv_sat2_UT,'--','Color', 'b', LineWidth=1.5)
plot(tt/cspice_spd(),Pv_sat2_MC,'Color', 'k', LineWidth=1.5)
xlabel('Time [MJD2000]',FontSize=15)
ylabel('$3\sigma_v$ [km/s]', 'Interpreter','latex',FontSize=16)
legend(["LinCov","UT","MC"])
title({'Tango Velocity Covariance'})
grid on

% 3.B Evolution of mean and covariance projected in the orbital plane
%%
% Satellite 1
satellite_1 = 'Mango';
cov_plot(xx_sat1_lin, P_sat1_lin, state_sat1_MC, P_sat1_MC, sigma1, xx_sat1_UT, P_sat1_UT, satellite_1)

% Satellite 2
satellite_2 = 'Tango';
cov_plot(xx_sat2_lin, P_sat2_lin, state_sat2_MC, P_sat2_MC, sigma2, xx_sat2_UT, P_sat2_UT, satellite_2)

%cspice_kclear()
%% Functions

% PLOT RESULTS PROJECTED IN ORBITAL PLANE 
%%

%%
function cov_plot(state_lin, P_lin, state_MC, P_MC, points_sat, state_UT, P_UT, satellite)

    % Loop over the selected revolutions
    n_rev = [2, 5, 10];
    for n=n_rev

        % Rotation Matrix from ECI to LVLH 
        eci2lvlh_ROTM = rotation_eci2lvlh(state_lin(1:3,n),state_lin(4:6,n));
        % Defing the distance from the origins
        delta = state_lin(1:3,n);

        % LINEAR COVARIANCE
        % Rotating linear mean and covariance to LVLH
        rr_lvlh_lin = eci2lvlh_ROTM * (state_lin(1:3,n) - delta); 
        P_rr_lvlh_lin = eci2lvlh_ROTM * P_lin(1:3,1:3,n) * eci2lvlh_ROTM';
        % Projecting mean and covariance on the orbital plane  
        r_sat_lin = rr_lvlh_lin(1:2);
        P_r_sat_lin = P_rr_lvlh_lin(1:2,1:2);
        % Get the LC ellipse
        [centre_x_lin,centre_y_lin] = ellipse(r_sat_lin, P_r_sat_lin);
        % UNSCENT TRANSFORM
        % Rotating Uscent Transform mean and covariance to LVLH
        rr_lvlh_sat_Ut = eci2lvlh_ROTM * (state_UT(1:3,n) - delta);
        P_rr_lvlh_sat_Ut = eci2lvlh_ROTM * P_UT(1:3,1:3,n) * eci2lvlh_ROTM';
        % Projecting Uscent Transform mean and covariance on the orbital plane
        r_sat_UT = rr_lvlh_sat_Ut(1:2);
        P_r_sat_UT = P_rr_lvlh_sat_Ut(1:2,1:2);
        % Get the Uscent Transform ellipse
        [centre_x_UT,centre_y_UT] = ellipse(r_sat_UT, P_r_sat_UT);
        % MONTE CARLO
        % Monte Carlo mean and covariance rotated to LVLH
        r_lvlh_sat_MC = eci2lvlh_ROTM * (state_MC(1:3,n) - delta);
        P_r_lvlh_sat_MC = eci2lvlh_ROTM * P_MC(1:3,1:3,n) * eci2lvlh_ROTM';
        % Projecting MC mean and covariance on the orbital plane
        r_sat_MC = r_lvlh_sat_MC(1:2);
        P_r_sat_MC = P_r_lvlh_sat_MC(1:2,1:2);
        % Get the MC ellipse
        [centre_x_MC,centre_Y_MC] = ellipse(r_sat_MC, P_r_sat_MC);
    
        % MC points from ECI to LVLH
        for i=1:size(points_sat,1) % over the different points
            points_sat_lvlh(i,:) = (eci2lvlh_ROTM * (points_sat(i,:,n)' - delta))';
        end
        % Project the MC points on the orbital plane
        points_sat_plane = points_sat_lvlh(:,1:2);
        
        figure()
        hold on
        plot(centre_x_lin,centre_y_lin, 'k', LineWidth=1.5)
        % Monte Carlo Points
        plot(points_sat_plane(:,1),points_sat_plane(:,2), 'b.', LineWidth=2)
        % Monte Carlo Mean
        plot(r_sat_MC(1),r_sat_MC(2), 'b*', LineWidth=1.5)
        % Monte Carlo Ellipse
        plot(centre_x_MC,centre_Y_MC, 'b', LineWidth=1.5)
        
        % Uscented Transpose Mean
        plot(r_sat_UT(1),r_sat_UT(2), 'r*', LineWidth=1.5)
        % Uscented Transpose Ellipse
        plot(centre_x_UT,centre_y_UT, 'r--', LineWidth=1.5)
        hold off

        axis([-8e-3 8e-3 -1.7 1.7])
        set(gca,'FontSize',17)
        xlabel('x-LVLH [km]')
        ylabel('y-LVLH [km]')
        legend(["Lin Cov","MC Samples","MC Mean","MC Cov","UT Mean","UT Cov"],'Location','best',FontSize=12)
        title_str = {"Uncertainty Propagation of:" + satellite + " Rev: " + num2str(n)};
        title(title_str)
        grid on
        
    end    
end

function [x_center,y_center] = ellipse(center, P)
    
    % Compute eig_val and eig_vect of the covariance matrix
    [eig_vect, eig_val] = eig(P);
    
    % From Covariance Matrix Extracting Semi Major and Semi Minor Axis
    SMA = 3*sqrt(eig_val(1, 1));
    sma = 3*sqrt(eig_val(2, 2));
    
    % Generate parametric param_angles
    param_angles = linspace(0, 2*pi, 1000);
    
    % Parametric equations for ellipse
    x_center = center(1) + SMA * cos(param_angles) * eig_vect(1, 1) + sma * sin(param_angles) * eig_vect(1, 2);
    y_center = center(2) + SMA * cos(param_angles) * eig_vect(2, 1) + sma * sin(param_angles) * eig_vect(2, 2);

end

function [state_MC, cov_MC, samples] = MCprop(et0, x0, P0, Nrev, T1, mu, Ns)

    % Initialization
    initial_state = x0;
    initial_cov = P0;

    % Samples generation
    S = mvnrnd(initial_state, initial_cov, Ns); % each row is a sample (Ns,6) Matrix
    samples(:,:,1) = S(:,1:3); % all rows/samples, taking only position

    % For each revolution (N)
    for i=1:Nrev

        % Define the integration time boundaries
        et_f = et0 + i*T1;

        % Propagation of the samples on the orbit
        store = zeros(size(S)); 
        for x=1:Ns
            [store(x,:), ~, ~, ~] = keplerian_propagator(et0, S(x,:)', et_f, mu);
        end
        samples(:,:,i+1) = store(:,1:3); % each sample is propagated and stored in a row

        % Computing sample mean and covariance of the mapping
        m = 1/Ns*sum(store,1); % row vector
        sample_covariance = zeros(6);
        for x=1:Ns
            sample_covariance = sample_covariance + 1/(Ns-1)*(store(x,:)'-m')*(store(x,:)'-m')';
        end

        % Store results of each time step
        state_MC(:,i) = m;
        cov_MC(:,:,i) = sample_covariance;

    end

    % Insert initial time values for mean and covariance
    state_MC = [x0, state_MC];
    cov_MC(:,:,2:end+1) = cov_MC;
    cov_MC(:,:,1) = P0;

end

function [state_UT,P_sat_UT] = UnscentedPropagation(x0, P0, et0, T, GM, N)

% Unscent Transform Tuning Parameters
alpha = 1e-1;
n = size(P0,1); % Dimension of the state vector
beta = 2;
k = 0;
lambda = (alpha^2)*(n+k) - n;

% Updates
current_state = x0; % Satellite state
current_cov = P0; % Satellite covariance

for i=1:N

    % Time integration boundaries
    et1 = et0 + (i-1)*T;
    et2 = et0 + i*T;

    % Compute 
    sqrt_sat = sqrtm((n+lambda)*current_cov);
    
    % Sigma samples composition
    sigma_sat = zeros(n,2*n+1);
    sigma_sat(:,1) = current_state;
    for x=1:n
        sigma_sat(:,x+1) = current_state + sqrt_sat(:,x);
    end
    for x=n+1:2*n
        sigma_sat(:,x+1) = current_state - sqrt_sat(:,x-n);
    end
    
    % Weights Computation
    Wm = zeros(2*n+1,1);
    Wm(1) = lambda/(n+lambda);
    Wm(2:end) = 1/(2*(n+lambda));
    Wc = zeros(2*n+1,1);
    Wc(1) = lambda/(n+lambda) + (1-alpha^2+beta);
    Wc(2:end) = 1/(2*(n+lambda));
    
    % Sigma samples Propagation
    Y_sat = zeros(n,2*n+1);
    for x=1:2*n+1
        [yf, ~, ~, ~] = keplerian_propagator(et1, sigma_sat(:,x), et2, GM);
        Y_sat(:,x) = yf;
    end
    
    % Compute the weighted sampled mean and covariance
    y_sat = 0;
    for x=1:2*n+1
        y_sat = y_sat + Wm(x)*Y_sat(:,x);
    end
    state_UT(:,i) = y_sat;
   
    Py_sat = zeros(n);
    for x=1:2*n+1
        Py_sat = Py_sat + Wc(x)*(Y_sat(:,x)-y_sat)*(Y_sat(:,x)-y_sat)';
    end
    P_sat_UT(:,:,i) = Py_sat;
    
    % Update the initial mean and covariance for the next propagation step
    current_state = y_sat;
    current_cov = Py_sat;
end

% Insert initial time values for mean and covariance
state_UT = [x0, state_UT];
P_sat_UT(:,:,2:end+1) = P_sat_UT;
P_sat_UT(:,:,1) = P0;

end



% Propagating uncertainties with Linearized Approach
function [xx_sat_lin, P_sat_lin] = LinPropagation(x0,P0,et0,T,GM,N)

% Updates
current_state = x0; % Satellite state
current_cov = P0; % Satellite covariance

for i=1:N
  
    % Time integration boundaries
    et1 = et0 + (i-1)*T;
    et2 = et0 + i*T;

    % Propagation of Dynamics and STM
    [xf,STMf,~,~,~]  = keplerianSTM_propagation(et1,current_state,et2,GM); % Compunting final state and covariance
    xx_sat_lin(:,i) = xf;
    P_sat_final = STMf*current_cov*STMf';
    P_sat_lin(:,:,i) = P_sat_final;
     
    % Updating initial state and covariance
    current_state = xf;
    current_cov = P_sat_final;

end

% Creating complete data structure
xx_sat_lin = [x0,xx_sat_lin];
P_sat_lin(:,:,2:end+1) = P_sat_lin;
P_sat_lin(:,:,1) = P0;

end

% Compute orbital period assuming keplerian motion
function [T] = Period(r0_sat, v0_sat, GM)
rho = norm(r0_sat);
vel = norm(v0_sat);

sma = -0.5*GM/( 0.5*vel.^2 - GM /rho ); % semi-maxor axis
T =  2*pi*sqrt(sma.^3./GM); % orbital period (knowing that is a closed orbit)
% T = 98.8*60 % Online found orbital period for satellite 1
end

function [xf, PHIf, tf, xx, tt]  = keplerianSTM_propagation(t0,x0,tf,mu)

    tof = tf - t0;

    % Initialize State Transition Matrix at t0
    Phi0 = eye(6);

    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0; Phi0(:)];
    
    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11; ones(36,1)*1e-12]);
    [tt, xx] = ode78(@(t,x) keplerSTM_rhs(t,x,mu), [0 tof], x0Phi0, options);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    PHIf = reshape(xx(end,7:end),6,6);
    tf = tt(end);

end

function [dxdt] = keplerSTM_rhs(t, x, mu)

    % Initialize right-hand-side
    dxdt = zeros(42,1);
    
    % Extract position
    rr = x(1:3);
    Phi = reshape(x(7:end),[6,6]);
    
    % Compute square distance and distance
    dist = norm(rr);
    
    % Compute the gravitational acceleration using Newton's law
    aa_grav = - mu*rr/dist^3;
    
    % Compute the derivative dfdv
    dfdv = 3*mu/dist^5*(rr*rr') - mu/dist^3*eye(3);
           
    % Assemble the matrix A(t)=dfdx
    A = [zeros(3), eye(3);...
         dfdv, zeros(3)];
    
    % Compute the derivative of the state transition matrix
    Phidot = A*Phi;
    
    % Extract the resultant dxdt[42,1]
    dxdt(1:3) = x(4:6);
    dxdt(4:6) = aa_grav;
    dxdt(7:end) = Phidot(:);

end

function [xf, PHIf, tf, xx, tt]  = numerical_propagate_STM(et0,x0,etf,GM)
[xf_nom, tt, xx] = keplerian_propagator(et0, x0, etf, GM);

% Initialize 
e = zeros(length(x0));
PHI_N = zeros(length(x0));
tof = etf - et0;

% Computation of STM by columns
for i = 1 : length(x0)
    e(i,i) = max(sqrt(eps), sqrt(eps)*abs(x0(i))) ;  
    xf_e = keplerian_propagator(0, x0+e(:,i), tof, GM);  
    PHI_N(:,i) = (xf_e - xf_nom) / e(i,i);
end
xf = xf_nom;
tf = tt(end);
PHIf = PHI_N;
end

function [xf, tf, xx, tt] = keplerian_propagator(et0, x0, etf, GM)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM

tof = etf-et0;  % computing time of flight in case initial state is passed with ephemeris time 

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]); % tolerance settend differently for position and velocity

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 tof], x0, options);

% Extract state vector 
xf = xx(end,1:6)';
tf = tt(end);
end

function [dxdt] = keplerian_rhs(t, x, GM) %, et0)

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3); % ECI

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is obxect's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2);% + rotm2(4:6,4:6)*a_x2;
end

function [ROTM_r_ECI2LVLH] = rotation_eci2lvlh(r_eci,rv_eci)

    % Position
    rr = r_eci;
    r = norm(rr);

    % Velocity
    vv = rv_eci;

    % Angular momentum
    hh = cross(rr,vv);
    h = norm(hh);

    % R
    i_lvlh_eci = rr/r; % unit vector of i_LVLH wrt ECI
    k_lvlh_eci = hh/h; % unit vector of j_LVLH wrt ECI
    j_lvlh_eci = cross(k_lvlh_eci,i_lvlh_eci); % unit vector of k_LVLH wrt ECI
    ROTM_r_ECI2LVLH = [i_lvlh_eci';j_lvlh_eci';k_lvlh_eci'];

end