% SGN - ASSIGNMENT #1
% Emanuele Marzorati 10724126 225268
% Exercise 3: Continuous Guidance

%clc; clear; close all;

% Load kernels
cspice_furnsh('ex02.tm');

% 3 Solving problem with first set of parameters

% Integration frame string 
center = 'SSB';
frame = 'ECLIPJ2000'; % Ecliptic plane
GM = cspice_bodvrd('Sun','GM',1); % Solar gravitational constant [km^3/s^2]

% Spacecraft parameters
m0 = 1000; % [kg]
Tmax = 0.8; % [N]
Isp = 3120; % [s]

% Reference parameters for adimensionalization
MU = m0; % [kg]
LU = cspice_convrt(1,'AU','km'); % [km]
tref = sqrt(LU^3/GM); % [s]

% Adimensional parameters
Adim.MU = MU;
Adim.LU = LU;
Adim.tref = tref;

Adim.m0 = m0/MU;
Adim.Tmax = Tmax*tref^2/(MU*LU*1e3);
Adim.Isp = Isp/tref;
Adim.mu = 1;
Adim.g0 = 9.80665*1e-3*tref^2/LU;

% Initial time
init_epoch_str = '2023-05-28-14:13:09.000 UTC';
t0 = cspice_str2et(init_epoch_str);
t0_adim = t0/tref;

% Initial state vector
x0 = zeros(7,1);
x0(1:6) = cspice_spkezr('Earth',t0, frame, 'NONE', center);
x0(7) = m0;

% Dimensionless initial state vector
x0_adim = zeros(7,1);
x0_adim(1:3) = x0(1:3)/LU;
x0_adim(4:6) = x0(4:6)*tref/LU;
x0_adim(7) = x0(7)/MU;

% While loop
exitflag = -1; % set a negative exitflag to enter the loop
varsol = zeros(8,1);
iter = 0;
while exitflag<=0

    iter = iter + 1

    % Guess on initial costate vector
    L0 = -20 + (20+20)*rand(7,1);

    % Guess on delta time
    Dt_adim = 2*pi*rand(1,1);
    
    % Unknowns: var=[lr0(3x1), lv0(3x1), lm0(1), tf(1x1)]
    var0 = [L0;Dt_adim];
    
    % Determine fun to solve
    fun = @(var)paramfun(var,x0_adim,Adim,frame,center,t0_adim);
    
    % Call fsolve with fun and var0
    [varsol, fval, exitflag] = fsolve(fun,var0);

end

% Solution for Lambda0
sol_L0 = varsol(1:7);

% Solution for tf

sol_Dt_adim = varsol(8);
sol_Dt_secs = sol_Dt_adim*tref;
sol_Dt_days = sol_Dt_secs/(3600*24);

sol_tf_adim = t0_adim + sol_Dt_adim;
sol_tf = sol_tf_adim*tref;
sol_tf_date = cspice_et2utc(sol_tf,'C',4)

% Propagate solution to check for H(t)=constant
[~,~,yy,~] = propagate_y(0,[x0_adim;sol_L0],sol_Dt_adim,Adim);
for i=1:size(yy,1)
    r = yy(i,1:3);
    v = yy(i,4:6);
    m = yy(i,7);
    Lr = yy(i,8:10);
    Lv = yy(i,11:13);
    Lm = yy(i,14);
    
    HH(i) = 1 + dot(Lr,v) - Adim.mu*dot(r,Lv)/norm(r)^3 + Adim.Tmax/(Adim.Isp*Adim.g0)*(-norm(Lv)*Adim.Isp*Adim.g0/m - Lm);
end

% Computing Delta Position and Velocity between Venus and Spacecraft
e_r = norm(fval(1:3))*LU
e_v = norm(fval(4:6))*LU/tref     

HH
sol_L0

% 4 

% Numerical Continuation
T0_var = linspace(0.75,0.5,6);

% Initial Guesses

for i = 1:length(T0_var)

    % Adimensional parameters
    Adim_step.MU = MU;
    Adim_step.LU = LU;
    Adim_step.tref = tref;

    Adim_step.m0 = m0/MU;
    Adim_step.Tmax = T0_var(i)*tref^2/(MU*LU*1e3);
    Adim_step.Isp = Isp/tref;
    Adim_step.mu = 1;
    Adim_step.g0 = 9.80665*1e-3*tref^2/LU;

    % Initial guess, equal to solution of the optimization for T = 0.8
    guess_costate = sol_L0'; % guess on costate
    guess_dtime = sol_Dt_adim; % guess on time

    % While loop
    exitflag = -1; % set a negative exitflag to enter the loop
    varsol = zeros(8,1);
    iter = 0;

    while exitflag<=0

        iter = iter + 1

        % Guess on initial costate vector, equal to inital costate vector
        % of solution for T = 0.8
        L0_step = guess_costate';

        % Guess on delta time, equal to delta time of solution for T = 0.8
        Dt_adim_step = guess_dtime;

        % Unknowns: var=[lr0(3x1), lv0(3x1), lm0(1), tf(1x1)]
        var0_step = [L0_step;Dt_adim_step];

        % Determine fun to solve
        fun = @(var)paramfun(var,x0_adim,Adim_step,frame,center,t0_adim);

        % Call fsolve with fun and var0
        [varsol, fval, exitflag] = fsolve(fun,var0_step);
    end

    % Update for initial costate
    guess_costate = varsol(1:7);
    
    % Update for initial dtime
    guess_dtime = varsol(8);

    % Solution for tf
    sol_Dt_adim = varsol(8);
    sol_Dt_secs = sol_Dt_adim*tref;
    sol_Dt_days_step = sol_Dt_secs/(3600*24) 
    
    sol_tf_adim = t0_adim + sol_Dt_adim;
    sol_tf = sol_tf_adim*tref;
    sol_tf_date_step = cspice_et2utc(sol_tf,'C',4)

end

sol_L0_step = varsol(1:7)
e_r_step = norm(fval(1:3))*LU
e_v_step = norm(fval(4:6))*LU/tref

%cspice_kclear

% Functions

% Zero finding problem Objective Function

function F = paramfun(var,x0_adim,Adim,frame,center,t0_adim)
    
    % Adimensional parameters
    MU = Adim.MU;
    LU = Adim.LU;
    tref = Adim.tref;
    m0 = Adim.m0;
    Tmax = Adim.Tmax;
    Isp = Adim.Isp;
    mu = Adim.mu;
    g0 = Adim.g0;

    % Variables
    L0 = var(1:7);
    Dt_adim = var(8);

    % Propagation
    [yf,~,~,~] = propagate_y(0,[x0_adim;L0],Dt_adim,Adim);
    rf = yf(1:3);
    vf = yf(4:6);
    mf = yf(7);
    Lrf = yf(8:10);
    Lvf = yf(11:13);
    Lmf = yf(14);
    
    % F(1:6): SC and Venus states coincide at tf
    % a - Calling cspice (YES)
    rvf_venus = cspice_spkezr('Venus',(t0_adim+Dt_adim)*tref,frame,'NONE',center);
    % b - Propagating
%     rv0_venus = cspice_spkezr('Venus',t0_adim*tref,frame,'NONE',center);
%     [rvf_venus, ~, ~, ~]  = propagate_tbp(0,rv0_venus,Dt_adim,frame,center);
    
    rf_venus_adim = rvf_venus(1:3)/LU;
    vf_venus_adim = rvf_venus(4:6)*tref/LU;

    F(1:3) = rf - rf_venus_adim;
    F(4:6) = vf - vf_venus_adim;

    % F(7) Lambda_m is zero at tf
    F(7) = Lmf;

    % F(8) Transversality condition
    Hf_adim = 1 + dot(Lrf,vf) - 1/norm(rf)^3*dot(rf,Lvf) + Tmax/(Isp*g0)*(-(norm(Lvf)/mf)*(Isp*g0)-Lmf);

    %%% a) Using cspice (YES)
    accf_venus_adim = -rf_venus_adim/norm(rf_venus_adim)^3;
    %%% b) Calling tbp dynamics
%     dxdt_venus = tbp_rhs(Dt_adim,rvf_venus,frame,center);
%     accf_venus_adim = dxdt_venus(4:6)*tref^2/LU;

    F(8) = Hf_adim - dot(Lrf,vf_venus_adim) - dot(Lvf,accf_venus_adim);
    
end


% Y=[X;Lambda] Propagator
function [yf,tf,yy,tt]  = propagate_y(t0,y0,tf,Adim)

    % Time span
    tof = tf - t0;

    % Perform integration
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, yy] = ode78(@(t,y) dynamics_y(t, y, Adim), [0 tof], y0, options_STM);

    % Extract state vector and State Transition Matrix
    yf = yy(end,:)';
    tf = tt(end);

end

% Y=[X;Lambda] RHS 
function [dydt] = dynamics_y(t, y, Adim)

    % Adimensional parameters
    MU = Adim.MU;
    LU = Adim.LU;
    tref = Adim.tref;
    m0 = Adim.m0;
    Tmax = Adim.Tmax;
    Isp = Adim.Isp;
    mu = Adim.mu;
    g0 = Adim.g0;

    % State
    xr=y(1:3); xv=y(4:6); xm=y(7);
    
    % Costate 
    lr=y(8:10); lv=y(11:13); lm=y(14);
    
    % Initialize right-hand-side 
    dydt = zeros(14,1);
    
    % Compute square distance and distance
    dist2 = dot(xr,xr);
    dist = sqrt(dist2);
    
    % Dynamics of the state
    dydt(1:3) = xv;
    dydt(4:6) = - mu*xr/dist^3 - (Tmax/xm)*(lv/norm(lv));
    dydt(7)   = - Tmax/(Isp*g0);
    
    % Dynamics of the costate
    dydt(8:10)  = - 3*mu*dot(xr,lv)*xr/dist^5 + mu*lv/dist^3;
    dydt(11:13) = - lr;
    dydt(14)    = - (norm(lv)*Tmax/xm^2);

end

% Venus RHS for Two-Body Problem
function [dxdt] = tbp_rhs(t,x,frame,center)

    % Initialize right-hand-side
    dxdt = zeros(6,1);

    % Extract the asteroid position
    rr_venus = x(1:3);

    % Retrieve position of the body
    rr_sun = cspice_spkpos('Sun',t,frame,'NONE',center);

    % Extract asteroid position wrt. body
    rr_venus_sun = rr_venus - rr_sun;

    % Compute square distance
    dist = norm(rr_venus_sun);

    % Compute the gravitational acceleration using Newton's law
    GM = cspice_bodvrd('Sun','GM',1);

    % Assemble right-hand side (state)
    dxdt(1:3) = x(4:6);
    dxdt(4:6) = - GM*(rr_venus_sun)/(dist^3);

end

% Venus Propagator for Two-Body Problem
function [xf, tf, xx, tt]  = propagate_tbp(t0,x0,tf,frame,center)

    tof = tf - t0;
    
    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);
    [tt, xx] = ode78(@(t,x) tbp_rhs(t,x,frame,center), [0 tof], x0, options);

    % Extract final state vector and time
    xf = xx(end,1:6)';
    tf = tt(end);

end