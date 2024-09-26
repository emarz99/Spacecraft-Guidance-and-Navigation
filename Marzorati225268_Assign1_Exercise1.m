clc, clearvars, close all
format long g

% SGN Assignment 1
% Emanuele Marzorati, 225268, 10724126
% Exercise 1

%% EXERCISE 1.1

% Circular Restricted Three-Body Problem Dynamics
% Using the symbolic toolkit for matlab in order to compute the lagrangian mechanics
% of the circular restricted three-Body problem dynamics
syms x y z mu

r1 = sqrt((x + mu)^2     + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);
dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
dUdz = -(1-mu)/r1^3*z - mu/r2^3*z;

% positioning on x axis
% using the command subs to evaluate in the specified y z mu values 
m=0.012150;
fun = subs(dUdx,[y z mu],[0 0 m]);
% giving initial guesses for the zero finding function
guess = [-1 0 1];

for i=1:3
    l{i}.coordinates.x = vpasolve(fun,x,guess(i));
    l{i}.coordinates.y = 0;
    l{i}.name = ["Lagrange Point ",num2str(i)];
end

% printing these results to make sure we have identified the correct
% lagrangian point
fplot(fun,'m','LineWidth',1.5)
grid on
hold on
plot3(1-m,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
plot3(-m,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
for i=1:3
     
     plot(l{i}.coordinates.x, l{i}.coordinates.y,".",'MarkerSize',25)
 end
%title('Collinear Lagrangian Points'),xlabel('X'),ylabel('Y'),ylim([-30 30])
hold off
l{2}.coordinates.x

%% 
% EX 1.2
% ITERATIVELY FINDING THE CORRECTION TO BE APPLIED

clc; clear; close all;

mu = 0.012150;

x0 = 1.08892819445324;
y0 = 0;
z0 = 0.0591799623455459;
vx0 = 0;
vy0 = 0.257888699435051;
vz0 = 0;

xx0 = [x0; y0; z0; vx0; vy0; vz0]; % first guess

% Propagation until time tf
tf = 2; 
%[~,tf,xx_F,~]  = propagate(0,xx0,tf,mu,true); % To update tf with the time of intersection with y=0 in the first propagation (unmodified initial state)
% Forward and backward propagation
[~,~,xx_F,~]  = propagate(0,xx0,tf,mu,true); 
[~,~,xx_B,~]  = propagate(0,xx0,-tf,mu,true);
figure
hold on
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
scatter3(xx_F(1,1),xx_F(1,2),xx_F(1,3),'k*')
% these to show position of the attractors
%plot3(1-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
%plot3(-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
xlabel('x [-]',FontSize=10)
ylabel('y [-]',FontSize=10)
zlabel('z [-]',FontSize=10)
legend(["Forward","Backward"],'Location','best',FontSize=10)

% Initialization of the differential correction while loop
Nmax    = 50;   % Maximum number of iterations to avoid it to get stuck
iter    = 0;    
tol     = 1e-8; % Error tolerance for the final velocities Vxf, Vzf

% We will impose the stopping conditions based upon how much we are near to
% the final state after the corrections, so the errors, we are setting these 
% to a high value not to stop at first iteration, higher than the tolerance
err_vxf = 1;   
err_vzf = 1;   

% Set the update to the first guess for our three control parameters
x0_new = x0;
vy0_new = vy0;

% This is the first term of the right side of the equation, this is the
% final state after the correction, this is what I'm aiming at
vxf_ref = 0;   % this could be dropped since the reference value we aim is zero 
vzf_ref = 0;

% Analytical while loop 
while (abs(err_vxf)>tol || abs(err_vzf)>tol) && iter < Nmax % setting the condition to exit from the while loop
    
    % Perform propagation up to event time
    [xf,PHI,te,~,~]  = propagate_STM(0,[x0_new;y0;z0;vx0;vy0_new;vz0],tf,mu,true);

    % Compute the deviation in the final state (x,z, velocity)
    err_vxf = xf(4) - vxf_ref;
    err_vzf = xf(6) - vzf_ref;

    % Compute the correction
    delta_vy = (-xf(6)+(PHI(6,1)*xf(4))/PHI(4,1))/(PHI(6,5)-(PHI(6,1)*PHI(4,5))/PHI(4,1));
    delta_x = -(PHI(4,5)*delta_vy+xf(4))/PHI(4,1);
    vy0_new = vy0_new + delta_vy;
    x0_new = x0_new + delta_x;

    % Update iteration counter
    iter = iter+1;

end

% Initializing the new starting conditions
xx0_new = [x0_new;y0;z0;vx0;vy0_new;vz0];

% Plot Forward-Backward propagation
[~,te,xx_F,tt_F] = propagate(0,xx0_new, te,mu,false); % Obtain 'te' as the final event time
[~,~,xx_B,tt_B] = propagate(0,xx0_new, -te,mu,false);
figure
hold on
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
scatter3(xx_F(1,1),xx_F(1,2),xx_F(1,3),'k*')
% these to show position of the attractors
% plot3(1-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
% plot3(-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
xlabel('x [-]',FontSize=10)
ylabel('y [-]',FontSize=10)
zlabel('z [-]',FontSize=10)
legend(["Forward","Backward"],'Location','best',FontSize=10)

% Plot the full forward propagation
[xf,te,xx,tt] = propagate(0,xx0_new,2*te,mu,false);
figure
hold on
plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',2)
scatter3(xx(1,1),xx(1,2),xx(1,3),'k*')
plot3(1-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
plot3(-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
xlabel('x [-]',FontSize=10)
ylabel('y [-]',FontSize=10)
zlabel('z [-]',FontSize=10)

% Compute the error of the final state w.r.t. the initial state (shall coincide to acomplish the periodicity)
err_x0 = (xf-xx0_new);
max_mean_error = abs(err_x0)./max(abs(xx0_new),abs(xf)); % Just to avoid divisions by zero if reference value is zero
fprintf('\nPosition error [km]: %+.3e %+.3e %+.3e\nVelocity error [km]: %+.3e %+.3e %+.3e\nRelative error [km]: %.3e %.3e %.3e %.3e %+.3e %+.3e\n', err_x0(1:3),err_x0(4:6),max_mean_error);
xx0_new

%% 1.3 Families of periodic halo orbits
clc; clear; close all;

mu = 0.012150;

% Initial guess
x0 = 1.0902780529;
y0 = 0;
z0_var = linspace(0.05917996234,0.034,6);
vx0 = 0;
vy0 = 0.260349391891;
vz0 = 0;

% Initialization
Nmax    = 50;   % Maximum number of iterations to avoid it to get stuck
tol     = 1e-8; % Error tolerance for the final velocities Vxf, Vzf

% Propagation time tf
tf = 2; % First guess shall exceed the time of intersection y=0

% % First guess of the initial variables
x0_new = x0;
vy0_new = vy0;

% Constraints of the final state (to create a periodic orbit)
vxf_ref = 0;
vzf_ref = 0;

% For loop to create family
figure
hold on
for i = 1:length(z0_var)

    % Select x0
    z0 = z0_var(i);

    % Restore initialization
    err_vxf = 1;    % High value not to stop at first iteration
    err_vzf = 1;    % High value not to stop at first iteration
    iter    = 0;    % Iteration counter

    % While loop
    while (abs(err_vxf)>tol || abs(err_vzf)>tol) && iter < Nmax
        
        % Perform propagation up to event time
        [xf,PHI,te,~,~]  = propagate_STM(0,[x0_new; y0; z0; vx0; vy0_new; vz0],tf,mu,true);
    
        % Compute the deviation in the final state (x,z, velocity)
        err_vxf = xf(4) - vxf_ref;
        err_vzf = xf(6) - vzf_ref;
    
        % Compute the correction
        Dvy0 = (-xf(6)+(PHI(6,1)*xf(4))/PHI(4,1))/(PHI(6,5)-(PHI(6,1)*PHI(4,5))/PHI(4,1));
        Dx0 = -(PHI(4,5)*Dvy0+xf(4))/PHI(4,1);
        vy0_new = vy0_new + Dvy0;
        x0_new = x0_new + Dx0;
        
        % Update iteration counter
        iter = iter+1;
    
    end
    
    xx0_new = [x0_new;y0;z0;vx0;vy0_new;vz0];
    
    % Plot the full forward propagation
    [xf,PHI,te,xx,tt] = propagate_STM(0,xx0_new,2*te,mu,false);
    
    plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',2)
    scatter3(xx(1,1),xx(1,2),xx(1,3),'k*')
    xlabel('x [-]',FontSize=10)
    ylabel('y [-]',FontSize=10)
    zlabel('z [-]',FontSize=10)

end
plot3(1-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
plot3(-mu,0,0,'o','Marker', 'o', 'MarkerFaceColor','k','MarkerSize',8)
hold off

%% FUNCTIONS


% Propagate State
function [xf, tf, xx, tt]  = propagate(t0,x0,tf,mu,varargin)

    if nargin>4
        evtFlag=varargin{1};
    else
        evtFlag=true;
    end

    tof = tf - t0;
    
    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(tt,xx) y_plane_crossing(tt,xx,evtFlag));
    [tt, xx] = ode78(@(t,x) xyzCR3BP(t,x,mu), [0 tof], x0, options);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    tf = tt(end);

end

function [dxdt] = xyzCR3BP(~,xx,mu)
    
    % Extract variables
    x = xx(1);
    y = xx(2);
    z = xx(3);
    vx = xx(4);
    vy = xx(5);
    vz = xx(6);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);

    % Assemble right-hand side
    dxdt = zeros(6,1);

    dUdx = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
    dUdy = y - y*(1-mu)/r1^3 - y*mu/r2^3;
    dUdz = -z*(1-mu)/r1^3 - z*mu/r2^3;

    dxdt(1:3) = xx(4:6);
    dxdt(4) = 2*vy + dUdx;
    dxdt(5) = -2*vx + dUdy;
    dxdt(6) = dUdz;
    
end

% Propagate State and STM
function [xf, PHIf, tf, xx, tt]  = propagate_STM(t0,x0,tf,mu,varargin)

    if nargin>4
        evtFlag=varargin{1};
    else
        evtFlag=true;
    end

    tof = tf - t0;

    % Initialize State Transition Matrix at t0
    Phi0 = eye(6);

    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0; Phi0(:)];
    
    % Perform integration
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(tt,xx) y_plane_crossing(tt,xx,evtFlag));
    [tt, xx] = ode78(@(t,x) xyzCR3BP_STM(t,x,mu), [0 tof], x0Phi0, options_STM);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    PHIf = reshape(xx(end,7:end),6,6);
    tf = tt(end);

end

% RHS and STM for 3D-CRTBP 
function [dxdt] = xyzCR3BP_STM(~,xx, mu)
    
    % Extract variables
    x = xx(1);
    y = xx(2);
    z = xx(3);
    vx = xx(4);
    vy = xx(5);
    vz = xx(6);

    % Put PHI in matrix form
    Phi = reshape(xx(7:end),6,6);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);

    % Assemble the matrix A(t)=dfdx 6x6 matrix
    G = zeros(3);
    G(1,1) = 1 - (1-mu)/r1^3 + 3*(1-mu)*(x+mu)^2/r1^5 - mu/r2^3 + 3*mu*(x+mu-1)^2/r2^5;
    G(2,2) = 1 - (1-mu)/r1^3 + 3*(1-mu)*y^2/r1^5 - mu/r2^3 + 3*mu*y^2/r2^5;
    G(3,3) = -(1-mu)/r1^3 + 3*(1-mu)*z^2/r1^5 - mu/r2^3 + 3*mu*z^2/r2^5;
    G(1,2) = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5; G(2,1) = G(1,2);
    G(1,3) = 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x+mu-1)*z/r2^5; G(3,1) = G(1,3);
    G(2,3) = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5; G(3,2) = G(2,3);

    H = zeros(3);
    H(1,2) = 2; H(2,1) = -H(1,2);

    dfdx = [zeros(3), eye(3); G, H];
    
    % Compute the derivative of the STM
    Phidot = dfdx*Phi;

    % Assemble right-hand side
    dxdt = zeros(42,1);

    dUdx = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
    dUdy = y - y*(1-mu)/r1^3 - y*mu/r2^3;
    dUdz = -z*(1-mu)/r1^3 - z*mu/r2^3;

    dxdt(1:3) = xx(4:6);
    dxdt(4) = 2*vy + dUdx;
    dxdt(5) = -2*vx + dUdy;
    dxdt(6) = dUdz;
    dxdt(7:end) = Phidot(:);
    
end

% Termination Event
function [value, isterminal, direction] = y_plane_crossing(~,xx,isTerminal)

    value = xx(2);
    isterminal = isTerminal;
    direction = 0;

end