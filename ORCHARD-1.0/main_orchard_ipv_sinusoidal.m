
% -----------------------------------------------------------------
%  main_orchard_ivp_sinusoidal.m
%
%  This is the main file for a program that simulates the nonlinear
%  dynamics of an orchad sprayer tower which modeled as an inverted 
%  pendulum over a vehicle suspension. The system dynamics evolves 
%  according to the following system of differential equations
%
%    [M]Qacce + [C]Qvelo + [N]Qvelo.^2 + [K]Qdisp = g - h,
%  
%  where
%  
%   [M]   - configuration dependent mass matrix
%   [C]   - configuration dependent damping matrix
%   [N]   - configuration dependent circulatory matrix
%   [K]   - configuration dependent stiffness matrix
%   g     - configuration dependent gravity vector
%   h     - configuration dependent harmonic external forcing
%   Qdisp - generalized displacement vector
%   Qvelo - generalized velocity vector
%   Qacee - generalized acceleration vector
%  
%  Reference:
%  
%  S. Sartori Junior, J. M. Balthazar and B. R. Pontes Junior
%  Non-linear dynamics of a tower orchard sprayer based on an
%  inverted pendulum model
%  Biosystems Engineering
%  vol. 103 p. 417-426, 2009
% ----------------------------------------------------------------- 
%  programmers: Americo Barbosa da Cunha Junior
%               americo.cunhajr@gmail.com
%
%  last update: Feb 15, 2017
% -----------------------------------------------------------------

clear all
close all
clc


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Orchad Sprayer Tower Dynamics                      ')
disp(' (nonlinear dynamics integration)                   ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
disp(' Jorge Luis Palacios Felix                          ')
disp(' jorge.felix@unipampa.edu.br                        ')
disp('                                                    ')
disp(' Jose Manoel Balthazar                              ')
disp(' jmbaltha@gmail.com                                 ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
%case_name = 'orchard_ivp_sinusoidal_A_100mm';
case_name = 'orchard_ivp_sinusoidal_A_150mm';
%case_name = 'orchard_ivp_sinusoidal_A_200mm';
%case_name = 'orchard_ivp_case1';
%case_name = 'orchard_ivp_case2';
%case_name = 'orchard_ivp_case3';
%case_name = 'orchard_ivp_case4';
%case_name = 'orchard_ivp_case5';
%case_name = 'orchard_ivp_case6';
%case_name = 'orchard_ivp_case7';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

kt  = 45.0e3;    % junction torcional stiffness (N/rad)
ct  = 50.0e3;    % junction torcional damping (N/rad/s)
k1  = 465.0e3;   % left  wheel stiffness (N/m)
k2  = 465.0e3;   % right wheel stiffness (N/m)
c1  = 5.6e3;     % left  wheel damping (N/m/s)
c2  = 5.6e3;     % right wheel damping (N/m/s)
B1  = 850.0e-3;  % left  wheel distance from center line stiffness (m)
B2  = 850.0e-3;  % right wheel distance from center line stiffness (m)
g   = 9.81;      % gravity acceleration (m/s^2)
m1  = 6500;      % suspension mass (kg)
m2  = 800;       % tower mass (kg)
L1  = 0.2;       % distance between P and suspension GC (m)
L2  = 2.4;       % tower length (m)
I1  = 6850;      % suspension moment of inertia (kg m^2)
I2  = 6250;      % tower moment of inertia (kg m^2)
A   = 150e-3;    % forcing amplitude (m)
w   = 2*pi*1.0;  % forcing frequency (rad/s)
rho = pi/4;      % forcing phase shift (rad)


% velocity of translation (km/h)
vtrans_km_h = 12.0;

% tire diameter (m)
Dtire = 1.0;

% displacement natural frequency (rad/s)
wn1 = sqrt((k1+k2)/(m1+m2));

% static equilibrium configuration
yst = -(m1+m2)*g/(k1+k2);

% initial suspension displacement (m)
y1_0 = yst;

% initial suspension rotation (rad)
phi1_0 = 0.0;

% initial tower rotation (rad)
phi2_0 = 0.0;

% initial suspension velocity (m/s)
y1dot_0 = 0.0;

% initial suspension angular velocity (rad/s)
phi1dot_0 = 0.0;

% initial tower angular velocity (rad/s)
phi2dot_0 = 0.0;
% -----------------------------------------------------------



% integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nonlinear dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');


% physical parameters vector
phys_param = [kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2 A w rho];

% initial conditions
IC = [y1_0; phi1_0; phi2_0; y1dot_0; phi1dot_0; phi2dot_0];

% initial time (s)
t0 = 0.0;

% final time (s)
t1 = t0 + 30.0;

% interval of analysis (s)
tspan = [t0 t1];

% optional parameter for ODE solver
% (RelTol = 1.0e-3 AbsTol = 1.0e-6)
%opt = odeset ('RelTol',1.0e-3,'AbsTol',1.0e-6);

% ODE solver Runge-Kutta45
[time,y]=ode45(@(t,y)orchard_rhs(t,y,phys_param),tspan,IC);

% time series of suspension displacement (m)
Qy1 = y(:,1);

% time series of suspension rotation (rad)
Qphi1 = y(:,2);

% time series of tower rotation (rad)
Qphi2 = y(:,3);

% time series of suspension velocity (m/s)
Qy1dot = y(:,4);

% time series of suspension angular velocity (rad/s)
Qphi1dot = y(:,5);

% time series of tower angular velocity (rad/s)
Qphi2dot = y(:,6);

% time series of tower horizontal displacement
Qx2 = - L1*sin(Qphi1) - L2*sin(Qphi2);

% time series of tower vertical displacement
Qy2 = Qy1 + L1*cos(Qphi1) + L2*cos(Qphi2);

% time series of tower horizontal velocity
Qx2dot = - L1*cos(Qphi1).*Qphi1dot - L2*cos(Qphi2).*Qphi2dot;

% time series of tower vertical velocity
Qy2dot = Qy1dot - L1*sin(Qphi1).*Qphi1dot - L2*sin(Qphi2).*Qphi2dot;

% number of time steps
Ndt = length(time);

% number of steps for steady state
Nss = round(0.8*Ndt);

toc
% -----------------------------------------------------------



% save simulation results
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------



% post processing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');

disp('trailer vertical displacement maximum amplitude')
max(abs(Qy1(Nss:Ndt)))
disp(' ')

disp('trailer rotation maximum amplitude')
max(abs(Qphi1(Nss:Ndt)))
disp(' ')

disp('tower rotation maximum amplitude')
max(abs(Qphi2(Nss:Ndt)))
disp(' ')

disp('tower horizontal displacement maximum amplitude')
max(abs(Qx2(Nss:Ndt)))
disp(' ')
 
disp('tower vertical displacement maximum amplitude')
max(abs(Qy2(Nss:Ndt)))
disp(' ')
 
%post_orchard_ivp;


% animate orchard sparayer dynamics
% ...........................................................
vtitle = 'Orchard Sprayer Tower Dynamics';
legend = 'harmonic excitation';
xmin   = -2;
xmax   =  2;
ymin   = -2;
ymax   =  4;
vname  = [num2str(case_name),'__animation'];
fig00  = post_orchard_animation(time,Qy1,Qphi1,Qphi2,...
                                B1,B2,L1,L2,Dtire,vtrans_km_h,...
                                vtitle,[],legend,xmin,xmax,ymin,ymax);
%close(fig00);
% ...........................................................

toc
% -----------------------------------------------------------



