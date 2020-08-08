
% -----------------------------------------------------------------
%  main_orchard_ps_sinusoidal.m
%
%  This script is the main file for a program that performs a
%  parametric study in the nonlinear dynamics of an orchad 
%  sprayer tower, modeled as an inverted pendulum over a vehicle 
%  suspension. The system deterministic dynamics evolves according
%  to the following system of differential equations
%
%    [M]Qacce + [C]Qvelo + [N]Qvelo^2 + [K]Qdisp = g - h,
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
%  vol. 103 p. 417?426, 2009
% ----------------------------------------------------------------- 
%  programmers: Americo Barbosa da Cunha Junior
%               americo.cunhajr@gmail.com
%
%  last update: Feb 14, 2016
% -----------------------------------------------------------------

clear all
close all
clc


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Orchad Sprayer Tower Dynamics                      ')
disp(' (parametric study)                                 ')
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
case_name = 'orchad_param_study_A_03';
%case_name = 'orchad_param_study_A_04';
%case_name = 'orchad_param_study_A_05';

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


A   = 0.5;       % forcing amplitude (m)
w   = 1.0;       % forcing frequency (rad/s)
rho = pi/9;      % forcing phase shift (rad)
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


y1_0      = 0.0; % initial suspension displacement
phi1_0    = 0.0; % initial suspension rotation
phi2_0    = 0.0; % initial tower rotation
y1dot_0   = 0.0; % initial suspension velocity
phi1dot_0 = 0.0; % initial suspension angular velocity
phi2dot_0 = 0.0; % initial tower angular velocity


% physical parameters vector
%phys_param = [A w rho kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];
% -----------------------------------------------------------



% -----------------------------------------------------------
k_t = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 ...
       150 200 250 300 350 400 450 500]*1.0e3;

w   = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];


% -----------------------------------------------------------





% integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- running the parametric study  --- ');
disp(' ');
disp('    ... ');
disp(' ');


t0 = 0.0;   % initial time (s)
t1 = 1.5e3; % final time (s)

% initial conditions
IC = [y1_0 phi1_0 phi2_0 y1dot_0 phi1dot_0 phi2dot_0]';

% ODE solver optional parameters
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

% preallocate memory for the system responde tables
y1   = zeros(length(k_t),length(w));
phi1 = zeros(length(k_t),length(w));
phi2 = zeros(length(k_t),length(w));


for iw = 1:length(w)

    for ikt = 1:length(k_t)
        
        % physical parameters vector
        phys_param = [A w(iw) rho k_t(ikt) ct k1 k2 c1 c2 ...
                      B1 B2 g m1 m2 L1 L2 I1 I2];

        % ODE solver Runge-Kutta45
        [time,Y]=ode45(@(t,x)orchard(t,x,phys_param),[t0 t1],IC);
        
        % number of time steps
        Ndt = length(time);

        % number of steps for steady state
        Nss = round(0.8*Ndt);
        
        % save the maximum amplitude values
          y1(ikt,iw) = max(abs(Y(Nss:Ndt,1)));
        phi1(ikt,iw) = max(abs(Y(Nss:Ndt,3)));
        phi2(ikt,iw) = max(abs(Y(Nss:Ndt,5)));
        
    end
end

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


% plot maximum amplitude of y1 vs kt
% ...........................................................
gtitle = ' suspension displacement';
leg1   = 'w = 1 rad/s';
leg2   = 'w = 2 rad/s';
leg3   = 'w = 3 rad/s';
leg4   = 'w = 4 rad/s';
xlab   = ' torsional stiffness (N m/rad)';
ylab   = ' maximum amplitude (m)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
ymax   = 1.0;
gname  = 'y1_vs_kt_A_03';
flag   = 'eps';
fig1a  = graph_type4(k_t,y1(:,1),...
                     k_t,y1(:,2),...
                     k_t,y1(:,3),...
                     k_t,y1(:,4),...
                     gtitle,leg1,leg2,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1a);
% ...........................................................


% plot maximum amplitude phi1 vs kt
% ...........................................................
gtitle = ' suspension rotation';
leg1   = 'w = 1 rad/s';
leg2   = 'w = 2 rad/s';
leg3   = 'w = 3 rad/s';
leg4   = 'w = 4 rad/s';
xlab   = ' torsional stiffness (N m/rad)';
ylab   = ' maximum amplitude (rad)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
ymax   = 0.4;
gname  = 'phi1_vs_kt_A_03';
flag   = 'eps';
fig2a  = graph_type4(k_t,phi1(:,1),...
                     k_t,phi1(:,2),...
                     k_t,phi1(:,3),...
                     k_t,phi1(:,4),...
                     gtitle,leg1,leg2,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2a);
% ...........................................................



% plot maximum amplitude of phi2 vs kt
% ...........................................................
gtitle = ' tower rotation';
leg1   = 'w = 1 rad/s';
leg2   = 'w = 2 rad/s';
leg3   = 'w = 3 rad/s';
leg4   = 'w = 4 rad/s';
xlab   = ' torsional stiffness (N m/rad)';
ylab   = ' maximum amplitude (rad)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
ymax   = 2.0;
gname  = 'phi2_vs_kt_A_03';
flag   = 'eps';
fig3a  = graph_type4(k_t,phi2(:,1),...
                     k_t,phi2(:,2),...
                     k_t,phi2(:,3),...
                     k_t,phi2(:,4),...
                     gtitle,leg1,leg2,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3a);
% ...........................................................



% plot resonance curves
% ...........................................................
gtitle = ' resonance curves';
leg1   = 'susp. displacement ';
leg2   = 'susp. rotation ';
leg3   = 'tower rotation ';
xlab   = ' excitation frequency (rad/s)';
ylab   = ' maximum amplitude';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
ymax   = 8.0;
gname  = 'ampl_vs_freq_A_03';
flag   = 'eps';
fig2  = graph_type3(w,  y1(8,:),...
                    w,phi1(8,:),...
                    w,phi2(8,:),...
                    gtitle,leg1,leg2,leg3,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2);
% ...........................................................



% resonance frequencies
% ...........................................................
[  y1_pks,  y1_freq] = findpeaks(  y1(8,:));
[phi1_pks,phi1_freq] = findpeaks(phi1(8,:));
[phi2_pks,phi2_freq] = findpeaks(phi2(8,:));

disp(' suspension displacement resonance frequencies')
disp(' ')
disp(w(y1_freq));
disp(' ')

disp(' suspension rotation resonance frequencies')
disp(' ')
disp(w(phi1_freq));
disp(' ')

disp(' tower rotation resonance frequencies')
disp(' ')
disp(w(phi2_freq));
disp(' ')
% ...........................................................

toc
% -----------------------------------------------------------

