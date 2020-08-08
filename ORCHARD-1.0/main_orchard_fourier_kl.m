
% -----------------------------------------------------------------
%  main_orchard_fourier_kl.m
%
%  This is the main file for a program to do spectral analysis
%  of the nonlinear dynamics of an orchad sprayer tower which 
%  is modeled as an inverted pendulum over a vehicle suspension. 
%  The system dynamics evolves according to the following system 
%  of differential equations
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
%  last update: Feb 16, 2017
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



% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

kt  = 100.0e3;   % junction torcional stiffness (N/rad)
ct  = 40.0e3;    % junction torcional damping (N/rad/s)
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

% natural frequencies
omega1 = sqrt((k1+k2)/(m1+m2));
omega2 = sqrt((kt)/(I1+m2*L1^2));
omega3 = sqrt((kt)/(I2+m2*L2^2));

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


% signal processing parameters
% -----------------------------------------------------------
disp(' ');
disp(' --- define signal processing parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');


% minimun frequency of the band (Hz)
freq_min = 0.0;

% maximum frequency of the band  (Hz)
freq_max = 5.0;

% Nyquist sampling frequency  (Hz)
freq_samp = 2*freq_max;

% time step or sampling period (s)
dt = 1/freq_samp;

% number of samples in a signal copy
%Nsamp = 256;
%Nsamp = 512;
Nsamp = 1024;
%Nsamp = 2048;
%Nsamp = 4096;

% period of a signal copy (s)
Tsignal = Nsamp*dt;

% upper bound of variance estimator
%eps_v = 0.01;

% number of signal copies
%Ncopies = round(2/eps_v)+1;
Ncopies = 200;

% initial time of analysis (s)
t0 = 0.0;

% final time of analysis (s)
t1 = t0 + (Ncopies+1)*Tsignal;

% number of time steps
Ndt = (Ncopies+1)*Nsamp;
% -----------------------------------------------------------



% generate external forcing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- generating external forcing --- ');
disp(' ');
disp('    ... ');
disp(' ');

% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30082016);
RandStream.setDefaultStream(rng_stream); % Matlab 2009
%RandStream.setGlobalStream(rng_stream); % Matlab 2013

% number of samples
Ns = 1;

% velocity of translation (m/s)
vtrans = vtrans_km_h*(1000/3600);

% translation domain upper limit (m)
xupp = 50.0;

% correlation length (m)
%corrlen = 0.1;
corrlen = 1.0;
%corrlen = 10.0;
%corrlen = 50.0;
%corrlen = 100.0;

% spacial mean of random process Ye1 and Ye2 (m)
mu_x_Ye1 = 0.5;
mu_x_Ye2 = 0.5;

% coeficient of variation
delta = 0.35;

% translation std. deviation (m)
sig_x = delta*mu_x_Ye1;

% number of eigenpairs to be computed
Neig = 500;

% computing the autocorr function eigenpairs
[lambda,phi,phi_dot,omega,...
    time_KL,iroot,ising,iter] = ...
            fredholm_expcorr_eig(xupp/vtrans,corrlen/vtrans,Neig,Ndt);

% generate the uncorrelated random variables
ksi1 = randn(Ns,Neig);
ksi2 = randn(Ns,Neig);

% multiply eigenvalues by variance
lambda = lambda*(sig_x/vtrans)^2;

% representation energy level
energy_level = 0.999;

% number of eigenpairs used in KL expansion
[~,Nkl] = max(cumsum(lambda)/sum(lambda) >= energy_level);

% KL representation of random process Ye1 (Ndt x Ns)
Ye1 = mu_x_Ye1/vtrans + ...
        phi(:,1:Nkl)*sqrt(diag(lambda(1:Nkl)))*ksi1(:,1:Nkl)';

% KL representation of random process Ye2 (Ndt x Ns)
Ye2 = mu_x_Ye2/vtrans + ...
        phi(:,1:Nkl)*sqrt(diag(lambda(1:Nkl)))*ksi2(:,1:Nkl)';

% KL representation of random process Ye1_dot (Ndt x Ns)
Ye1dot = phi_dot(:,1:Nkl)*sqrt(diag(lambda(1:Nkl)))*ksi1(:,1:Nkl)';

% KL representation of random process Ye2_dot (Ndt x Ns)
Ye2dot = phi_dot(:,1:Nkl)*sqrt(diag(lambda(1:Nkl)))*ksi2(:,1:Nkl)';

toc
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
phys_param = [kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];

% initial conditions
IC = [y1_0; phi1_0; phi2_0; y1dot_0; phi1dot_0; phi2dot_0];

% time instants of analysis (s)
tspan = linspace(t0,t1,Ndt);

% optional parameter for ODE solver
% (RelTol = 1.0e-3 AbsTol = 1.0e-6)
opt = odeset ('RelTol',1.0e-3,'AbsTol',1.0e-6);

% ODE solver Runge-Kutta45
[time,y] = ode45(@(t,y)orchard_rhs_kl(t,y,phys_param,...
                                      Ye1(:,1),Ye2(:,1),...
                                      Ye1dot(:,1),Ye2dot(:,1),tspan),...
                                      tspan,IC);

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

% time series of tower horizontal displacement/velocity
Qx2    = - L1*sin(Qphi1) ...
         - L2*sin(Qphi2);
Qx2dot = - L1*cos(Qphi1).*Qphi1dot ...
         - L2*cos(Qphi2).*Qphi2dot;

% time series of tower vertical displacement/velocity
Qy2    = Qy1 ...
         + L1*cos(Qphi1) ...
         + L2*cos(Qphi2);
Qy2dot = Qy1dot ...
         - L1*sin(Qphi1).*Qphi1dot ...
         - L2*sin(Qphi2).*Qphi2dot;

% number of steps for steady state
Nss = round(0.5*Ndt);

time_elapsed = toc
% -----------------------------------------------------------


% spectral analysis of the signals
% -----------------------------------------------------------
Nfft = 2^(nextpow2(Nsamp));

%[psd_Qx2,freq] = pwelch(  Qx2(Nsamp+1:end),Nsamp,Nsamp/2,Nfft,freq_samp);
%[psd_Qy1     ] = pwelch(  Qy1(Nsamp+1:end),[],[],Nfft,freq_samp);
%[psd_Qy2     ] = pwelch(  Qy2(Nsamp+1:end),[],[],Nfft,freq_samp);
%[psd_Qphi1   ] = pwelch(Qphi1(Nsamp+1:end),[],[],Nfft,freq_samp);
%[psd_Qphi2   ] = pwelch(Qphi2(Nsamp+1:end),[],[],Nfft,freq_samp);

% power spectral density
[psd_Qx2  ,freq] = signal_psd(  Qx2(Nsamp+1:end),freq_max,freq_samp,Nfft,Nsamp,Ncopies);
[psd_Qy1  ,freq] = signal_psd(  Qy1(Nsamp+1:end),freq_max,freq_samp,Nfft,Nsamp,Ncopies);
[psd_Qy2  ,freq] = signal_psd(  Qy2(Nsamp+1:end),freq_max,freq_samp,Nfft,Nsamp,Ncopies);
[psd_Qphi1,freq] = signal_psd(Qphi1(Nsamp+1:end),freq_max,freq_samp,Nfft,Nsamp,Ncopies);
[psd_Qphi2,freq] = signal_psd(Qphi2(Nsamp+1:end),freq_max,freq_samp,Nfft,Nsamp,Ncopies);

% % Savitzky?Golay filtered power spectral density
% psd_Qy1_sg   = sgolayfilt(psd_Qy1,  5,41);
% psd_Qx2_sg   = sgolayfilt(psd_Qx2,  5,41);
% psd_Qy2_sg   = sgolayfilt(psd_Qy2,  5,41);
% psd_Qphi1_sg = sgolayfilt(psd_Qphi1,5,41);
% psd_Qphi2_sg = sgolayfilt(psd_Qphi2,5,41);
% -----------------------------------------------------------


% simulation information
% -----------------------------------------------------------
%case_name = 'orchard_fourier_kl_corr_01m_v_12kmph';
%case_name = 'orchard_fourier_kl_corr_1m_v_12kmph';
%case_name = 'orchard_fourier_kl_corr_10m_v_12kmph';
%case_name = 'orchard_fourier_kl_corr_50m_v_12kmph';
%case_name = 'orchard_fourier_kl_corr_100m_v_12kmph';

case_name  = ['orchard_fourier_kl',...
              '_corr_',num2str(corrlen),'m',...
              '_v_',num2str(vtrans_km_h),'kmph'];

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
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

 
post_orchard_fourier_kl;


toc
% -----------------------------------------------------------



