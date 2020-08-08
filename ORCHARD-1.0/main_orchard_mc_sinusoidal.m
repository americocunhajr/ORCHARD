
% -----------------------------------------------------------------
%  main_orchard_mc_sinusoidal.m
%
%  This script is the main file for a program that simulates
%  the stochastic nonlinear dynamics, using Monte Carlo method,
%  of an orchad sprayer tower, modeled as an inverted pendulum 
%  over a vehicle suspension. The system deterministic dynamics
%  evolves according to the following system of differential 
%  equations
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
%  vol. 103 p. 417-426, 2009
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 15, 2017
% -----------------------------------------------------------------

clc
clear all
close all



% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Orchad Sprayer Tower Stochastic Dynamics           ')
disp(' (uncertanty propagation calculation)               ')
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
%case_name = 'orchard_mc_sinusoidal_A';
%case_name = 'orchard_mc_sinusoidal_w';
%case_name = 'orchard_mc_sinusoidal_rho';
case_name = 'orchard_mc_sinusoidal_A_w';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define deterministic physical parameters
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
A   = 0.15;      % forcing amplitude (m)
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



% define stochastic parameters
% -----------------------------------------------------------
mu_A    = A;      % displacement amplitude mean (m)
sigma_A = 0.2*A;  % displacement amplitude std (m)
delta_A = sigma_A/mu_A;  % displacement amplitude dispersion

w1 = 2*pi*0.0;  % support of w lower bound (rad/s)
w2 = 2*pi*2.0;  % support of w upper bound (rad/s)

mu_w = 0.5*(w1+w2);  % tire displacement frequency mean (rad/s)
w    = mu_w;         % tire displacement frequency (rad/s)
%delta_w = 0.025;    % tire displacement frequency dispersion

rho1 = 0;     % support of rho lower bound (rad)
rho2 = pi/2;  % support of rho upper bound (rad)

% phase shift mean (rad)
mu_rho = 0.5*(rho1+rho2);
rho    = mu_rho;
% -----------------------------------------------------------



% nominal model
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nominal dynamics --- ');
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

% steady state index
Nss = round(0.8*Ndt);

toc
% -----------------------------------------------------------



% generate random parameters
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- generating random parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');


% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% number of samples
%Ns = 65536;
%Ns = 16384;
%Ns = 4096;
Ns = 1024;
%Ns = 256;

% number of samples in each direction
ns = sqrt(Ns);


if strcmp(case_name,'orchard_mc_sinusoidal_A')
    
    % displacement amplitude (m)
    %A = exprnd(mu_A,[Ns,1]);
    A = gamrnd(1/delta_A^2,mu_A*delta_A^2,[Ns,1]);
    
elseif strcmp(case_name,'orchard_mc_sinusoidal_w')


    % tire displacement frequency (rad/s)
    w = w1 + (w2-w1)*rand([Ns,1]);
    %w = exprnd(mu_w,[Ns,1]);
    %w = gamrnd(1/delta_w^2,mu_w*delta_w^2,[Ns,1]);
    
elseif strcmp(case_name,'orchard_mc_sinusoidal_rho')
    
    % phase shift (rad)
    %rho = rho1 + (rho2-rho1)*rand([Ns,1]);
    rho = exprnd(mu_rho,[Ns,1]);
    %rho = gamrnd(1/delta_rho^2,mu_rho*delta_rho^2,[Ns,1]);
    
elseif strcmp(case_name,'orchard_mc_sinusoidal_A_w')
    
    % displacement amplitude (m)
    %A = exprnd(mu_A,[ns,1]);
    A = gamrnd(1/delta_A^2,mu_A*delta_A^2,[ns,1]);
    
    % tire displacement frequency (rad/s)
    w = w1 + (w2-w1)*rand([ns,1]);
    %w = exprnd(mu_w,[ns,1]);
    %w = gamrnd(1/delta_w^2,mu_w*delta_w^2,[ns,1]);
    
    % replicate the random samples to access the sample space
    A = kron(A,ones(ns,1));
    w = repmat(w,ns,1);
end

toc
% -----------------------------------------------------------

    



% Monte Carlo simulation
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Monte Carlo Simulation --- ');
disp(' ');


% preallocate memory for DoFs
  y1 = zeros(Ns,Ndt);
phi1 = zeros(Ns,Ndt);
phi2 = zeros(Ns,Ndt);

% preallocate memory for DoFs time derivatives
  y1dot = zeros(Ns,Ndt);
phi1dot = zeros(Ns,Ndt);
phi2dot = zeros(Ns,Ndt);

% preallocate memory for MC conv metric
normL2_Q = zeros(1,Ns);


for imc=1:Ns
    
    if mod(imc,sqrt(Ns)) == 0
        disp('')
        disp(imc)
    end

    % physical paramters vector
	if strcmp(case_name,'orchard_mc_sinusoidal_A')
        
        phys_param = [A(imc) w rho kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];
        
    elseif strcmp(case_name,'orchard_mc_sinusoidal_w')
            
        phys_param = [A w(imc) rho kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];
        
    elseif strcmp(case_name,'orchard_mc_sinusoidal_rho')
        
        phys_param = [A w rho(imc) kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];
        
    elseif strcmp(case_name,'orchard_mc_sinusoidal_A_w')
        
        phys_param = [A(imc) w(imc) rho kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];
	end
    
    
    % ODE solver Runge-Kutta45
    [time_MC,Y_MC]=ode45(@(t,y)orchard_rhs(t,y,phys_param),tspan,IC);
    
    % recover the solution at the desired instants
    Y_MC = interp1(time_MC,Y_MC,time);
    
    % sample value value for the DoFs
          y1(imc,:) = Y_MC(:,1);
        phi1(imc,:) = Y_MC(:,2);
        phi2(imc,:) = Y_MC(:,3);
    
    % sample value for the DoFs time derivatives
      y1dot(imc,:) = Y_MC(:,4);
    phi1dot(imc,:) = Y_MC(:,5);
    phi2dot(imc,:) = Y_MC(:,6);
    
    % compute MC convergence metric
    normL2_Q(1,imc) = trapz(time,  y1(imc,:).^2 + ...
                                 phi1(imc,:).^2 + ...
                                 phi2(imc,:).^2);

end

% tower horizontal displacement/velocity
x2    = - L1*sin(phi1) ...
        - L2*sin(phi2);
x2dot = - L1*cos(phi1).*phi1dot...
        - L2*cos(phi2).*phi2dot;

% tower vertical displacement
y2   =   y1 ...
       + L1*cos(phi1) ...
       + L2*cos(phi2);
y2dot =   y1dot ...
        - L1*sin(phi1).*phi1dot ...
        - L2*sin(phi2).*phi2dot;

toc
% -----------------------------------------------------------


% compute the statistics
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- computing statistics --- ');
disp(' ');
disp('    ... ');
disp(' ');

% sample average of the random processes
     x2_smp_avg = mean(x2);
     y1_smp_avg = mean(y1);
     y2_smp_avg = mean(y2);
   phi1_smp_avg = mean(phi1);
   phi2_smp_avg = mean(phi2);
   
  x2dot_smp_avg = mean(x2dot);
  y1dot_smp_avg = mean(y1dot);
  y2dot_smp_avg = mean(y2dot);
phi1dot_smp_avg = mean(phi1dot);
phi2dot_smp_avg = mean(phi2dot);

% time average of the steady state random processes
     x2_time_avg = mean(  x2,2);
     y1_time_avg = mean(  y1,2);
     y2_time_avg = mean(  y2,2);
   phi1_time_avg = mean(phi1,2);
   phi2_time_avg = mean(phi2,2);
   
  x2dot_time_avg = mean(  x2dot,2);
  y1dot_time_avg = mean(  y1dot,2);
  y2dot_time_avg = mean(  y2dot,2);
phi1dot_time_avg = mean(phi1dot,2);
phi2dot_time_avg = mean(phi2dot,2);

% standard deviation of the random processes
     x2_std = std(x2);
     y1_std = std(y1);
     y2_std = std(y2);
   phi1_std = std(phi1);
   phi2_std = std(phi2);
   
  x2dot_std = std(x2dot);
  y1dot_std = std(y1dot);
  y2dot_std = std(y2dot);
phi1dot_std = std(phi1dot);
phi2dot_std = std(phi2dot);

% skewness of the random processes
     x2_skew = skewness(x2);
     y1_skew = skewness(y1);
     y2_skew = skewness(y2);
   phi1_skew = skewness(phi1);
   phi2_skew = skewness(phi2);
   
  x2dot_skew = skewness(x2dot);
  y1dot_skew = skewness(y1dot);
  y2dot_skew = skewness(y2dot);
phi1dot_skew = skewness(phi1dot);
phi2dot_skew = skewness(phi2dot);

% kurtosis of the random processes
     x2_kurt = kurtosis(x2)-3;
     y1_kurt = kurtosis(y1)-3;
     y2_kurt = kurtosis(y2)-3;
   phi1_kurt = kurtosis(phi1)-3;
   phi2_kurt = kurtosis(phi2)-3;
   
  x2dot_kurt = kurtosis(x2dot)-3;
  y1dot_kurt = kurtosis(y1dot)-3;
  y2dot_kurt = kurtosis(y2dot)-3;
phi1dot_kurt = kurtosis(phi1dot)-3;
phi2dot_kurt = kurtosis(phi2dot)-3;


% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence bands upper bounds
     x2_upp = prctile(      x2,r_plus);
     y1_upp = prctile(      y1,r_plus);
     y2_upp = prctile(      y2,r_plus);
   phi1_upp = prctile(    phi1,r_plus);
   phi2_upp = prctile(    phi2,r_plus);
   
  x2dot_upp = prctile(  x2dot,r_plus);
  y1dot_upp = prctile(  y1dot,r_plus);
  y2dot_upp = prctile(  y2dot,r_plus);
phi1dot_upp = prctile(phi1dot,r_plus);
phi2dot_upp = prctile(phi2dot,r_plus);

% confidence bands lower bounds
     x2_low = prctile(      x2,r_minus);
     y1_low = prctile(      y1,r_minus);
     y2_low = prctile(      y2,r_minus);
   phi1_low = prctile(    phi1,r_minus);
   phi2_low = prctile(    phi2,r_minus);
   
  x2dot_low = prctile(  x2dot,r_minus);
  y1dot_low = prctile(  y1dot,r_minus);
  y2dot_low = prctile(  y2dot,r_minus);
phi1dot_low = prctile(phi1dot,r_minus);
phi2dot_low = prctile(phi2dot,r_minus);

% compute MC convergence metric
MC_conv = randvar_mc_conv(normL2_Q,Ns);
    
% number of bins
Nbins = round(sqrt(Ns));

% random processes distribution
[  x2_bins_t,  x2_freq_t] = randvar_pdf(randvar_normalize(    x2),Nbins);
[  y1_bins_t,  y1_freq_t] = randvar_pdf(randvar_normalize(    y1),Nbins);
[  y2_bins_t,  y2_freq_t] = randvar_pdf(randvar_normalize(    y2),Nbins);
[phi1_bins_t,phi1_freq_t] = randvar_pdf(randvar_normalize(  phi1),Nbins);
[phi2_bins_t,phi2_freq_t] = randvar_pdf(randvar_normalize(  phi2),Nbins);


% number of KDS points
Nksd = 100;

% random processes KSD distribution
[  x2_ksd_t,  x2_supp_t] = randvar_ksd(randvar_normalize(    x2),Nksd);
[  y1_ksd_t,  y1_supp_t] = randvar_ksd(randvar_normalize(    y1),Nksd);
[  y2_ksd_t,  y2_supp_t] = randvar_ksd(randvar_normalize(    y2),Nksd);
[phi1_ksd_t,phi1_supp_t] = randvar_ksd(randvar_normalize(  phi1),Nksd);
[phi2_ksd_t,phi2_supp_t] = randvar_ksd(randvar_normalize(  phi2),Nksd);


% histograms estimations (steady) PDFs
[  x2_bins,  x2_freq] = randvar_pdf(randvar_normalize(  x2_time_avg),Nbins);
[  y1_bins,  y1_freq] = randvar_pdf(randvar_normalize(  y1_time_avg),Nbins);
[  y2_bins,  y2_freq] = randvar_pdf(randvar_normalize(  y2_time_avg),Nbins);
[phi1_bins,phi1_freq] = randvar_pdf(randvar_normalize(phi1_time_avg),Nbins);
[phi2_bins,phi2_freq] = randvar_pdf(randvar_normalize(phi2_time_avg),Nbins);

% KSD estimations (steady) PDFs
[  x2_ksd,  x2_supp] = ksdensity(randvar_normalize(  x2_time_avg));
[  y1_ksd,  y1_supp] = ksdensity(randvar_normalize(  y1_time_avg));
[  y2_ksd,  y2_supp] = ksdensity(randvar_normalize(  y2_time_avg));
[phi1_ksd,phi1_supp] = ksdensity(randvar_normalize(phi1_time_avg));
[phi2_ksd,phi2_supp] = ksdensity(randvar_normalize(phi2_time_avg));

% random variables entropy
%  Sx2 = shannon_entropy(  x2_bins,  x2_freq);
%  Sy1 = shannon_entropy(  y1_bins,  y1_freq);
%  Sy2 = shannon_entropy(  y2_bins,  y2_freq);
%Sphi1 = shannon_entropy(phi1_bins,phi1_freq);
%Sphi2 = shannon_entropy(phi2_bins,phi2_freq);

% random variables statistics
% x2_time_avg_mean = mean(x2_time_avg)
% x2_time_avg_std  =  std(x2_time_avg)
% x2_time_avg_skew = skewness(x2_time_avg)
% x2_time_avg_kurt = kurtosis(x2_time_avg)-3
% 
% y1_time_avg_mean = mean(y1_time_avg)
% y1_time_avg_std  =  std(y1_time_avg)
% y1_time_avg_skew = skewness(y1_time_avg)
% y1_time_avg_kurt = kurtosis(y1_time_avg)-3
% 
% y2_time_avg_mean = mean(y2_time_avg)
% y2_time_avg_std  =  std(y2_time_avg)
% y2_time_avg_skew = skewness(y2_time_avg)
% y2_time_avg_kurt = kurtosis(y2_time_avg)-3
% 
% phi1_time_avg_mean = mean(phi1_time_avg)
% phi1_time_avg_std  =  std(phi1_time_avg)
% phi1_time_avg_skew = skewness(phi1_time_avg)
% phi1_time_avg_kurt = kurtosis(phi1_time_avg)-3
% 
% phi2_time_avg_mean = mean(phi2_time_avg)
% phi2_time_avg_std  =  std(phi2_time_avg)
% phi2_time_avg_skew = skewness(phi2_time_avg)
% phi2_time_avg_kurt = kurtosis(phi2_time_avg)-3


% probability of random variables be <= given value
Px2   = (0.1*B2-x2_smp_avg)./x2_std;
Py1   = 0.0*ones(1,Ndt);
Py2   = 0.0*ones(1,Ndt);
Pphi1 = 0.0*ones(1,Ndt);
Pphi2 = 0.0*ones(1,Ndt);

prob_x2   = randvar_probval(  x2_bins_t,  x2_freq_t,Px2);
prob_y1   = randvar_probval(  y1_bins_t,  y1_freq_t,Py1);
prob_y2   = randvar_probval(  y2_bins_t,  y2_freq_t,Py2);
prob_phi1 = randvar_probval(phi1_bins_t,phi1_freq_t,Pphi1);
prob_phi2 = randvar_probval(phi2_bins_t,phi2_freq_t,Pphi2);

toc
% -----------------------------------------------------------




% save simulation data into a file
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- Saving Workspace --- ');
disp(' ');
disp('    ... ');
disp(' ');


save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------


