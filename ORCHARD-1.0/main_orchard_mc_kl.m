
% -----------------------------------------------------------------
%  main_orchard_mc_kl.m
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
disp(' jorge.felix@uffs.edu.br                            ')
disp('                                                    ')
disp(' Jose Manoel Balthazar                              ')
disp(' jmbaltha@gmail.com                                 ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------


% simulation information
% -----------------------------------------------------------
%case_name = 'orchad_mc_kl_test';

case_name = 'orchad_mc_kl_corr_1m_v_12kmph';
%case_name = 'orchad_mc_kl_corr_10m_v_12kmph';

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

% velocity of translation (km/h)
vtrans_kmph = 12.0;

% velocity of translation (m/s)
vtrans = vtrans_kmph*(1000/3600);

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


% frequency band (Hz)
freq_min = 0.0;
freq_max = 5.0;

% Nyquist sampling frequency (Hz)
freq_samp = 2*freq_max;

% time step or sampling period (s)
dt = 1/freq_samp;

% number of samples for signal copy
Nsamp = 300;
%Nsamp = 1200;
%Nsamp = 3000;
%Nsamp = 6000;
%Nsamp = 7500;

% period of a signal copy (s)
Tsignal = Nsamp*dt;

% initial time of analysis (s)
t0 = 0.0;

% final time of analysis (s)
t1 = t0 + Tsignal;

% number of time steps
Ndt = Nsamp;

% number of time steps for steady state
Nss = round(0.5*Ndt);
% -----------------------------------------------------------



% generating external forcing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- generating random parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setDefaultStream(rng_stream); % Matlab 2009
%RandStream.setGlobalStream(rng_stream); % Matlab 2013

% number of random samples
%Ns = 4;
%Ns = 16;
%Ns = 64;
Ns = 256;
%Ns = 1024;
%Ns = 4096;
%Ns = 16384;
%Ns = 65536;

% normalized translation domain upper limit (m)
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
% (divide by vtrans to convert from meter to second)
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



% Monte Carlo simulation
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Monte Carlo Simulation --- ');
disp(' ');


% preallocate memory for DoFs
MC_y1   = zeros(Ns,Ndt);
MC_phi1 = zeros(Ns,Ndt);
MC_phi2 = zeros(Ns,Ndt);

% preallocate memory for DoFs time derivatives
MC_y1dot   = zeros(Ns,Ndt);
MC_phi1dot = zeros(Ns,Ndt);
MC_phi2dot = zeros(Ns,Ndt);

% preallocate memory for MC conv metric
MC_conv = zeros(1,Ns);

% physical parameters vector
phys_param = [kt ct k1 k2 c1 c2 B1 B2 g m1 m2 L1 L2 I1 I2];

% initial conditions
IC = [y1_0; phi1_0; phi2_0; y1dot_0; phi1dot_0; phi2dot_0];

% instants of analysis (s)
time = linspace(t0,t1,Ndt);
%tspan = [t0 t1];

% optional parameter for ODE solver
% (RelTol = 1.0e-3 AbsTol = 1.0e-6)
opt = odeset ('RelTol',1.0e-3,'AbsTol',1.0e-6);


for imc=1:Ns
    
    if mod(imc,16) == 0
        disp('')
        disp(imc)
    end
    
    % ODE solver Runge-Kutta45
    [MC_time,MC_Y] = ode45(@(t,y)orchard_rhs_kl(t,y,phys_param,...
                                                  Ye1(:,imc),...
                                                  Ye2(:,imc),...
                                                  Ye1dot(:,imc),...
                                                  Ye2dot(:,imc),...
                                                  time),...
                                                [t0 t1],IC);
    
    % recover the solution at the desired instants
    MC_Y = interp1(MC_time,MC_Y,time);
    
    % time series for generalized displacements
    MC_y1(imc,:)   = MC_Y(:,1);
    MC_phi1(imc,:) = MC_Y(:,2);
    MC_phi2(imc,:) = MC_Y(:,3);
    
    % time series for generalized velicities
    MC_y1dot(imc,:)   = MC_Y(:,4);
    MC_phi1dot(imc,:) = MC_Y(:,5);
    MC_phi2dot(imc,:) = MC_Y(:,6);
    
    % compute MC convergence metric
    MC_conv(1,imc) = trapz(time,  MC_y1(imc,:).^2 + ...
                                MC_phi1(imc,:).^2 + ...
                                MC_phi2(imc,:).^2);

end

% compute MC convergence metric
MC_conv = randvar_mc_conv(MC_conv,Ns);

% tower horizontal displacement/velocity
MC_x2    = - L1*sin(MC_phi1) ...
           - L2*sin(MC_phi2);
MC_x2dot = - L1*cos(MC_phi1).*MC_phi1dot...
           - L2*cos(MC_phi2).*MC_phi2dot;

% tower vertical displacement
MC_y2    = MC_y1 ...
           + L1*cos(MC_phi1) ...
           + L2*cos(MC_phi2);
MC_y2dot = MC_y1dot ...
           - L1*sin(MC_phi1).*MC_phi1dot ...
           - L2*sin(MC_phi2).*MC_phi2dot;

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
     x2_smp_avg = mean(MC_x2  );
     y1_smp_avg = mean(MC_y1  );
     y2_smp_avg = mean(MC_y2  );
   phi1_smp_avg = mean(MC_phi1);
   phi2_smp_avg = mean(MC_phi2);
   
  x2dot_smp_avg = mean(MC_x2dot  );
  y1dot_smp_avg = mean(MC_y1dot  );
  y2dot_smp_avg = mean(MC_y2dot  );
phi1dot_smp_avg = mean(MC_phi1dot);
phi2dot_smp_avg = mean(MC_phi2dot);

% time average of the steady state random processes
     x2_time_avg = mean(MC_x2  ,2);
     y1_time_avg = mean(MC_y1  ,2);
     y2_time_avg = mean(MC_y2  ,2);
   phi1_time_avg = mean(MC_phi1,2);
   phi2_time_avg = mean(MC_phi2,2);
   
  x2dot_time_avg = mean(MC_x2dot  ,2);
  y1dot_time_avg = mean(MC_y1dot  ,2);
  y2dot_time_avg = mean(MC_y2dot  ,2);
phi1dot_time_avg = mean(MC_phi1dot,2);
phi2dot_time_avg = mean(MC_phi2dot,2);

% standard deviation of the random processes
     x2_std = std(MC_x2  );
     y1_std = std(MC_y1  );
     y2_std = std(MC_y2  );
   phi1_std = std(MC_phi1);
   phi2_std = std(MC_phi2);
   
  x2dot_std = std(MC_x2dot  );
  y1dot_std = std(MC_y1dot  );
  y2dot_std = std(MC_y2dot  );
phi1dot_std = std(MC_phi1dot);
phi2dot_std = std(MC_phi2dot);

% skewness of the random processes
     x2_skew = skewness(MC_x2  );
     y1_skew = skewness(MC_y1  );
     y2_skew = skewness(MC_y2  );
   phi1_skew = skewness(MC_phi1);
   phi2_skew = skewness(MC_phi2);
   
  x2dot_skew = skewness(MC_x2dot  );
  y1dot_skew = skewness(MC_y1dot  );
  y2dot_skew = skewness(MC_y2dot  );
phi1dot_skew = skewness(MC_phi1dot);
phi2dot_skew = skewness(MC_phi2dot);

% kurtosis of the random processes
     x2_kurt = kurtosis(MC_x2  )-3;
     y1_kurt = kurtosis(MC_y1  )-3;
     y2_kurt = kurtosis(MC_y2  )-3;
   phi1_kurt = kurtosis(MC_phi1)-3;
   phi2_kurt = kurtosis(MC_phi2)-3;
   
  x2dot_kurt = kurtosis(MC_x2dot  )-3;
  y1dot_kurt = kurtosis(MC_y1dot  )-3;
  y2dot_kurt = kurtosis(MC_y2dot  )-3;
phi1dot_kurt = kurtosis(MC_phi1dot)-3;
phi2dot_kurt = kurtosis(MC_phi2dot)-3;


% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence bands upper bounds
     x2_upp = prctile(MC_x2  ,r_plus);
     y1_upp = prctile(MC_y1  ,r_plus);
     y2_upp = prctile(MC_y2  ,r_plus);
   phi1_upp = prctile(MC_phi1,r_plus);
   phi2_upp = prctile(MC_phi2,r_plus);
   
  x2dot_upp = prctile(MC_x2dot  ,r_plus);
  y1dot_upp = prctile(MC_y1dot  ,r_plus);
  y2dot_upp = prctile(MC_y2dot  ,r_plus);
phi1dot_upp = prctile(MC_phi1dot,r_plus);
phi2dot_upp = prctile(MC_phi2dot,r_plus);

% confidence bands lower bounds
     x2_low = prctile(MC_x2  ,r_minus);
     y1_low = prctile(MC_y1  ,r_minus);
     y2_low = prctile(MC_y2  ,r_minus);
   phi1_low = prctile(MC_phi1,r_minus);
   phi2_low = prctile(MC_phi2,r_minus);
   
  x2dot_low = prctile(MC_x2dot  ,r_minus);
  y1dot_low = prctile(MC_y1dot  ,r_minus);
  y2dot_low = prctile(MC_y2dot  ,r_minus);
phi1dot_low = prctile(MC_phi1dot,r_minus);
phi2dot_low = prctile(MC_phi2dot,r_minus);
    
% number of bins
Nbins = round(sqrt(Ns));

% number of KSD points
Nksd = 100;

% time-average histograms estimation
[  x2_bins,  x2_freq] = randvar_pdf(randvar_normalize(  x2_time_avg),Nbins);
[  y1_bins,  y1_freq] = randvar_pdf(randvar_normalize(  y1_time_avg),Nbins);
[  y2_bins,  y2_freq] = randvar_pdf(randvar_normalize(  y2_time_avg),Nbins);
[phi1_bins,phi1_freq] = randvar_pdf(randvar_normalize(phi1_time_avg),Nbins);
[phi2_bins,phi2_freq] = randvar_pdf(randvar_normalize(phi2_time_avg),Nbins);

% time-average PDFs curves estimation
[  x2_ksd,  x2_supp] = randvar_ksd(randvar_normalize(  x2_time_avg),Nksd);
[  y1_ksd,  y1_supp] = randvar_ksd(randvar_normalize(  y1_time_avg),Nksd);
[  y2_ksd,  y2_supp] = randvar_ksd(randvar_normalize(  y2_time_avg),Nksd);
[phi1_ksd,phi1_supp] = randvar_ksd(randvar_normalize(phi1_time_avg),Nksd);
[phi2_ksd,phi2_supp] = randvar_ksd(randvar_normalize(phi2_time_avg),Nksd);

% random processes histograms estimation
[  x2_bins_t,  x2_freq_t] = randvar_pdf(randvar_normalize(MC_x2  ),Nbins);
[  y1_bins_t,  y1_freq_t] = randvar_pdf(randvar_normalize(MC_y1  ),Nbins);
[  y2_bins_t,  y2_freq_t] = randvar_pdf(randvar_normalize(MC_y2  ),Nbins);
[phi1_bins_t,phi1_freq_t] = randvar_pdf(randvar_normalize(MC_phi1),Nbins);
[phi2_bins_t,phi2_freq_t] = randvar_pdf(randvar_normalize(MC_phi2),Nbins);

% random processes marginal PDFs curves estimation
[  x2_ksd_t,  x2_supp_t] = randvar_ksd(randvar_normalize(MC_x2  ),Nksd);
[  y1_ksd_t,  y1_supp_t] = randvar_ksd(randvar_normalize(MC_y1  ),Nksd);
[  y2_ksd_t,  y2_supp_t] = randvar_ksd(randvar_normalize(MC_y2  ),Nksd);
[phi1_ksd_t,phi1_supp_t] = randvar_ksd(randvar_normalize(MC_phi1),Nksd);
[phi2_ksd_t,phi2_supp_t] = randvar_ksd(randvar_normalize(MC_phi2),Nksd);


% random variables entropy
%  Sx2 = shannon_entropy(  x2_bins,  x2_freq);
%  Sy1 = shannon_entropy(  y1_bins,  y1_freq);
%  Sy2 = shannon_entropy(  y2_bins,  y2_freq);
%Sphi1 = shannon_entropy(phi1_bins,phi1_freq);
%Sphi2 = shannon_entropy(phi2_bins,phi2_freq);


% random variables statistics
x2_time_avg_mean =     mean(x2_time_avg)
x2_time_avg_std  =      std(x2_time_avg)
x2_time_avg_skew = skewness(x2_time_avg)
x2_time_avg_kurt = kurtosis(x2_time_avg)-3

y1_time_avg_mean =     mean(y1_time_avg)
y1_time_avg_std  =      std(y1_time_avg)
y1_time_avg_skew = skewness(y1_time_avg)
y1_time_avg_kurt = kurtosis(y1_time_avg)-3

y2_time_avg_mean =     mean(y2_time_avg)
y2_time_avg_std  =      std(y2_time_avg)
y2_time_avg_skew = skewness(y2_time_avg)
y2_time_avg_kurt = kurtosis(y2_time_avg)-3

phi1_time_avg_mean =     mean(phi1_time_avg)
phi1_time_avg_std  =      std(phi1_time_avg)
phi1_time_avg_skew = skewness(phi1_time_avg)
phi1_time_avg_kurt = kurtosis(phi1_time_avg)-3

phi2_time_avg_mean =     mean(phi2_time_avg)
phi2_time_avg_std  =      std(phi2_time_avg)
phi2_time_avg_skew = skewness(phi2_time_avg)
phi2_time_avg_kurt = kurtosis(phi2_time_avg)-3

% probability of |random variables(t)| <= threshold
x2_plus   = ( 0.3*B2-mean(MC_x2))./std(MC_x2);
x2_minus  = (-0.3*B2-mean(MC_x2))./std(MC_x2);

P_x2_plus   = randvar_probval(x2_bins_t,x2_freq_t,x2_plus );
P_x2_minus  = randvar_probval(x2_bins_t,x2_freq_t,x2_minus);
P_x2        = P_x2_plus - P_x2_minus;

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




% post processing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


post_orchard_mc_kl;


% animate histogram
% ...........................................................
vtitle = 'horizontal displacement distribution';
xlab   = ' (normalized) displacement ';
ylab   = ' probability density function';
xmin   = -4;
xmax   =  4;
ymin   =  0;
ymax   =  1;
vname  = [num2str(case_name),'__pdf_video_x2'];
%fig11  = post_histogram_evolution(time,x2_bins_t,x2_freq_t,...
%                                  vtitle,xlab,ylab,...
%                                  xmin,xmax,ymin,ymax,vname);
% ...........................................................
                                    
toc
% -----------------------------------------------------------


