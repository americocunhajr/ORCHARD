
% -----------------------------------------------------------------
%  main_orchard_randfield.m
%
%  This is the main file for a program that geneates a real-valued
%  random field {S(x), x \in D}. The representation of this random
%  process is performed through the Karhunen-Lo?ve decomposition.
%  
%  References:
%  R. Ghanem and P. Spanos
%  Stochastic Finite Element Method: A Spectral Approch
%  Dover Publications, 2003
%  Pages 28-33
%  
%  D. Xiu
%  Numerical Methods for Stochastic Computations
%  Princeton University Press, 2009
%  Pages 47-49
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Mar 31, 2016
% -----------------------------------------------------------------


clear all
close all
clc


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Random Field Generation                            ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'random_field';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% normalized domain upper limit (m)
b = 1.0;

% correlation length (m)
a = 0.3;

% random field std. deviation (m)
sig = 0.025;

% number of mesh points for domain discretization
Nx = 2000;

% number of eigenpairs to be computed
Neig = 200;

% number of random samples
Ns = 2^10;

% number of eigenfunctions to display
Ndisp = 4;
% -----------------------------------------------------------



% autocorrelation function eigenpairs
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- computing the autocorr function eigenpairs --- ');
disp(' ');
disp('    ... ');
disp(' ');

[lambda,phi,omega,xmesh,iroot,ising,iter] = ...
                  fredholm_expcorr_eig(b,a,Neig,Nx);
    
toc
% -----------------------------------------------------------


% representaton of random field S
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- computing the random field representation --- ');
disp(' ');
disp('    ... ');
disp(' ');

% domain limits D = [x_min,x_max]
x_min = 0.0;
x_max = 100.0;

% new domain mesh
xmesh_new = linspace(x_min,x_max,Nx);

% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setDefaultStream(rng_stream);

% generate the uncorrelated random variables
ksi = randn(Ns,Neig);
%ksi = -1 + rand(Ns,Neig);

% multiply eigenvalues by variance
lambda = sig*sig*lambda;

% representation energy level
energy_level = 0.99;

% number of eigenpairs used in KL expansion
[~,Nkl] = max(cumsum(lambda)/sum(lambda) >= energy_level);

% KL representation of the random field
S = phi(:,1:Nkl)*sqrt(diag(lambda(1:Nkl)))*ksi(:,1:Nkl)';
    
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


figure(1)
eig_index = 1:1:Neig;
fig1 = semilogy(eig_index,lambda,'o--');
xlabel('eigenvalue index')
ylabel('egenvalue')


figure(2)
fig2 = plot(xmesh_new,phi(:,1:Ndisp));
xlabel('domain parameter')
ylabel('eigenvector amplitude')


figure(3)
fig3 = plot(xmesh_new,S(:,1:Ndisp));
xlabel('domain parameter')
ylabel('random field realizations')
% -----------------------------------------------------------
