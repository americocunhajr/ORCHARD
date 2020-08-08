
% ----------------------------------------------------------------- 
%  orchard_rhs.m
%
%  This function defines the following nonlinear system of 
%  ordinary differential equations
%
%    [M]Qacce + [C]Qvelo + [N]Qvelo.^2 + [K]Qdisp + g = h,
%  
%  where
%  
%   [M]   - configuration dependent mass matrix
%   [C]   - configuration dependent damping matrix
%   [N]   - configuration dependent circulatory matrix
%   [K]   - configuration dependent stiffness matrix
%   g     - configuration dependent gravity vector
%   h     - configuration dependent external excitation
%   Qdisp - generalized displacement vector
%   Qvelo - generalized velocity vector
%   Qacee - generalized acceleration vector
%
%  
%  This dynamical system is implemented in state space form
%  
%    dq/dt = f(q,t),
%  
%  where
%  
%    q = (y1 phi1 phi2 y1_dot phi1_dot phi2_dot)^T.
%
%  Reference:
%  S. Sartori Junior, J. M. Balthazar and B. R. Pontes Junior
%  Non-linear dynamics of a tower orchard sprayer based on an
%  inverted pendulum model
%  Biosystems Engineering
%  vol. 103 p. 417-426, 2009
%  
%  
%  Input:
%  t          - time (s)
%  q          - state space vector (6 x 1)
%  phys_param - physical parameters vector
%  
%  
%  
%  
%  
%  Output:
%  qdot - right hand side function (6 x 1)
% ----------------------------------------------------------------- 
%  programmers: Americo Barbosa da Cunha Junior
%               americo.cunhajr@gmail.com
%
%               Jorge Luis Palacios Felix
%               jorge.felix@uffs.edu.br
%
%  last update: Feb 15, 2017
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function qdot = orchard_rhs(t,q,phys_param)

% physical parameters
kt  = phys_param(1);
ct  = phys_param(2);
k1  = phys_param(3);
k2  = phys_param(4);
c1  = phys_param(5);
c2  = phys_param(6);
B1  = phys_param(7);
B2  = phys_param(8);
g   = phys_param(9);
m1  = phys_param(10);
m2  = phys_param(11);
L1  = phys_param(12);
L2  = phys_param(13);
I1  = phys_param(14);
I2  = phys_param(15);
A   = phys_param(16);
w   = phys_param(17);
rho = phys_param(18);


% define state variables
  y1     = q(1);
phi1     = q(2);
phi2     = q(3);
  y1_dot = q(4);
phi1_dot = q(5);
phi2_dot = q(6);


% preallocate memory for state vector time derivative
qdot = zeros(6,1);

% preallocate memory for mass matrix
M = zeros(3,3);

% preallocate memory for damping matrix
C = zeros(3,3);

% preallocate memory for circulatory matrix
N = zeros(3,3);

% preallocate memory for stiffness matrix
K = zeros(3,3);

% preallocate memory for gravity vector
G = zeros(3,1);

% preallocate memory for external excitation vector
H = zeros(3,1);


% external excitation
ye1 = A*sin(w*t);
ye2 = A*sin(w*t + rho);

% external excitation time derivative
ye1_dot = A*w*cos(w*t);
ye2_dot = A*w*cos(w*t + rho);


% auxiliar variables 
% (to save floating point operations)
sin_phi1 = sin(phi1);
cos_phi1 = cos(phi1);
sin_phi2 = sin(phi2);
cos_phi2 = cos(phi2);

sin_phi2_minus_phi1 = sin(phi2-phi1);
cos_phi2_minus_phi1 = cos(phi2-phi1);


% mass matrix entries
M(1,1) = m1 + m2;
M(1,2) = -m2*L1*sin_phi1;
M(1,3) = -m2*L2*sin_phi2;
M(2,1) = M(1,2);
M(2,2) = I1 + m2*L1*L1;
M(2,3) = m2*L1*L2*cos_phi2_minus_phi1;
M(3,1) = M(1,3);
M(3,2) = M(2,3);
M(3,3) = I2 + m2*L2*L2;

% damping matrix entries
C(1,1) = c1 + c2;
C(1,2) = (c2*B2-c1*B1)*cos_phi1;
C(2,1) = C(1,2);
C(2,2) = ct + (c1*B1*B1+c2*B2*B2)*cos_phi1*cos_phi1;
C(2,3) = -ct;
C(3,2) = C(2,3);
C(3,3) = ct;

% circulatory matrix entries
N(1,2) = -m2*L1*cos_phi1;
N(1,3) = -m2*L2*cos_phi2;
N(2,3) = -m2*L1*L2*sin_phi2_minus_phi1;
N(3,2) = N(2,3);

% stiffness matrix entries
K(1,1) = k1 + k2;
K(2,1) = (k2*B2-k1*B1)*cos_phi1;
K(2,2) =  kt;
K(2,3) = -kt;
K(3,2) = K(2,3);
K(3,3) =  kt;

% gravity vector entries
G(1,1) = (k2*B2-k1*B1)*sin_phi1 + (m1+m2)*g;
G(2,1) = ((k1*B1*B1+k2*B2*B2)*cos_phi1 - m2*g*L1)*sin_phi1;
G(3,1) = -m2*g*L2*sin_phi2;

% external excitation vector entries
H(1,1) = k1*ye1 + ...
         k2*ye2 + ...
         c1*ye1_dot + ...
         c2*ye2_dot;

H(2,1) = (B2*(k2*ye2 + c2*ye2_dot) - ...
          B1*(k1*ye1 + c1*ye1_dot))*cos_phi1;


% define state space system of equations
qdot(1:3) = [y1_dot; phi1_dot; phi2_dot];
qdot(4:6) = -M\(N*[y1_dot; phi1_dot; phi2_dot].^2 + ...
                C*[y1_dot; phi1_dot; phi2_dot] + ...
                K*[y1; phi1; phi2] - H + G);

end
% -----------------------------------------------------------------
