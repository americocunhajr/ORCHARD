
% -----------------------------------------------------------------
%  post_orchard_ivp.m
%
%  This script performs the post processing of simulation data
%  that comes from orchard nonlinear dynamics.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 15, 2017
% -----------------------------------------------------------------



% plot trailer vertical displacement time series
% ...........................................................
gtitle = '';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_y1_disp'];
flag   = 'eps';
fig1a  = graph_type1(time,Qy1,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1a);
% ...........................................................


% plot trailer angular displacement time series
% ...........................................................
gtitle = '';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_phi1_disp'];
flag   = 'eps';
fig1b  = graph_type1(time,Qphi1,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1b);
% ...........................................................


% plot tower angular displacement time series
% ...........................................................
gtitle = '';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_phi2_disp'];
flag   = 'eps';
fig1c  = graph_type1(time,Qphi2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1c);
% ...........................................................



% plot tower horizontal displacement time series
% ...........................................................
gtitle = '';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = -0.4;
ymax   =  0.4;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_x2_disp'];
flag   = 'eps';
fig1d  = graph_type1(time,Qx2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1d);
% ...........................................................


% plot tower horizontal displacement time series
% ...........................................................
gtitle = '';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_y2_disp'];
flag   = 'eps';
fig1e  = graph_type1(time,Qy2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1e);
% ...........................................................


% plot trailer vertical phase space trajectory
% ...........................................................
gtitle = '';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_y1'];
flag   = 'eps';
fig2a  = graph_type1(Qy1(Nss:Ndt),Qy1dot(Nss:Ndt),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2a);
% ...........................................................


% plot trailer angular phase space trajectory
% ...........................................................
gtitle = '';
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_phi1'];
flag   = 'eps';
fig2b  = graph_type1(Qphi1(Nss:Ndt),Qphi1dot(Nss:Ndt),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2b);
% ...........................................................


% plot tower angular phase space trajectory
% ...........................................................
gtitle = '';
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_phi2'];
flag   = 'eps';
fig2c  = graph_type1(Qphi2(Nss:Ndt),Qphi2dot(Nss:Ndt),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2c);
% ...........................................................


% plot tower horizontal phase space trajectory
% ...........................................................
gtitle = '';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_x2'];
flag   = 'eps';
fig2d  = graph_type1(Qx2(Nss:Ndt),Qx2dot(Nss:Ndt),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2d);
% ...........................................................


% plot tower vertical phase space trajectory
% ...........................................................
gtitle = '';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_y2'];
flag   = 'eps';
fig2e  = graph_type1(Qy2(Nss:Ndt),Qy2dot(Nss:Ndt),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2e);
% ...........................................................


% plot trailer vertical atractor
% ...........................................................
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
zlab   = ' time (s)';
gtitle = '';
xmin   = min(Qy1);
xmax   = max(Qy1);
ymin   = min(Qy1dot);
ymax   = max(Qy1dot);
zmin   = min(time(Nss:Ndt));
zmax   = max(time(Nss:Ndt));
gname  = [num2str(case_name),'__atractor_y1'];
flag   = 'eps';
fig3a  = graph_3d_atractor(   Qy1(Nss:Ndt),...
                          Qy1dot(Nss:Ndt),...
                          time(Nss:Ndt),...
                          time(Nss:Ndt),...
                          gtitle,xlab,ylab,zlab,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3a);
% ...........................................................


% plot trailer angular atractor
% ...........................................................
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
zlab   = ' time (s)';
gtitle = '';
xmin   = min(Qphi1);
xmax   = max(Qphi1);
ymin   = min(Qphi1dot);
ymax   = max(Qphi1dot);
zmin   = min(time(Nss:Ndt));
zmax   = max(time(Nss:Ndt));
gname  = [num2str(case_name),'__atractor_phi1'];
flag   = 'eps';
fig3b  = graph_3d_atractor(   Qphi1(Nss:Ndt),...
                          Qphi1dot(Nss:Ndt),...
                          time(Nss:Ndt),...
                          time(Nss:Ndt),...
                          gtitle,xlab,ylab,zlab,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3b);
% ...........................................................


% plot tower angular atractor
% ...........................................................
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
zlab   = ' time (s)';
gtitle = '';
xmin   = min(Qphi2);
xmax   = max(Qphi2);
ymin   = min(Qphi2dot);
ymax   = max(Qphi2dot);
zmin   = min(time(Nss:Ndt));
zmax   = max(time(Nss:Ndt));
gname  = [num2str(case_name),'__atractor_phi2'];
flag   = 'eps';
fig3c  = graph_3d_atractor(   Qphi2(Nss:Ndt),...
                           Qphi2dot(Nss:Ndt),...
                           time(Nss:Ndt),...
                           time(Nss:Ndt),...
                           gtitle,xlab,ylab,zlab,...
                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3c);
% ...........................................................


% plot tower horizontal atractor
% ...........................................................
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
zlab   = ' time (s)';
gtitle = '';
xmin   = -0.5;
xmax   = 0.5;
ymin   = -2;
ymax   = 2;
%xmin   = min(Qx2);
%xmax   = max(Qx2);
%ymin   = min(Qx2dot);
%ymax   = max(Qx2dot);
zmin   = min(time(Nss:Ndt));
zmax   = max(time(Nss:Ndt));
gname  = [num2str(case_name),'__atractor_x2'];
flag   = 'eps';
fig3d  = graph_3d_atractor(   Qx2(Nss:Ndt),...
                           Qx2dot(Nss:Ndt),...
                           time(Nss:Ndt),...
                           time(Nss:Ndt),...
                           gtitle,xlab,ylab,zlab,...
                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3d);
% ...........................................................


% plot tower vertical atractor
% ...........................................................
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
zlab   = ' time (s)';
gtitle = '';
xmin   = min(Qy2);
xmax   = max(Qy2);
ymin   = min(Qy2dot);
ymax   = max(Qy2dot);
zmin   = min(time(Nss:Ndt));
zmax   = max(time(Nss:Ndt));
gname  = [num2str(case_name),'__atractor_y2'];
flag   = 'eps';
fig3e  = graph_3d_atractor(   Qy2(Nss:Ndt),...
                           Qy2dot(Nss:Ndt),...
                           time(Nss:Ndt),...
                           time(Nss:Ndt),...
                           gtitle,xlab,ylab,zlab,...
                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3e);
% ...........................................................
